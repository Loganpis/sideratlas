"""Conservative radial-sector leaky diffusion post-processing.

Gas trajectories never change. Finite cooling photons and intercepted stellar
photons enter a radiation reservoir; diffusion escape and approximate radiation
expansion work are recorded. No new recombination, heating or hydro solution.
"""
import numpy as np
from scipy.spatial import cKDTree,ConvexHull
from numba import njit
from src.detail.fields import coarse_sigma
from src.stellar_structure import MSUN,RSUN,G
SIGMA=5.670374419e-8;C=299792458.;RM=RSUN/100;MM=MSUN/1000;CV=1.5*1.380649e-23/(1.3*1.67262192369e-27)

def sphere(n=1536):
    i=np.arange(n);z=1-2*(i+.5)/n;a=i*2.399963229728653;v=np.column_stack([np.sqrt(1-z*z)*np.cos(a),np.sqrt(1-z*z)*np.sin(a),z]);faces=ConvexHull(v).simplices
    a,b,c=v[faces[:,0]],v[faces[:,1]],v[faces[:,2]];omega=2*np.arctan2(abs(np.sum(a*np.cross(b,c),axis=1)),1+np.sum(a*b+b*c+c*a,axis=1));weight=np.zeros(n)
    for j in range(3):np.add.at(weight,faces[:,j],omega/3)
    return v,faces,weight

@njit(cache=True)
def radial_scatter(r,sigma,near,aw,mass,T,centers,edges,n):
    m=np.zeros((n,len(centers)));heat=m.copy()
    for j in range(len(r)):
        h=max(.025,sigma[j]);norm=0.
        for k in range(len(centers)):
            d=(centers[k]-r[j])/h
            if abs(d)<4:norm+=np.exp(-.5*d*d)*(edges[k+1]-edges[k])
        if norm==0:continue
        for k in range(len(centers)):
            d=(centers[k]-r[j])/h
            if abs(d)>4:continue
            w=mass[j]*np.exp(-.5*d*d)*(edges[k+1]-edges[k])/norm
            for l in range(near.shape[1]):
                value=w*aw[j,l];m[near[j,l],k]+=value;heat[near[j,l],k]+=value*T[j]
    return m,heat

def source_state(m,t,q):
    s=m.scene;report=m.report;dt=report['physical_time_unit_seconds'];energy=G*MSUN**2/RSUN/1e7;lum=float(np.interp(t,m.times,m.data['gas_luminosity_code']+m.data['remnant_luminosity_code']))*energy/dt
    heat=float(np.interp(t,m.times,m.data['heat_code']))*energy;R=q['radius'];base=report['donor']['Teff']*(s.R/R)**.5;a,dis,alive,_=s.state(t)
    v=np.clip((1.05*R-a)/(.65*R),0,1);engulf=v*v*(3-2*v);cv=1.5*1.380649e-23/(.62*1.67262192369e-27)
    response=base+.06*heat/(cv*(s.M-.179547731)*MM)*engulf;area=4*np.pi*(R*RM)**2;power=min(.5*lum*engulf,area*SIGMA*max(0,response**4-base**4));Teff=(base**4+power/(area*SIGMA))**.25
    return area*SIGMA*Teff**4,Teff,max(0,lum-power)

def columns(m,t,dirs,omega,edges,tree):
    q=m.at(float(t));ids=np.flatnonzero(q['live']);rr=(edges[:-1]+edges[1:])/2;nr=len(rr);n=len(dirs);dr=np.diff(edges)
    if len(ids):
        x=q['x'][ids]-q['donor'];r=np.linalg.norm(x,axis=1);unit=x/np.maximum(r[:,None],1e-12);sigma=coarse_sigma(x,q['radius'],np.zeros(3));dist,near=tree.query(unit,k=min(n,max(24,int(np.ceil(n*.10)))),workers=2);angle=2*np.arcsin(np.minimum(1,dist/2));sa=np.maximum(sigma/np.maximum(r,1e-9),.025)
        z=-.5*(angle/sa[:,None])**2;z-=z.max(axis=1,keepdims=True);w=np.exp(z)*omega[near];w[angle>4*sa[:,None]]=0;empty=w.sum(axis=1)<=0;w[empty,0]=1;w/=w.sum(axis=1,keepdims=True)
        mass,heat=radial_scatter(r,sigma,near,w,np.full(len(ids),m.mass),q['T'][ids],rr,edges,n)
    else:mass=np.zeros((n,nr));heat=mass.copy()
    vol=omega[:,None]*(edges[1:]**3-edges[:-1]**3)[None,:]/3;dtau=.01*MM/RM**2*mass/vol*dr;outside=np.cumsum(dtau[:,::-1],axis=1)[:,::-1]-.5*dtau;tau0=dtau.sum(axis=1)
    rph=np.full(n,q['radius']);opaque=tau0>2/3
    for i in np.flatnonzero(opaque):
        cumulative=0.
        for k in range(nr-1,-1,-1):
            if cumulative+dtau[i,k]>=2/3:
                rph[i]=max(q['radius'],edges[k+1]-(2/3-cumulative)/max(dtau[i,k],1e-40)*dr[k]);break
            cumulative+=dtau[i,k]
    gasQ=heat*MM*CV/(m.report['config']['gas_cooling_time_code']*m.report['physical_time_unit_seconds']);trapped=gasQ*(-np.expm1(-np.minimum(outside,100)));free=(gasQ-trapped).sum(axis=1)
    Lstar,Tstar,Lavailable=source_state(m,t,q);star_trap=Lstar*omega/(4*np.pi)*(-np.expm1(-np.minimum(tau0,100)));Q=trapped.sum(axis=1)+star_trap
    # Diffusion-time estimate integrated over outward optical depth. This
    # approximation is finite and tends to a light crossing time when thin.
    delay=3*np.cumsum((outside*dr*RM/C)[:,::-1],axis=1)[:,::-1]
    td=(np.sum(trapped*delay,axis=1)+star_trap*delay[:,0])/np.maximum(Q,1e-100);td=np.maximum(rph*RM/C,td)
    return dict(q=q,mass=mass,heat=heat,dtau=dtau,outside=outside,rph=rph,Q=Q,free=free,tdiff=td,Lstar=Lstar,Tstar=Tstar,available=Lavailable,gas_cooling=float(gasQ.sum()),mass_error=float(mass.sum()-len(ids)*m.mass),direct_stellar_power=float(Lstar-star_trap.sum()))

@njit(cache=True)
def render_photosphere(vertices,faces,power,colors,right,up,eye,center,W,H,pixel,cy):
    rgb=np.zeros((H,W,3),np.float32);depth=np.full((H,W),-np.inf,np.float32);mask=np.zeros((H,W),np.float32)
    area=np.zeros(len(vertices));normals=np.zeros_like(vertices)
    for f in faces:
        a,b,c=vertices[f[0]],vertices[f[1]],vertices[f[2]];v=np.cross(b-a,c-a);ar=.5*np.sqrt((v*v).sum())
        if np.dot(v,a-center)<0:v=-v
        for j in f:area[j]+=ar/3;normals[j]+=v
    for i in range(len(vertices)):normals[i]/=max(1e-30,np.sqrt((normals[i]*normals[i]).sum()))
    # Here power/area is the emergent bolometric flux. A normalized linear
    # limb law is an imaging/atmosphere closure, not extra stellar energy.
    flux=power/np.maximum(area,1e-30)
    projected=np.empty_like(vertices)
    for i in range(len(vertices)):
        projected[i,0]=W/2+np.dot(vertices[i],right)/pixel-.5;projected[i,1]=cy-np.dot(vertices[i],up)/pixel-.5;projected[i,2]=np.dot(vertices[i],eye)
    for f in faces:
        a,b,c=projected[f[0]],projected[f[1]],projected[f[2]];den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1])
        if abs(den)<1e-10:continue
        lo=max(0,int(np.floor(min(a[0],b[0],c[0]))));hi=min(W-1,int(np.ceil(max(a[0],b[0],c[0]))));ly=max(0,int(np.floor(min(a[1],b[1],c[1]))));hy=min(H-1,int(np.ceil(max(a[1],b[1],c[1]))))
        for y in range(ly,hy+1):
            for x in range(lo,hi+1):
                w0=((b[1]-c[1])*(x-c[0])+(c[0]-b[0])*(y-c[1]))/den;w1=((c[1]-a[1])*(x-c[0])+(a[0]-c[0])*(y-c[1]))/den;w2=1-w0-w1
                if min(w0,w1,w2)<-1e-8:continue
                z=w0*a[2]+w1*b[2]+w2*c[2]
                if z<=depth[y,x]:continue
                normal=w0*normals[f[0]]+w1*normals[f[1]]+w2*normals[f[2]];mu=max(0.,np.dot(normal,eye)/max(1e-30,np.sqrt((normal*normal).sum())))
                if mu<=0:continue
                depth[y,x]=z;mask[y,x]=1
                for ch in range(3):rgb[y,x,ch]=((w0*flux[f[0]]*colors[f[0],ch]+w1*flux[f[1]]*colors[f[1],ch]+w2*flux[f[2]]*colors[f[2],ch]))*(.4+.6*mu)/.8
    return rgb,depth,mask,area

def angular_operator(dirs,omega,tdiff,dt):
    """Positive conservative angular diffusion approximation.

    U=energy per solid angle. A symmetric conductance graph approximates
    diffusion on a common inner radiation reservoir. 1/tdiff is a prescribed
    effective angular diffusivity, not a multidimensional FLD solution.
    """
    from scipy.sparse import coo_matrix,diags
    tree=cKDTree(dirs);dist,near=tree.query(dirs,k=9)
    i=np.repeat(np.arange(len(dirs)),8);j=near[:,1:].ravel();lo=np.minimum(i,j);hi=np.maximum(i,j);pairs=np.unique(np.column_stack([lo,hi]),axis=0);i,j=pairs.T
    angle=2*np.arcsin(np.minimum(1,np.linalg.norm(dirs[i]-dirs[j],axis=1)/2));D=1/np.maximum(tdiff,1.);conduct=np.sqrt(omega[i]*omega[j])*(2*D[i]*D[j]/np.maximum(D[i]+D[j],1e-100))/(8*np.maximum(angle,.005)**2)
    rows=np.r_[i,j,i,j];cols=np.r_[i,j,j,i];values=np.r_[conduct,conduct,-conduct,-conduct];lap=coo_matrix((values,(rows,cols)),shape=(len(dirs),len(dirs))).tocsr()
    return diags(omega)+dt*lap

def angular_exchange(E,dirs,omega,tdiff,dt):
    from scipy.sparse.linalg import spsolve
    if dt<=0:return E.copy()
    u=spsolve(angular_operator(dirs,omega,tdiff,dt).tocsc(),E);return omega*u
