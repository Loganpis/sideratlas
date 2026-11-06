"""Mass-carrying coarse kernels and trajectory-advected emissivity reconstruction.

Fine kernels are an explicit sub-resolution VISUAL interpolation, not extra
mass, simulated turbulence, or new fluid trajectories. No random pattern.
"""
import numpy as np
from scipy.spatial import cKDTree
from numba import njit

Y=np.array([.2126729,.7151522,.0721750])

def warm_chromaticity(T):
    # Fixed monotone compression of physical temperature into a warm display
    # color-temperature range. This is a grade, not a temperature measurement.
    global _WARM
    try:temps,colors=_WARM
    except NameError:
        from pathlib import Path
        cmf=np.loadtxt(Path(__file__).resolve().parents[2]/'data/spectral/CIE_xyz_1931_2deg.csv',delimiter=',');wave=cmf[:,0]*1e-9
        temps=np.geomspace(30,1e7,1024);mapped=1700+3500/(1+np.exp(-.9*np.log(temps/12000)))
        B=2*6.62607015e-34*299792458.**2/wave[None,:]**5/np.expm1(6.62607015e-34*299792458./(wave[None,:]*1.380649e-23*mapped[:,None]))
        xyz=np.trapz(B[:,:,None]*cmf[None,:,1:4],wave,axis=1);matrix=np.array([[3.2404542,-1.5371385,-.4985314],[-.969266,1.8760108,.041556],[.0556434,-.2040259,1.0572252]])
        colors=np.maximum(0,xyz@matrix.T);colors/=np.maximum(colors@Y,1e-30)[:,None];_WARM=(temps,colors)
    return np.stack([np.interp(np.log(np.maximum(T,30)),np.log(temps),colors[:,j]) for j in range(3)],-1)

def coarse_sigma(x,radius,center):
    if len(x)<2:return np.full(len(x),.1*radius)
    dist=cKDTree(x).query(x,k=min(17,len(x)),workers=2)[0];h=dist[:,-1] if dist.ndim==2 else np.full(len(x),.1*radius)
    # Gaussian sigma, not full support. Uses 16 neighbors, vs previous 32.
    return np.maximum(.025,np.minimum(h/3.,.13*np.maximum(np.linalg.norm(x-center,axis=1),radius)))

def fine_covariances(motion,t,q,ids,sigma):
    """Finite deformation fit to same neighboring material at t-6 and t.

    A material-phase filter rejects unrelated launch families and velocity
    branches. Neighbor identities never form connecting segments.
    """
    p=q['x'][ids];v=q['v'][ids];n=len(p)
    if n<2:return np.eye(3)[None]*np.maximum(sigma[:,None,None],.025)**2,dict(fitted=0,parcels=n)
    tau=6.;old=motion.at(max(0,t-tau));op=old['x'][ids];ov=old['v'][ids]
    k=min(32,n);distance,near=cKDTree(p).query(p,k=k,workers=2)
    if k==1:near=near[:,None];distance=distance[:,None]
    dx=op[near]-op[:,None,:];dy=p[near]-p[:,None,:];dv=ov[near]-ov[:,None,:]
    birth=motion.birth[ids];same=(abs(birth[near]-birth[:,None])<8.)&old['live'][ids[near]]&old['live'][ids,None]
    velocity_scale=np.maximum(.04,np.linalg.norm(ov,axis=1)*.45)
    w=np.exp(-distance**2/(2*np.maximum(3*sigma[:,None],.03)**2))*same*np.exp(-np.sum(dv*dv,axis=-1)/(2*velocity_scale[:,None]**2));w[:,0]=0
    cov0=np.einsum('nk,nki,nkj->nij',w,dx,dx);cross=np.einsum('nk,nki,nkj->nij',w,dy,dx);ridge=(np.trace(cov0,axis1=1,axis2=2)*.005+1e-8)
    F=np.einsum('nij,njk->nik',cross+ridge[:,None,None]*np.eye(3),np.linalg.inv(cov0+ridge[:,None,None]*np.eye(3)))
    fit=(np.sum(w>.05,axis=1)>=5)&(birth<t-tau);F[~fit]=np.eye(3)
    U,s,V=np.linalg.svd(F);s=np.clip(s,.3,3.);limited=U*s[:,None,:]
    C=np.einsum('nij,nkj->nik',limited,limited)
    # Current velocity supplies only a material orientation fallback, never a
    # rigid source/star angle or a long trail. Fine axes have capped aspect.
    tangent=v/np.maximum(np.linalg.norm(v,axis=1,keepdims=True),1e-12)
    ev,evec=np.linalg.eigh(C);across=evec[:,:,0]-np.sum(evec[:,:,0]*tangent,axis=1)[:,None]*tangent
    bad=np.linalg.norm(across,axis=1)<1e-6
    across[bad]=np.cross(tangent[bad],np.array([.12,.98,.15]));across/=np.maximum(np.linalg.norm(across,axis=1,keepdims=True),1e-12)
    other=np.cross(tangent,across);basis=np.stack([tangent,across,other],axis=-1)
    narrow=np.maximum(.025,sigma*.30*np.sqrt(np.clip(ev[:,0],.4,2)))
    major=np.maximum(narrow,np.minimum(sigma*np.sqrt(np.clip(ev[:,-1],.8,4)),narrow*4))
    axes=np.stack([major,narrow,narrow*np.sqrt(np.clip(ev[:,1]/np.maximum(ev[:,0],1e-10),1,2))],axis=1)
    result=np.einsum('nij,nkj->nik',basis*axes[:,None,:],basis*axes[:,None,:])
    return result,dict(fitted=int(fit.sum()),parcels=n,median_axis_ratio=float(np.median(major/narrow)),fine_transverse_display_factor=.30,max_axis_ratio=4.)

@njit(cache=True)
def splat_emission(p,cov,power,colors,depth,W,H,pixel,cy,occlude=True):
    """Analytic orthographic Gaussian column, EWA pixel response.

    Per-source discrete footprint is normalized before viewport clipping.
    Luminosity is preserved even for tiny native footprints.
    """
    out=np.zeros((H,W,3),np.float32);in_power=0.
    for j in range(len(p)):
        a=cov[j,0,0]/pixel**2+.25;b=cov[j,0,1]/pixel**2;c=cov[j,1,1]/pixel**2+.25;det=max(a*c-b*b,1e-12);ia=c/det;ib=-b/det;ic=a/det
        cx=W/2+p[j,0]/pixel-.5;yy=cy-p[j,1]/pixel-.5;sx=3.5*np.sqrt(a);sy=3.5*np.sqrt(c);lo=int(np.floor(cx-sx));hi=int(np.ceil(cx+sx));ly=int(np.floor(yy-sy));hy=int(np.ceil(yy+sy));norm=0.
        if hi<0 or lo>=W or hy<0 or ly>=H:continue
        for y in range(ly,hy+1):
            for x in range(lo,hi+1):
                dx=x-cx;dy=-(y-yy);e=ia*dx*dx+2*ib*dx*dy+ic*dy*dy
                if e<12.25:norm+=np.exp(-.5*e)
        if norm<=0:continue
        for y in range(max(0,ly),min(H-1,hy)+1):
            for x in range(max(0,lo),min(W-1,hi)+1):
                dx=x-cx;dy=-(y-yy);e=ia*dx*dx+2*ib*dx*dy+ic*dy*dy
                if e>=12.25:continue
                if occlude and p[j,2]<depth[y,x]:continue
                value=power[j]*np.exp(-.5*e)/norm
                for ch in range(3):out[y,x,ch]+=value*colors[j,ch]
                in_power+=value
    return out,in_power

@njit(cache=True)
def optical_depth(p,sigma,mass,depth,W,H,pixel,cy,kappa_column):
    import_dummy=0
    result=np.zeros((H,W),np.float32)
    # erfc is imported at module scope below; center-depth conditional for
    # isotropic mass kernels makes foreground fraction analytic.
    for j in range(len(p)):
        s=max(.025,sigma[j]);s2=s*s+pixel*pixel/12;cx=W/2+p[j,0]/pixel-.5;yy=cy-p[j,1]/pixel-.5;r=4*np.sqrt(s2)/pixel
        for y in range(max(0,int(yy-r)),min(H-1,int(yy+r+1))+1):
            for x in range(max(0,int(cx-r)),min(W-1,int(cx+r+1))+1):
                if not np.isfinite(depth[y,x]):continue
                dx=(x-cx)*pixel;dy=(y-yy)*pixel;e=(dx*dx+dy*dy)/s2
                if e>16:continue
                fraction=.5*erfc((depth[y,x]-p[j,2])/(np.sqrt(2)*s));result[y,x]+=mass[j]*np.exp(-.5*e)/(2*np.pi*s2)*kappa_column*fraction
    return result
from math import erfc

def gas_display_grade(rgb):
    """No.001-inspired warm brightness palette, retaining thermal ordering.

    This fixed display operation is separate from the physical photon ledger.
    At fixed radiance, 25% of thermal chromaticity is retained. Colors cannot
    be interpreted as measured temperatures.
    """
    y=np.maximum(0,rgb@Y);s=y/(y+.04);stops=np.array([0,.06,.18,.35,.55,.76,.92,1.]);colors=np.array([[0,0,0],[17,0,0],[73,3,0],[151,15,0],[225,48,0],[255,121,7],[255,200,46],[255,246,177]],float)/255
    c=np.stack([np.interp(s,stops,colors[:,j]) for j in range(3)],-1);c=np.where(c<=.04045,c/12.92,((c+.055)/1.055)**2.4);c/=np.maximum(c@Y,1e-30)[...,None]
    thermal=rgb/np.maximum(y,1e-30)[...,None];hue=.75*c+.25*thermal;shown=.20*(y/.20)**1.3
    return hue*shown[...,None]
