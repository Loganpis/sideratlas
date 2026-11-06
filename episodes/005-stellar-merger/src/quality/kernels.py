"""Rendering kernels adapted from No.004; no fluid or magnetic evolution."""
import numpy as np
from numba import njit,prange

@njit(cache=True,inline='always')
def bilinear(a,it,im,ia,x,y):
    ny,nx=a.shape[-2:];ix=int(np.floor(x));iy=int(np.floor(y));fx=x-ix;fy=y-iy
    return ((1-fy)*((1-fx)*a[it,im,ia,iy%ny,ix%nx]+fx*a[it,im,ia,iy%ny,(ix+1)%nx])+fy*((1-fx)*a[it,im,ia,(iy+1)%ny,ix%nx]+fx*a[it,im,ia,(iy+1)%ny,(ix+1)%nx]))/4096.

@njit(parallel=True,cache=True)
def surface_sample(n,observer,mu,fields,times,physical_time,scale,mu_grid,meso,spot_fields,spot_times,spot_time,spot_center,east,north):
    """No.004 overlapping chart sampler, scalar normalized morphology only."""
    out=np.zeros((len(n),3),np.float32);ny,nx=fields.shape[-2:];pitch=19.2;half=12.
    for j in prange(len(n)):
        weights=np.abs(n[j])**24;weights/=weights.sum();weights[weights<1e-5]=0.;weights/=weights.sum();value=0.;large=0.
        nod=np.dot(n[j],observer)
        for ax in range(3):
            if weights[ax]<1e-5:continue
            p=(ax+1)%3;q=(ax+2)%3;face=2*ax+int(n[j,ax]<0);sgn=1. if n[j,ax]>=0 else -1.
            u=np.arctan2(n[j,p],abs(n[j,ax]))*scale;v=np.arctan2(n[j,q],abs(n[j,ax]))*scale
            norm=np.sqrt(max(1e-30,1-n[j,p]**2));ou=(observer[p]-n[j,p]*nod)/norm;ov=sgn*(n[j,ax]*observer[q]-n[j,q]*observer[ax])/norm
            ph=np.arctan2(ov,ou);im=max(0,min(len(mu_grid)-2,np.searchsorted(mu_grid,mu[j])-1));fm=min(1.,max(0.,(mu[j]-mu_grid[im])/(mu_grid[im+1]-mu_grid[im])))
            cellx=int(np.floor(u/pitch));celly=int(np.floor(v/pitch));acc=0.;lacc=0.;total=0.
            for tx in range(cellx,cellx+2):
                du=u-tx*pitch;wu=1. if abs(du)<=7.2 else (.5*(1+np.cos(np.pi*(abs(du)-7.2)/4.8)) if abs(du)<half else 0.)
                if wu<=0:continue
                for ty in range(celly,celly+2):
                    dv=v-ty*pitch;wv=1. if abs(dv)<=7.2 else (.5*(1+np.cos(np.pi*(abs(dv)-7.2)/4.8)) if abs(dv)<half else 0.)
                    w=wu*wv
                    if w<=0:continue
                    h=(((tx+32768)*73856093)^((ty+32768)*19349663)^((face+41004)*83492791))&2147483647;rot=h%4
                    uu=du;vv=dv
                    if rot==1:uu=-dv;vv=du
                    elif rot==2:uu=-du;vv=-dv
                    elif rot==3:uu=dv;vv=-du
                    x=(uu+half)/24*nx+(h//4)%nx;y=(vv+half)/24*ny+(h//1024)%ny
                    t=physical_time+((h//1048576)%6)*20.
                    it=max(0,min(len(times)-2,np.searchsorted(times,t)-1));ft=min(1.,max(0.,(t-times[it])/(times[it+1]-times[it])))
                    af=((ph+rot*np.pi/2)%(2*np.pi))/(np.pi/2);ia=int(af)%4;fa=af-np.floor(af)
                    for dt in range(2):
                        wt=(1-ft) if dt==0 else ft
                        for dm in range(2):
                            wm=(1-fm) if dm==0 else fm
                            for da in range(2):acc+=w*wt*wm*((1-fa) if da==0 else fa)*bilinear(fields,it+dt,im+dm,(ia+da)%4,x,y)
                    # A lower-frequency component of the same source sequence.
                    ix=int(np.floor(x));iy=int(np.floor(y));fx=x-ix;fy=y-iy
                    for dt in range(2):
                        a=meso[it+dt];val=(1-fy)*((1-fx)*a[iy%ny,ix%nx]+fx*a[iy%ny,(ix+1)%nx])+fy*((1-fx)*a[(iy+1)%ny,ix%nx]+fx*a[(iy+1)%ny,(ix+1)%nx])
                        lacc+=w*((1-ft) if dt==0 else ft)*val
                    total+=w
            if total>0:value+=weights[ax]*acc/total;large+=weights[ax]*lacc/total
        spot=1.;front=np.dot(n[j],spot_center)
        if front>0:
            u=np.arctan2(np.dot(n[j],east),front);v=np.arctan2(np.dot(n[j],north),front)
            if abs(u)<.32 and abs(v)<.16:
                wx=1. if abs(u)<.26 else .5*(1+np.cos(np.pi*(abs(u)-.26)/.06))
                wy=1. if abs(v)<.12 else .5*(1+np.cos(np.pi*(abs(v)-.12)/.04))
                yy=(v/.32+.5)*(spot_fields.shape[1]-1);xx=(u/.64+.5)*(spot_fields.shape[2]-1)
                iy=max(0,min(spot_fields.shape[1]-2,int(yy)));ix=max(0,min(spot_fields.shape[2]-2,int(xx)));fy=yy-iy;fx=xx-ix
                it=max(0,min(len(spot_times)-2,np.searchsorted(spot_times,spot_time)-1));ft=min(1.,max(0.,(spot_time-spot_times[it])/(spot_times[it+1]-spot_times[it])))
                ratio=0.
                for dt in range(2):
                    a=spot_fields[it+dt];val=(1-fy)*((1-fx)*a[iy,ix]+fx*a[iy,ix+1])+fy*((1-fx)*a[iy+1,ix]+fx*a[iy+1,ix+1])
                    ratio+=((1-ft) if dt==0 else ft)*val
                spot=1+wx*wy*(min(1.025,max(.78,max(.05,ratio)**.25))-1)
        out[j,0]=value;out[j,1]=large;out[j,2]=spot
    return out

@njit(cache=True,inline='always')
def compact(q):
    if q>=1:return 0.
    return (1-q)**6*(1+6*q+35*q*q/3)

@njit(cache=True,inline='always')
def lattice_value(field,x,y,z):
    n=field.shape[0];x=x*n/12;y=y*n/12;z=z*n/12
    ix=int(np.floor(x));iy=int(np.floor(y));iz=int(np.floor(z));fx=x-ix;fy=y-iy;fz=z-iz;val=0.
    for k in range(2):
        for j in range(2):
            for i in range(2):val+=field[(iz+k)%n,(iy+j)%n,(ix+i)%n]*(fz if k else 1-fz)*(fy if j else 1-fy)*(fx if i else 1-fx)
    return val

@njit(cache=True,inline='always')
def material_eta(c0,c1,c2,phase,field):
    if field is None:
        return (np.cos(1.7*c0+2.3*c1+.7*c2)*np.cos(phase)+.6*np.cos(3.2*c0-.8*c1+2.7*c2)*np.cos(2*phase+.71))
    # Frozen, correlated material field; even pairing retains the packet's
    # continuous centroid. Neighbor launch phases vary continuously.
    x=phase*.91;y=phase*.47;z=phase*1.13
    return .7071067811865476*(lattice_value(field,x+c0,y+c1,z+c2)+lattice_value(field,x-c0,y-c1,z-c2))

@njit(cache=True)
def deposit_packets(positions,axes,sigmas,masses,temperatures,clump_sigma,phase,shape,origin,spacing,zcenters,zsizes,support,material_axes=None,clump_field=None):
    """Conservative compact flow packets; fixed material-space clumping.

    Full stencil normalization includes virtual off-camera cells, preserving
    source mass without renormalizing clipped material into the picture.
    """
    mass=np.zeros(shape,np.float32);heat=np.zeros(shape,np.float32);N=len(positions)
    inside=0.
    for j in range(N):
        p=positions[j];sig=sigmas[j];basis=axes[j]
        # Shape eigenvectors can flip signs or exchange order when eigenvalues
        # cross. Statistical structure instead follows continuous flow axes.
        material=basis if material_axes is None else material_axes[j]
        correlation=max(.06,min(sig[0],min(sig[1],sig[2])))
        # Covariance voxel broadening is an explicit rendering-resolution floor.
        extent=np.sqrt(np.sum((basis*sig[None,:])**2,axis=1))*support
        lower=np.floor((p-extent-origin)/spacing-.5).astype(np.int64);upper=np.ceil((p+extent-origin)/spacing-.5).astype(np.int64)
        lower[2]=np.searchsorted(zcenters,p[2]-extent[2])-1;upper[2]=np.searchsorted(zcenters,p[2]+extent[2])
        if p[2]-extent[2]<zcenters[0]:lower[2]=int(np.floor((p[2]-extent[2]-zcenters[0])/zsizes[0]))
        if p[2]+extent[2]>zcenters[-1]:upper[2]=len(zcenters)-1+int(np.ceil((p[2]+extent[2]-zcenters[-1])/zsizes[-1]))
        if upper[0]<0 or upper[1]<0 or upper[2]<0 or lower[0]>=shape[2] or lower[1]>=shape[1] or lower[2]>=shape[0]:continue
        norm=0.;tnorm=0.
        for z in range(lower[2],upper[2]+1):
            for y in range(lower[1],upper[1]+1):
                for x in range(lower[0],upper[0]+1):
                    zz=zcenters[z] if 0<=z<len(zcenters) else (zcenters[0]+z*zsizes[0] if z<0 else zcenters[-1]+(z-len(zcenters)+1)*zsizes[-1])
                    dz=zsizes[z] if 0<=z<len(zcenters) else (zsizes[0] if z<0 else zsizes[-1])
                    dx=origin[0]+(x+.5)*spacing[0]-p[0];dy=origin[1]+(y+.5)*spacing[1]-p[1];ddz=zz-p[2]
                    u0=(dx*basis[0,0]+dy*basis[1,0]+ddz*basis[2,0])/sig[0]
                    u1=(dx*basis[0,1]+dy*basis[1,1]+ddz*basis[2,1])/sig[1]
                    u2=(dx*basis[0,2]+dy*basis[1,2]+ddz*basis[2,2])/sig[2]
                    q=np.sqrt(u0*u0+u1*u1+u2*u2)/support
                    if q>=1:continue
                    c0=(dx*material[0,0]+dy*material[1,0]+ddz*material[2,0])/correlation
                    c1=(dx*material[0,1]+dy*material[1,1]+ddz*material[2,1])/correlation
                    c2=(dx*material[0,2]+dy*material[1,2]+ddz*material[2,2])/correlation
                    eta=material_eta(c0,c1,c2,phase[j],clump_field)
                    factor=np.exp(clump_sigma[j]*eta)
                    w=compact(q)*factor*dz;norm+=w;tnorm+=w*factor**.2
        if norm<=0:continue
        thermal_norm=norm/max(tnorm,1e-30)
        for z in range(max(0,lower[2]),min(shape[0]-1,upper[2])+1):
            for y in range(max(0,lower[1]),min(shape[1]-1,upper[1])+1):
                for x in range(max(0,lower[0]),min(shape[2]-1,upper[0])+1):
                    zz=zcenters[z] if 0<=z<len(zcenters) else (zcenters[0]+z*zsizes[0] if z<0 else zcenters[-1]+(z-len(zcenters)+1)*zsizes[-1])
                    dz=zsizes[z] if 0<=z<len(zcenters) else (zsizes[0] if z<0 else zsizes[-1])
                    dx=origin[0]+(x+.5)*spacing[0]-p[0];dy=origin[1]+(y+.5)*spacing[1]-p[1];ddz=zz-p[2]
                    u0=(dx*basis[0,0]+dy*basis[1,0]+ddz*basis[2,0])/sig[0]
                    u1=(dx*basis[0,1]+dy*basis[1,1]+ddz*basis[2,1])/sig[1]
                    u2=(dx*basis[0,2]+dy*basis[1,2]+ddz*basis[2,2])/sig[2]
                    q=np.sqrt(u0*u0+u1*u1+u2*u2)/support
                    if q>=1:continue
                    c0=(dx*material[0,0]+dy*material[1,0]+ddz*material[2,0])/correlation
                    c1=(dx*material[0,1]+dy*material[1,1]+ddz*material[2,1])/correlation
                    c2=(dx*material[0,2]+dy*material[1,2]+ddz*material[2,2])/correlation
                    eta=material_eta(c0,c1,c2,phase[j],clump_field)
                    factor=np.exp(clump_sigma[j]*eta);w=masses[j]*compact(q)*factor*dz/norm
                    mass[z,y,x]+=w;heat[z,y,x]+=w*temperatures[j]*factor**.2*thermal_norm;inside+=w
    return mass,heat,inside

@njit(cache=True,inline='always')
def interpolate_xy(a,z,x,y):
    iy=int(np.floor(y));ix=int(np.floor(x));fy=y-iy;fx=x-ix
    if ix<0 or iy<0 or ix>=a.shape[2]-1 or iy>=a.shape[1]-1:return 0.
    return (1-fy)*((1-fx)*a[z,iy,ix]+fx*a[z,iy,ix+1])+fy*((1-fx)*a[z,iy+1,ix]+fx*a[z,iy+1,ix+1])

@njit(parallel=True,cache=True)
def raycast(mass,heat,star,depth,origin,spacing,pixel,halfwidth,center_y,zcenters,zsizes,kappa_mass,lut_logT,lut_rgb,donorT,dilution):
    """Native image rays through the physical 3D field, never image upscaling."""
    H,W=depth.shape;out=star.copy();surviving=np.zeros((H,W),np.float32)
    nz=mass.shape[0];first=lut_logT[0];step=lut_logT[1]-first
    for iy in prange(H):
        v=(center_y-iy-.5)*pixel;gy=(v-origin[1])/spacing[1]-.5
        for ix in range(W):
            u=(ix+.5-W/2)*pixel;gx=(u-origin[0])/spacing[0]-.5;total=0.;foreground=0.
            for z in range(nz):
                dm=interpolate_xy(mass,z,gx,gy);tau=dm*kappa_mass;total+=tau
                zz=zcenters[z];f=min(1.,max(0.,(zz+.5*zsizes[z]-depth[iy,ix])/zsizes[z]));foreground+=tau*f
            scale=min(foreground,.16*(foreground/1e5)**.38)/max(foreground,1e-40)
            g0=0.;g1=0.;g2=0.;p0=0.;p1=0.;p2=0.
            for z in range(nz):
                dm=interpolate_xy(mass,z,gx,gy)
                if dm<=1e-25:continue
                zz=zcenters[z];f=min(1.,max(0.,(zz+.5*zsizes[z]-depth[iy,ix])/zsizes[z]));true_tau=dm*kappa_mass*f;tau=true_tau*scale
                if tau<=1e-10:continue
                a=-np.expm1(-tau);temp=max(300.,interpolate_xy(heat,z,gx,gy)/dm*dilution**.25)
                ti=(np.log(temp)-first)/step;k=max(0,min(len(lut_logT)-2,int(ti)));ft=min(1.,max(0.,ti-k))
                c0=(1-ft)*lut_rgb[k,0]+ft*lut_rgb[k+1,0];c1=(1-ft)*lut_rgb[k,1]+ft*lut_rgb[k+1,1];c2=(1-ft)*lut_rgb[k,2]+ft*lut_rgb[k+1,2]
                yy=.2126729*c0+.7151522*c1+.0721750*c2
                physical_alpha=-np.expm1(-min(100.,true_tau))
                p0=p0*(1-physical_alpha)+physical_alpha*c0
                p1=p1*(1-physical_alpha)+physical_alpha*c1
                p2=p2*(1-physical_alpha)+physical_alpha*c2
                target=.12*np.arcsinh(6*(temp/donorT)**4/.12);s=target/max(yy,1e-60)
                g0=g0*(1-a)+a*c0*s;g1=g1*(1-a)+a*c1*s;g2=g2*(1-a)+a*c2*s
            # Physical grey continuum is always included with exactly the
            # opacity used to extinguish stellar light. The enhanced gas-only
            # display is an additive diagnostic of cold bolometric radiation.
            trans=np.exp(-min(100.,foreground));out[iy,ix,0]=star[iy,ix,0]*trans+p0+g0;out[iy,ix,1]=star[iy,ix,1]*trans+p1+g1;out[iy,ix,2]=star[iy,ix,2]*trans+p2+g2
            surviving[iy,ix]=(star[iy,ix,0]+star[iy,ix,1]+star[iy,ix,2])*trans
    return out,surviving

@njit(cache=True)
def bolometric_upper(mass,heat,kappa_mass,pixel_m):
    total=0.
    for y in range(mass.shape[1]):
        for x in range(mass.shape[2]):
            intensity=0.
            for z in range(mass.shape[0]):
                m=mass[z,y,x]
                if m<=1e-25:continue
                alpha=-np.expm1(-min(80.,m*kappa_mass));temp=max(300.,heat[z,y,x]/m)
                intensity=intensity*(1-alpha)+(5.670374419e-8/np.pi)*temp**4*alpha
            total+=intensity
    return total*4*np.pi*pixel_m**2
