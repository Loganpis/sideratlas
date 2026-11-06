"""Position-based mass estimator and energy-normalized cooling display.

No velocity-aligned arms, age-growth, procedural texture or class-based coloring.
"""
import numpy as np
from scipy.spatial import cKDTree
from numba import njit,prange
from src.quality.kernels import interpolate_xy

def absorbed_fraction(mass,origin,spacing,source,radius,kappa_mass):
    """Finite ray quadrature for stellar radiation intercepted by this grid."""
    from scipy.ndimage import map_coordinates
    n=384;i=np.arange(n);z=1-2*(i+.5)/n;phi=i*2.399963229728653;directions=np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
    corner=np.maximum(abs(origin-source),abs(origin+spacing*np.array(mass.shape[::-1])-source));rmax=max(radius*1.01,np.linalg.norm(corner));edges=np.geomspace(radius*1.002,rmax,129);r=(edges[:-1]+edges[1:])/2
    points=source+directions[:,None,:]*r[None,:,None];coords=(points-origin)/spacing-.5
    dm=map_coordinates(mass,[coords[...,2],coords[...,1],coords[...,0]],order=1,mode='constant');tau=np.sum(dm*kappa_mass/spacing[2]*np.diff(edges)[None,:],axis=1)
    return float(np.mean(-np.expm1(-np.minimum(tau,100))))

@njit(cache=True)
def reprocessed_upper(mass,origin,spacing,kappa_mass,source,radius,luminosity,pixel_m,rsun_m):
    """One-view isotropic-equivalent photospheric reprocessing bound."""
    total=0.
    for y in range(mass.shape[1]):
        dy=origin[1]+(y+.5)*spacing[1]-source[1]
        for x in range(mass.shape[2]):
            dx=origin[0]+(x+.5)*spacing[0]-source[0];I=0.
            for z in range(mass.shape[0]):
                alpha=-np.expm1(-min(100.,mass[z,y,x]*kappa_mass));dz=origin[2]+(z+.5)*spacing[2]-source[2];r2=max(radius*radius,dx*dx+dy*dy+dz*dz)*rsun_m**2
                S=luminosity/(4*np.pi*np.pi*r2);I=I*(1-alpha)+alpha*S
            total+=I
    return total*4*np.pi*pixel_m**2

@njit(cache=True)
def deposit(x,m,T,h,shape,origin,spacing):
    mass=np.zeros(shape,np.float32);heat=np.zeros(shape,np.float32);inside=0.
    for j in range(len(x)):
        lo=np.floor((x[j]-h[j]-origin)/spacing-.5).astype(np.int64);hi=np.ceil((x[j]+h[j]-origin)/spacing-.5).astype(np.int64);norm=0.
        for z in range(lo[2],hi[2]+1):
            for y in range(lo[1],hi[1]+1):
                for k in range(lo[0],hi[0]+1):
                    dx=origin[0]+(k+.5)*spacing[0]-x[j,0];dy=origin[1]+(y+.5)*spacing[1]-x[j,1];dz=origin[2]+(z+.5)*spacing[2]-x[j,2];q=np.sqrt(dx*dx+dy*dy+dz*dz)/h[j]
                    if q<1:norm+=(1-q)**4*(1+4*q)
        if norm<=0:continue
        for z in range(max(0,lo[2]),min(shape[0]-1,hi[2])+1):
            for y in range(max(0,lo[1]),min(shape[1]-1,hi[1])+1):
                for k in range(max(0,lo[0]),min(shape[2]-1,hi[0])+1):
                    dx=origin[0]+(k+.5)*spacing[0]-x[j,0];dy=origin[1]+(y+.5)*spacing[1]-x[j,1];dz=origin[2]+(z+.5)*spacing[2]-x[j,2];q=np.sqrt(dx*dx+dy*dy+dz*dz)/h[j]
                    if q<1:
                        w=m[j]*(1-q)**4*(1+4*q)/norm;mass[z,y,k]+=w;heat[z,y,k]+=w*T[j];inside+=w
    return mass,heat,inside

def adaptive_grid(points,masses,temperatures,halfwidth,gw=128,nz=96,center_y=.64,eye=None,right=None,up=None,reference_center=None,reference_radius=3.72):
    if eye is None:eye=np.array([0.,0.,1.]);right=np.array([1.,0.,0.]);up=np.array([0.,1.,0.])
    # Support is the current 32nd-neighbor distance, not an orbital phase or
    # a trail. The same physical scalar support is used for every camera.
    k=min(33,len(points));dist=cKDTree(points).query(points,k=k,workers=2)[0]
    support=dist[:,-1] if dist.ndim==2 else np.zeros(len(points))
    center=np.zeros(3) if reference_center is None else reference_center
    support=np.minimum(support,.35*np.maximum(np.linalg.norm(points-center,axis=1),reference_radius))
    camera=np.column_stack([right,up,eye]);p=points@camera;gh=round(gw*16/9);pad=4;dx=2*halfwidth/gw
    Z=max(halfwidth,float(np.quantile(abs(p[:,2]),.995))+dx) if len(p) else halfwidth;dz=2*Z/nz
    # A voxel-response floor prevents missing a sub-resolution source. It is
    # explicitly numerical and changes only the estimator, not parcel motion.
    support=np.maximum(support,2.2*max(dx,dz));origin=np.array([-halfwidth-pad*dx,-(1-center_y)*gh*dx-pad*dx,-Z-pad*dz]);spacing=np.array([dx,dx,dz]);shape=(nz+2*pad,gh+2*pad,gw+2*pad)
    m,h,inside=deposit(p.astype('f8'),np.asarray(masses,'f8'),np.maximum(0.,temperatures).astype('f8'),support.astype('f8'),shape,origin,spacing)
    zc=origin[2]+(np.arange(shape[0])+.5)*dz;zs=np.full(shape[0],dz)
    return m,h,origin,spacing,zc,zs,dict(support_median_Rsun=float(np.median(support)),support_max_Rsun=float(np.max(support)),mass_inside=float(inside))

@njit(parallel=True,cache=True)
def cooling_rays(mass,heat,star,depth,origin,spacing,pixel,L,center_y,zc,zs,kappa_mass,lut_logT,lut_rgb,donorT,cooling_s,grey_kappa_SI,cv,display_gain,source_center,source_radius,stellar_luminosity,reprocess_scale,rsun_m,reprocess_display_gain=10.):
    H,W=depth.shape;out=star.copy();surviving=np.zeros((H,W),np.float32);first=lut_logT[0];step=lut_logT[1]-first
    for y in prange(H):
        gy=((center_y-y-.5)*pixel-origin[1])/spacing[1]-.5
        for x in range(W):
            gx=((x+.5-W/2)*pixel-origin[0])/spacing[0]-.5;fore=0.
            for z in range(mass.shape[0]):
                dm=interpolate_xy(mass,z,gx,gy);f=min(1.,max(0.,(zc[z]+.5*zs[z]-depth[y,x])/zs[z]));fore+=dm*kappa_mass*f
            scale=min(fore,.22*(fore/1e5)**.38)/max(fore,1e-40);g=np.zeros(3);rep=np.zeros(3);physical=np.zeros(3)
            for z in range(mass.shape[0]):
                dm=interpolate_xy(mass,z,gx,gy)
                if dm<=1e-25:continue
                f=min(1.,max(0.,(zc[z]+.5*zs[z]-depth[y,x])/zs[z]));tau=dm*kappa_mass*f
                if tau<=1e-12:continue
                trueT=max(0.,interpolate_xy(heat,z,gx,gy)/dm);T=max(300.,trueT);ti=(np.log(T)-first)/step;k=max(0,min(len(lut_logT)-2,int(ti)));ft=min(1.,max(0.,ti-k));c=(1-ft)*lut_rgb[k]+ft*lut_rgb[k+1];Y=.2126729*c[0]+.7151522*c[1]+.0721750*c[2]
                a=-np.expm1(-min(100.,tau));ad=-np.expm1(-min(100.,tau*scale))
                # Finite volume cooling, j_bol = rho cv T / (4 pi t_cool).
                # Thermal spectral shape is normalized to that cooling power;
                # a mass sample is never treated as a blackbody surface.
                source_factor=cv*trueT/(4*grey_kappa_SI*cooling_s*5.670374419e-8*T**4)
                shown=display_gain*cv*trueT/(4*grey_kappa_SI*cooling_s*5.670374419e-8*donorT**4)
                dx=(x+.5-W/2)*pixel-source_center[0];dy=(center_y-y-.5)*pixel-source_center[1];dz=zc[z]-source_center[2];r2=max(source_radius**2,dx*dx+dy*dy+dz*dz)*rsun_m**2
                # Effective emitting layer in irradiated optically thick gas.
                # It redistributes intercepted stellar power, not new heat.
                irrT=max(300.,(stellar_luminosity/(4*np.pi*5.670374419e-8*r2))**.25);it=(np.log(irrT)-first)/step;ik=max(0,min(len(lut_logT)-2,int(it)));iff=min(1.,max(0.,it-ik));irgb=(1-iff)*lut_rgb[ik]+iff*lut_rgb[ik+1];iY=.2126729*irgb[0]+.7151522*irgb[1]+.0721750*irgb[2]
                irrshown=reprocess_scale*stellar_luminosity/(4*np.pi*5.670374419e-8*r2*donorT**4)
                for ch in range(3):
                    physical[ch]=physical[ch]*(1-a)+a*(c[ch]*source_factor+reprocess_scale*irgb[ch])
                    # The scientific overlay is a transparent integral of
                    # finite cooling emissivity, not an opaque kernel surface.
                    # Tone compression is applied once to the integrated ray.
                    g[ch]+=tau*shown*c[ch]/max(Y,1e-60)
                    # Same thermal layer shown bolometrically, including IR;
                    # fixed one-times gain, distinct from optical photometry.
                    rep[ch]=rep[ch]*(1-a)+a*reprocess_display_gain*irrshown*irgb[ch]/max(iY,1e-60)
            trans=np.exp(-min(100.,fore))
            # One scalar luminance map preserves temperature chromaticity.
            luminance=.2126729*g[0]+.7151522*g[1]+.0721750*g[2]
            mapped=.06*np.arcsinh(luminance/.06);gain=mapped/max(luminance,1e-60)
            for ch in range(3):out[y,x,ch]=star[y,x,ch]*trans+physical[ch]+gain*g[ch]+rep[ch]
            surviving[y,x]=star[y,x].sum()*trans
    return out,surviving
