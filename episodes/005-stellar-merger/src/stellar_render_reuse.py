"""Generic rendering adapted from SiderAtlas No.004 (MIT; retained source/license).

The six-chart coordinates, fixed-Y color grade, Airy OTF, Moffat transfer,
FFT convolution and smooth tone mapping reuse No.004 algorithms. This module
provides generic imaging; the selected QualityRenderer separately uses explicitly
disclosed solar-derived dimensionless surface morphology, not solar thermodynamics.
"""
from pathlib import Path
import numpy as np
from scipy.fft import rfft2,irfft2,rfftfreq,fftfreq
ROOT=Path(__file__).resolve().parents[1]
YROW=np.array([.2126729,.7151522,.0721750])

def grade_photosphere(physical,reference,cfg):
    """No.004 golden_orbit/render.py fixed grade, preserving linear CIE Y."""
    y=physical@YROW
    relative=np.divide(y,reference,out=np.zeros_like(y),where=reference>0)
    knots=np.array(cfg['relative_intensity_knots']);gains=np.array(cfg['linear_rgb_gains_at_knots'])
    grade=np.stack([np.interp(relative,knots,gains[:,ch]) for ch in range(3)],axis=-1)
    rgb=physical*grade;gy=rgb@YROW
    rgb*=np.divide(y,gy,out=np.zeros_like(y),where=gy>0)[...,None]
    return rgb

def sample_six_charts(normals,fields,times,physical_time,scale,delay_phase=0.):
    """Vectorized adaptation of No.004's overlapping angular chart sampler.

    Input is a NEW statistically constrained scalar intensity atlas, not solar
    MHD/temperature data. Fixed tile transforms and independent phases are the
    same sampling approach used in No.004's photosphere/surface.py.
    """
    n=np.asarray(normals);result=np.zeros(len(n));weights=abs(n)**24;weights/=np.maximum(weights.sum(axis=1,keepdims=True),1e-30);weights[weights<=1e-5]=0.;weights/=np.maximum(weights.sum(axis=1,keepdims=True),1e-30)
    ny,nx=fields.shape[1:];pitch=8.;half=5.
    for ax in range(3):
        selected=weights[:,ax]>1e-5
        if not np.any(selected):continue
        m=n[selected];p=(ax+1)%3;q=(ax+2)%3;face=2*ax+(m[:,ax]<0)
        u=np.arctan2(m[:,p],abs(m[:,ax]))*scale;v=np.arctan2(m[:,q],abs(m[:,ax]))*scale
        cellx=np.floor(u/pitch).astype(np.int64);celly=np.floor(v/pitch).astype(np.int64)
        acc=np.zeros(len(m));total=np.zeros(len(m))
        for dx in [0,1]:
            for dy in [0,1]:
                tx=cellx+dx;ty=celly+dy;du=u-tx*pitch;dv=v-ty*pitch
                window=lambda a:np.where(abs(a)<=3.,1.,np.where(abs(a)<half,.5*(1+np.cos(np.pi*(abs(a)-3.)/2.)),0.))
                w=window(du)*window(dv);use=w>0
                if not np.any(use):continue
                key=(((tx+32768)*73856093)^((ty+32768)*19349663)^((face+41005)*83492791))&2147483647
                rot=key%4;uu=np.where(rot==1,-dv,np.where(rot==2,-du,np.where(rot==3,dv,du)))
                vv=np.where(rot==1,du,np.where(rot==2,-dv,np.where(rot==3,-du,dv)))
                x=(uu+half)/10*nx+(key//4)%nx;y=(vv+half)/10*ny+(key//1024)%ny
                ix=np.floor(x).astype(int);iy=np.floor(y).astype(int);fx=x-ix;fy=y-iy
                t=np.clip(physical_time+((key//1048576)%6)*.45+delay_phase,times[0],times[-1])
                it=np.clip(np.searchsorted(times,t)-1,0,len(times)-2);ft=(t-times[it])/(times[it+1]-times[it])
                value=np.zeros(len(m))
                for dt in [0,1]:
                    i0=ix%nx;i1=(ix+1)%nx;j0=iy%ny;j1=(iy+1)%ny
                    # Gather only the four needed texels per ray, never an
                    # entire atlas plane for every visible surface point.
                    val=(1-fy)*((1-fx)*fields[it+dt,j0,i0]+fx*fields[it+dt,j0,i1])+fy*((1-fx)*fields[it+dt,j1,i0]+fx*fields[it+dt,j1,i1])
                    value+=val*((1-ft) if dt==0 else ft)
                acc+=w*value;total+=w
        result[selected]+=weights[selected,ax]*np.divide(acc,total,out=np.ones_like(acc),where=total>0)
    return result

def airy_otf(frequency,first_zero):
    rho=np.minimum(1,frequency*(np.pi/3.8317059702075125)*first_zero)
    return 2/np.pi*(np.arccos(rho)-rho*np.sqrt(np.maximum(0,1-rho*rho)))

class ReusedOptics:
    def __init__(self,shape,teff=5018.53,exposure=.85):
        self.h,self.w=shape;self.pad=round(100*self.w/1080);self.shape=(self.h+2*self.pad,self.w+2*self.pad)
        freq=np.hypot(fftfreq(self.shape[0])[:,None],rfftfreq(self.shape[1])[None,:]);scale=self.w/1080
        wave=np.linspace(380,800,33);cmf=np.loadtxt(ROOT/'data/spectral/CIE_xyz_1931_2deg.csv',delimiter=',')
        lam=wave*1e-9;B=2*6.62607015e-34*299792458.**2/lam**5/np.expm1(6.62607015e-34*299792458./(lam*1.380649e-23*teff))
        weights=B*np.interp(wave,cmf[:,0],cmf[:,2]);weights[[0,-1]]*=.5;weights/=weights.sum()
        self.core=sum(w*airy_otf(freq,1.3*scale*l/550) for l,w in zip(wave,weights))*np.exp(-2*np.pi**2*(.25*scale)**2*freq**2)
        q=2*np.pi*75*scale*freq;self.wing=(1+q)*np.exp(-q);self.exposure=exposure;self.flux_errors=[];self.crop_fractions=[]
    def convolve(self,plane,otf):
        return irfft2(rfft2(np.pad(plane,self.pad),workers=2)*otf,s=self.shape,workers=2)
    def process(self,rgb):
        if rgb.sum()<=0:return np.zeros_like(rgb)
        core=np.stack([self.convolve(rgb[...,c],self.core) for c in range(3)],-1)
        halo=np.stack([self.convolve(rgb[...,c],self.wing) for c in range(3)],-1)
        eta=.045;physical=(1-eta)*core+eta*halo
        self.flux_errors.append(float(abs(physical.sum()/rgb.sum()-1)));p=self.pad
        self.crop_fractions.append(float(physical[p:p+self.h,p:p+self.w].sum()/physical.sum()))
        # No.004's enhanced imaging halo is retained as a display aid.
        hdr=np.maximum(((1-eta)*core+eta*2.2*halo)[p:p+self.h,p:p+self.w],0)
        y=hdr@YROW;l=-np.expm1(-self.exposure*y);out=hdr*np.divide(l,y,out=np.zeros_like(y),where=y>1e-15)[...,None]
        return out/np.maximum(1.,out.max(axis=-1))[...,None]
