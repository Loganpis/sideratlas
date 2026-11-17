"""Approximate IMF/Planck population builder and accepted frozen lookups."""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from . import ROOT, CONFIG
from .spectrum import make_lut,color,XYZ_RGB,H,C,K
SIGMA=5.670374419e-8
LSUN=3.828e33
YROW=np.array([.2126729,.7151522,.0721750])
def smooth(u):
    u=np.clip(u,0,1); return u*u*(3-2*u)

class Populations:
    def __init__(self):
        self.lut=make_lut()
        solar=color(5772.,self.lut)/(SIGMA*5772.**4/np.pi)
        self.norm=float(solar@YROW)
        cmf=np.loadtxt(ROOT/'data/spectral/CIE_xyz_1931_2deg.csv',delimiter=',')
        wave=cmf[:,0]*1e-9
        B=2*H*C*C/wave**5/np.expm1(H*C/(wave*K*5772.))
        self.y_per_bol=float(np.trapz(B*cmf[:,2],wave)/(SIGMA*5772.**4/np.pi))
        self.line_rgb={}
        for name,wavelength in [('ha',656.3),('hb',486.1)]:
            xyz=np.array([np.interp(wavelength,cmf[:,0],cmf[:,j]) for j in (1,2,3)])
            self.line_rgb[name]=np.maximum(0,XYZ_RGB@xyz)/self.y_per_bol
        # Same log-log interpolation of the actual No.006 coefficients;
        # no nebular geometry, wind model or opacity code is imported.
        self.alpha_B=2.59e-13
        self.bulge_rgb=self.rgb_per_Lbol(4500.)

        m=np.geomspace(.08,100.,600)
        imf=np.where(m<.5,m**(-1.3),.5*m**(-2.3))
        counts=imf*np.gradient(m); counts/=np.sum(counts*m)
        L=np.where(m<.43,.23*m**2.3,np.where(m<2,m**4,
            np.where(m<20,1.5*m**3.5,1.5*20**3.5*(m/20)**1.8)))
        radius=np.where(m<1,m**.8,m**.6)
        T=5772.*(L/radius**2)**.25
        life=np.maximum(3.,1e4*m/L*(1+.8*smooth((m-2)/8)))
        self.ages=np.unique(np.r_[np.linspace(0,15,151),np.linspace(15,100,171),np.linspace(100,1000,361)])
        alive=1/(1+np.exp(np.clip((self.ages[:,None]-life[None,:])/(.05*life[None,:]),-700,700)))
        spectral=self.rgb_per_Lbol(T)
        self.rgb_age=(alive*(counts*L)[None,:])@spectral
        self.L_age=alive@(counts*L)
        # Blackbody photon count per emitted bolometric erg. x integral as
        # a rapidly convergent exponential series; no arbitrary HII color.
        x0=13.6*1.602176634e-12/(1.380649e-16*T)
        integral=np.zeros_like(x0)
        for n in range(1,80):
            integral+=np.exp(-n*x0)*(x0*x0/n+2*x0/n**2+2/n**3)
        q_per_erg=15/(np.pi**4*1.380649e-16*T)*integral
        self.Q_age=alive@(counts*L*LSUN*q_per_erg)

    def rgb_per_Lbol(self,T):
        T=np.asarray(T)
        return color(T,self.lut)/(SIGMA*T[...,None]**4/np.pi)/self.norm

class Populations2(Populations):
    def __init__(self):
        super().__init__()
        # Convolve the SAME instantaneous IMF-integrated SSP with a normalized
        # 2-Myr association formation history. Protostellar association mass
        # exists before all its optical emission; no random flashing source.
        raw_rgb=self.rgb_age.copy();raw_Q=self.Q_age.copy();raw_L=self.L_age.copy()
        nodes,weights=np.polynomial.legendre.leggauss(12)
        age=self.ages;duration=np.minimum(age,2.)
        birth_delay=.5*duration[:,None]*(nodes[None,:]+1)
        aa=age[:,None]-birth_delay; ww=.5*duration[:,None]*weights[None,:]/2
        self.rgb_age=np.stack([np.sum(np.interp(aa,age,raw_rgb[:,j])*ww,axis=1) for j in range(3)],axis=-1)
        self.Q_age=np.sum(np.interp(aa,age,raw_Q)*ww,axis=1)
        self.L_age=np.sum(np.interp(aa,age,raw_L)*ww,axis=1)
        m=np.geomspace(.08,100.,600)
        imf=np.where(m<.5,m**(-1.3),.5*m**(-2.3));count=imf*np.gradient(m);count/=np.sum(count*m)
        lum=np.where(m<.43,.23*m**2.3,np.where(m<2,m**4,np.where(m<20,1.5*m**3.5,1.5*20**3.5*(m/20)**1.8)))
        radius=np.where(m<1,m**.8,m**.6);T=5772*(lum/radius**2)**.25
        life=np.maximum(3.,1e4*m/lum*(1+.8*smooth((m-2)/8)))
        self.cohort_age=np.linspace(1200,12000,109)
        alive=1/(1+np.exp(np.clip((self.cohort_age[:,None]-life)/(life*.05),-700,700)))
        spec=(alive*(count*lum)[None,:])@self.rgb_per_Lbol(T)
        power=alive@(count*lum)
        self.ms_color=spec/power[:,None]
        self.old_giants=.75*self.rgb_per_Lbol(4600.)+.25*self.rgb_per_Lbol(3800.)
        self.mid_giants=.8*self.rgb_per_Lbol(4900.)+.2*self.rgb_per_Lbol(4100.)
        atomic=np.load(ROOT/'data/atomic/pyneb-v3.npz')
        from scipy.interpolate import RegularGridInterpolator
        self.line_ne=np.geomspace(10,200,150)
        self.line_coeff={k:np.exp(RegularGridInterpolator(
            (np.log(atomic['T']),np.log(atomic['ne'])),np.log(atomic[k]),bounds_error=True)
            (np.c_[np.full(150,np.log(10000.)),np.log(self.line_ne)])) for k in ('ha','hb')}
        self.filling_factor=.01

    def young(self,age,mass):
        age=np.asarray(age);mass=np.asarray(mass);a=np.clip(age,0,1000)
        rgb=np.stack([np.interp(a,self.ages,self.rgb_age[:,j]) for j in range(3)],axis=-1)*mass[...,None]
        return rgb,np.interp(a,self.ages,self.Q_age)*mass,np.interp(a,self.ages,self.L_age)*mass

    def cohort(self,cohort,t):
        age=(8000 if cohort==0 else 2000)+(t-400)
        ms=np.array([np.interp(age,self.cohort_age,self.ms_color[:,j]) for j in range(3)])
        fraction=.55 if cohort==0 else .35
        rgb=(1-fraction)*ms+fraction*(self.old_giants if cohort==0 else self.mid_giants)
        ml=(3.2*(age/8000)**.65 if cohort==0 else 1.6*(age/2000)**.65)
        return rgb/ml,ml,age

    def gas_limited_lines(self,Q,gas_column_Msun_kpc2,available_mass=None,cell_ids=None):
        # This is a subgrid ionized-clump closure, not a gas-temperature solve.
        # Adopt clumps 4x the local mean disk density, bounded to 10–200 cm^-3,
        # and an observationally plausible, prescribed volume filling factor.
        sigma=np.maximum(gas_column_Msun_kpc2,0)
        mean_n=sigma*1.98847e33/(3.085677581491367e21)**2/(np.sqrt(2*np.pi)*.08*3.085677581491367e21*1.4*1.6735575e-24)
        ne=np.clip(4*mean_n,10,200)
        if available_mass is None:available_mass=sigma*np.pi*.1**2*.3
        max_recomb=self.alpha_B*ne*available_mass*1.98847e33/(1.4*1.6735575e-24)
        # The geometrical region cannot exceed 150 pc. Extra photons are not
        # reassigned as optical line power. No minimum displayed physical size.
        geometry_cap=4*np.pi/3*(150*3.085677581491367e18)**3*self.alpha_B*ne**2*self.filling_factor
        used=np.minimum(.65*Q,geometry_cap)
        if cell_ids is None:used=np.minimum(used,max_recomb)
        else:
            # Share one available gas budget among all associations in a cell;
            # repeated formation samples must not each ionize the same mass.
            demanded_mass=used/self.alpha_B/ne*(1.4*1.6735575e-24)/1.98847e33
            totals=np.bincount(cell_ids,weights=demanded_mass,minlength=int(np.max(cell_ids))+1)
            used*=np.minimum(1,available_mass/np.maximum(totals[cell_ids],1e-90))
        lum={key:used*np.interp(ne,self.line_ne,coefficient)/self.alpha_B for key,coefficient in self.line_coeff.items()}
        rgb=sum(lum[k][...,None]/LSUN*self.line_rgb[k] for k in lum)
        radius=(3*used/(4*np.pi*self.alpha_B*ne**2*self.filling_factor))**(1/3)/3.085677581491367e18
        self.last_ne=ne;self.last_ionized_mass=used/self.alpha_B/ne*(1.4*1.6735575e-24)/1.98847e33
        return rgb,lum,radius,used

CACHE=ROOT/'data/optical.npz'

class FrozenPop:
    young=Populations2.young
    cohort=Populations2.cohort
    gas_limited_lines=Populations2.gas_limited_lines
    def __init__(self):
        if not CACHE.exists():raise RuntimeError('Required data/optical.npz is missing')
        with np.load(CACHE) as a:
            for k in ['ages','rgb_age','Q_age','L_age','cohort_age','ms_color','old_giants','mid_giants','line_ne','bulge_rgb','solar_rgb','sample_colors']:
                setattr(self,k,a[k].copy())
            self.line_rgb={k:a['line_rgb_'+k].copy() for k in ['ha','hb']}
            self.line_coeff={k:a['line_coeff_'+k].copy() for k in ['ha','hb']}
            self.alpha_B=float(a['alpha_B']);self.filling_factor=float(a['filling_factor'])
        self.report=CONFIG['presentation']['radiation']


def verify_tables():
    """Compare rebuilt small lookups without overwriting accepted data."""
    p=Populations2();f=FrozenPop()
    keys=['ages','rgb_age','Q_age','L_age','cohort_age','ms_color','old_giants','mid_giants','line_ne','bulge_rgb']
    errors={k:float(np.max(np.abs(getattr(p,k)-getattr(f,k)))/max(np.max(np.abs(getattr(f,k))),1e-90)) for k in keys}
    print(json.dumps(errors,indent=2))
    return errors

if __name__=='__main__':
    verify_tables()
