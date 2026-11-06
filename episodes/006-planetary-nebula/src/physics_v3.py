"""Conservative advected AGB cohorts and atomic nebular line emission.

Radial pressure evolution is calculated on a coarse angular grid. Subsector
launch inhomogeneities are passively carried on interpolated radial trajectories;
their extra lateral forces, ablation and nonlinear instabilities are NOT solved.
No emissivity texture, extra mass, fake knots or fixed contour curves are added.
"""
from pathlib import Path
import numpy as np
from numba import njit
from .physics_v2 import Nebula as PriorNebula
from .evolution_v3 import Track
from .evolution import PC,MSUN,MH,RSUN,LSUN,YR,KM_S,smooth
from .atomic_v3 import Radiation
from .ion_history_v3 import calculate as calculate_ions,filling as ion_filling
from .interface_v3 import Interface
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]

@njit(cache=True)
def remap_faces(faces,extensive,edges):
    """Exact conservative overlap of piecewise constant material-cell fields.

    Input extensive values already include their material volumes. Interpolation
    in r^3 preserves the spherical volume measure, not an image-space kernel.
    """
    n,nf=faces.shape;nr=len(edges)-1;nk=extensive.shape[-1]
    out=np.zeros((n,nr,nk))
    for a in range(n):
        j=0
        for k in range(nr):
            lo=edges[k];hi=edges[k+1]
            while j<nf-1 and faces[a,j+1]<=lo:j+=1
            jj=j
            while jj<nf-1 and faces[a,jj]<hi:
                l=max(lo,faces[a,jj]);h=min(hi,faces[a,jj+1])
                if h>l:
                    frac=(h**3-l**3)/(faces[a,jj+1]**3-faces[a,jj]**3)
                    for b in range(nk):out[a,k,b]+=extensive[a,jj,b]*frac
                jj+=1
    return out

class Nebula(PriorNebula):
    def __init__(self,nmu=48,nphi=96,nr=384,dynamics='radial-v3t-12x24-144.npz'):
        # Existing source/wind startup model remains authoritative before100yr.
        super().__init__(nmu=nmu,nphi=nphi,nr=nr)
        self.track=Track();self.radiation=Radiation(self.track)
        self.cache=np.load(ROOT/'cache'/dynamics)
        self.dynamic_t=self.cache['t'];self.dynamic_r=self.cache['r'];self.dynamic_v=self.cache['v'];self.dynamic_energy=self.cache['energy']
        name=dynamics.split('-')[2];cu,cp=map(int,name.split('x'))
        self.coarse_shape=(cu,cp)
        fu=(self.mu.ravel()+1)*cu/2-.5;fp=self.phi.ravel()*cp/(2*np.pi)-.5
        iu=np.floor(fu).astype(int);ip=np.floor(fp).astype(int);du=fu-iu;dp=fp-ip
        self.angular=[]
        for di,dj in [(0,0),(0,1),(1,0),(1,1)]:
            idx=np.clip(iu+di,0,cu-1)*cp+(ip+dj)%cp
            weight=(du if di else 1-du)*(dp if dj else 1-dp)
            self.angular.append((idx,weight))
        cells=self.cache['m'].shape[-1]
        self.cohort_ages=np.r_[0.,np.geomspace(30,300000,cells)]
        labels=np.broadcast_to(-self.cohort_ages,(self.mu.size,cells+1))
        self.material_mass=np.maximum(0,-np.diff(self.cohort(labels),axis=-1))
        self.cold_mass=float(self.material_mass.sum())
        if 'v3t' not in dynamics:raise ValueError('quasi-equilibrium runs failed finite-photon validation; usev3t')
        self.interface=Interface(self)
        fingerprint=hashlib.sha256((ROOT/'src/ion_history_v3.py').read_bytes()+(ROOT/'src/interface_v3.py').read_bytes()+(ROOT/'src/physics_v3.py').read_bytes()+(ROOT/'cache'/dynamics).read_bytes()+self.material_mass.tobytes()).hexdigest()
        photo=ROOT/'cache'/f'photo-{Path(dynamics).stem}-{nmu}x{nphi}-{fingerprint[:12]}.npz'
        if photo.exists():a=np.load(photo)
        else:
            begin=time.monotonic();tt,hh,bb=calculate_ions(self)
            closure=np.abs(bb[:,0]-bb[:,1:].sum(axis=-1)-hh.sum(axis=-1))/np.maximum(bb[:,0],1)
            assert closure.max()<1e-10
            np.savez_compressed(photo,t=tt,ions=hh,budget=bb)
            (ROOT/'validation'/f'{photo.stem}.json').write_text(json.dumps(dict(elapsed_s=time.monotonic()-begin,steps=len(tt)-1,max_photon_integral_relative_error=float(closure.max()),fingerprint=fingerprint,scope='direct grey dusty screen,12percent He allowance,sharp Case-B H front; no diffuse radiation'),indent=2)+'\n')
            a=np.load(photo)
        self.photo_path=photo;self.photo_t=a['t'];self.photo_ions=a['ions'];self.photo_budget=a['budget']
    def base_material(self,t):
        j=int(np.clip(np.searchsorted(self.dynamic_t,t)-1,0,len(self.dynamic_t)-2))
        a=float(np.clip((t-self.dynamic_t[j])/(self.dynamic_t[j+1]-self.dynamic_t[j]),0,1))
        r=(1-a)*self.dynamic_r[j].astype(float)+a*self.dynamic_r[j+1]
        v=(1-a)*self.dynamic_v[j].astype(float)+a*self.dynamic_v[j+1]
        rf=sum(r[idx]*weight[:,None] for idx,weight in self.angular)
        vf=sum(v[idx]*weight[:,None] for idx,weight in self.angular)
        return rf,vf
    def material(self,t,velocities=True):
        rf,vf=self.base_material(t);offset=self.interface.offset(t,rf)
        if velocities:
            h=.5
            ahead=self.interface.offset(t+h,rf+h*vf)
            behind=self.interface.offset(t-h,rf-h*vf)
            vf=vf+(ahead-behind)/(2*h)
        rf=rf+offset
        assert np.all(np.diff(rf,axis=-1)>0)
        return rf,vf
    def state(self,t):
        if t<100:return super().state(t)
        s=self.track.state(t);faces,velocity=self.material(t);n,nf=faces.shape;nc=nf-1
        rad=(.5*(faces[:,:-1]**3+faces[:,1:]**3))**(1/3)
        dr=np.diff(faces,axis=-1);vol=self.omega*(faces[:,1:]**3-faces[:,:-1]**3)/3*PC**3
        mass=self.material_mass;rho=mass*MSUN/vol;nH=rho*self.track.X/MH
        Td=(s['L']/(16*np.pi*5.670374419e-5*(rad*PC)**2))**.25
        kap=100*smooth((1400-Td)/300)
        shp=(*self.mu.shape,nc)
        j=int(np.clip(np.searchsorted(self.photo_t,t)-1,0,len(self.photo_t)-2))
        mix=float(np.clip((t-self.photo_t[j])/(self.photo_t[j+1]-self.photo_t[j]),0,1))
        ions=(1-mix)*self.photo_ions[j]+mix*self.photo_ions[j+1]
        pb=(1-mix)*self.photo_budget[j]+mix*self.photo_budget[j+1]
        derivative=(self.photo_ions[j+1].sum()-self.photo_ions[j].sum())/((self.photo_t[j+1]-self.photo_t[j])*YR)
        F=ion_filling(mass*MSUN*self.track.X/MH,ions)
        line,x,ne,rr,of,nfion=self.radiation.solve_timed(s,rad.reshape(shp),vol.reshape(shp),nH.reshape(shp),rho.reshape(shp),kap.reshape(shp),dr.reshape(shp),self.omega,F.reshape(shp))
        rr.update(photon_integral_relative_closure=float(abs(pb[0]-pb[1:].sum()-ions.sum())/max(pb[0],1)),stored_H_atoms=float(ions.sum()),emitted_ionizing_photons_since100yr=float(pb[0]),independent_H_storage_derivative_s=float(derivative),instantaneous_storage_rate_disagreement_over_Q=float(abs(rr['Q_net_ionizations_s']-derivative)/max(rr['Q_H_s'],1)))
        # Fixed H-alpha red, H-beta blue, OIII cyan, NII amber narrowband composite.
        # NII/Halpha wavelengths are both red in true color; their separation here
        # is explicitly false color. Gains do not alter luminosity diagnostics.
        colors=np.array([[1,.012,.004],[.03,.12,1],[.015,.80,1],[1,.25,.015]])
        colors/=(colors@np.array([.2126729,.7151522,.0721750]))[:,None]
        colors*=np.array([1.,.35,1.2,1.6])[:,None]
        line_rgb=(line@colors).reshape(n,nc,3)/(4*np.pi)/self.track.Yref
        rgb=line_rgb.copy()
        dtau=kap*rho*dr*PC;tau=np.cumsum(dtau,axis=-1);tin=tau-dtau
        bb=self.track.rgb(s['T'])*self.track.Yref;Lvis=4*np.pi**2*(s['R']*PC)**2*bb
        rgb+=kap[...,None]*rho[...,None]*.5*Lvis/(16*np.pi**2*(rad[...,None]*PC)**2)*np.exp(-tau[...,None])/self.track.Yref
        absorbed=.5*s['L']*self.omega/(4*np.pi)*np.exp(-tin)*(-np.expm1(-dtau))
        tout=.12*(tau[:,-1:]-tau+.5*dtau);thick=.12*tau[:,-1]>.67
        asec=absorbed.sum(axis=-1)
        T4=.75*asec[:,None]/(self.omega*5.670374419e-5*(rad*PC)**2)*(tout+2/3);Tgray=np.maximum(T4,1)**.25
        escape=np.divide(-np.expm1(-tout),tout,out=np.ones_like(tout),where=tout>1e-8)
        jbol=.12*kap*rho*5.670374419e-5*Tgray**4/np.pi
        escaping=4*np.pi*np.sum(jbol*escape*vol,axis=-1)
        jbol*=np.divide(asec,escaping,out=np.zeros_like(asec),where=escaping>0)[:,None]
        jbol=np.where(thick[:,None],jbol,absorbed/vol/(4*np.pi))
        Tir=np.where(thick[:,None],Tgray,Td)
        irfrac=np.interp(np.log(np.maximum(Tir,2)),np.log(self.ir_temperatures),self.ir_fractions)
        jir=jbol*irfrac/self.track.Yref
        physical_fields=np.concatenate([rgb,(kap*rho*PC)[...,None],jir[...,None],(.12*kap*rho*PC)[...,None]],axis=-1)
        # Remap integrated radiation/mass, not density-squared after averaging.
        # Thin material cells therefore retain their calculated line luminosity.
        edges=np.geomspace(faces[:,0].min()*.85,faces[:,-1].max()*1.001,self.nr+1)
        r=np.sqrt(edges[:-1]*edges[1:]);gridvol=self.omega*np.diff(edges**3)/3*PC**3
        extensive=np.concatenate([physical_fields*vol[...,None],mass[...,None],(4*np.pi*jir*self.track.Yref*escape*vol)[...,None],line_rgb*vol[...,None]],axis=-1)
        packed=remap_faces(faces,extensive,edges)
        field=(packed[...,:6]/gridvol[None,:,None]).reshape(*self.mu.shape,self.nr,6).astype('f4')
        escaped=packed[...,7].reshape(*self.mu.shape,self.nr)
        linefield=(packed[...,8:11]/gridvol[None,:,None]).reshape(*self.mu.shape,self.nr,3).astype('f4')
        gridmass=packed[...,6].sum();lost=float(self.track.massfn(self.track.rel[0])-self.track.massfn(t))
        x=x.reshape(n,nc);ne=ne.reshape(n,nc)
        energy=np.array([np.interp(t,self.dynamic_t,self.dynamic_energy[:,i]) for i in range(3)])
        istart=self.dynamic_energy[0,0]
        closure=abs(energy[0]-istart-energy[1]+energy[2])/max(istart+energy[1],1e-30)
        report=dict(t_relative_yr=float(t),T_eff_K=s['T'],L_Lsun=s['L']/LSUN,M_star_Msun=s['M'],R_Rsun=s['R']*PC/RSUN,
                    total_lost_Msun=lost,mass_grid_Msun=float(gridmass),outside_or_hot_reservoir_Msun=lost-float(gridmass),represented_cold_mass_expected_Msun=self.cold_mass,
                    wind_shell_mass_relative_error=abs(gridmass-self.cold_mass)/self.cold_mass,
                    line_luminosity_erg_s=float(sum(rr['line_luminosities_erg_s'].values())),
                    ionized_mass_Msun=float(np.sum(mass*x)),max_density_cm3=float(nH.max()),
                    shell_radius_pc_mean=float(faces[:,0].mean()),shell_speed_km_s_mean=float(velocity[:,0].mean()/KM_S),
                    radial_tau_V_max=float(tau.max()),recombination_time_yr_median=float(np.median(1/(2.59e-13*ne[ne>1e-3]*YR))) if np.any(ne>1e-3) else None,
                    dust_IR_tau_median=float(np.median(.12*tau[:,-1])),dust_IR_photosphere_pc_median=0.,dust_IR_thick_solid_angle_fraction=float(thick.mean()),
                    dust_absorbed_direct_luminosity_erg_s=float(absorbed.sum()),dust_NIR_luminosity_erg_s=float(escaped.sum()),
                    bubble_energy_relative_closure=float(closure),material_face_speed_km_s_p05_p50_p95=np.quantile(velocity/KM_S,[.05,.5,.95]).tolist(),
                    min_material_cell_width_pc=float(dr.min()),coarse_pressure_grid=list(self.coarse_shape),passive_subsector_grid=[self.nmu,self.nphi],
                    morphology_scope='radial pressure/shock evolution, passive subsector loading, bounded RT-informed contact closure; no nonlinear3D instability or Helix knots',**rr)
        report['interface_model']=self.interface.report
        result=dict(star=s,r=r,edges=edges,fields=field,line_fields=linefield,report=report,escaped_ir_power=escaped,IR_photosphere=np.zeros_like(self.mu),
                    material_radii=rad,material_x=x,material_O2=of,material_N1=nfion,material_line=line,material_mass=mass,material_vol=vol)
        return result
