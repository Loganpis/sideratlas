"""Low-cost proof of concept: conservative wind + one-zone hot bubble/sector shell.

Asymmetry, clumping, wind speeds, shell thickness and Te are prescribed. No 3D
fluid instabilities or complete photoionization/thermal balance are calculated.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ndtr
from .evolution import Track,MSUN,PC,YR,MH,KM_S,LSUN,KB,smooth

class Nebula:
    def __init__(self,nmu=32,nphi=64,nr=192):
        self.ir_temperatures=np.geomspace(2,100000,600)
        xlo=1.438776877/(40e-4*self.ir_temperatures);xhi=1.438776877/(2e-4*self.ir_temperatures)
        xx=xlo[:,None]+(xhi-xlo)[:,None]*np.linspace(0,1,256)
        self.ir_fractions=15/np.pi**4*np.trapz(xx**3*np.exp(-xx)/(-np.expm1(-xx)),xx,axis=-1)
        self.track=Track();self.nmu=nmu;self.nphi=nphi;self.nr=nr
        mu=-1+(np.arange(nmu)+.5)*2/nmu;phi=(np.arange(nphi)+.5)*2*np.pi/nphi
        self.mu,self.phi=np.meshgrid(mu,phi,indexing='ij');self.omega=4*np.pi/(nmu*nphi)
        # Assumed mild, even-parity equatorial enhancement; no unmodeled companion.
        self.A=1+.6*(1-self.mu**2)+.12*np.cos(2*self.phi)*(1-self.mu**2)
        self.A/=self.A.mean();self.vslow=12*KM_S*self.A**(-.12)
        # Fixed initial-condition clumping. Coherent launch-time variation, not
        # fluid turbulence or decorative emissivity added after projection.
        self.k=np.array([[3,4,.8],[7,12,1.9],[13,20,4.1],[21,32,2.8]])
        self.r0=4e-6 # pc ~0.8 AU, below final AGB radii; dust radii handled later
        self.birth_edges=np.unique(np.r_[-np.geomspace(300000,.1,440),self.track.rel[(self.track.rel>=-300000)&(self.track.rel<=0)],0.])
        self.birth_mid=(self.birth_edges[:-1]+self.birth_edges[1:])*.5
        self.source_f=self.clump(self.birth_mid).reshape(self.mu.size,-1)
        dm=-np.diff(self.track.massfn(self.birth_edges))
        self.source_C=np.c_[np.zeros(self.mu.size),np.cumsum(self.source_f*dm[None,:]*self.omega/(4*np.pi),axis=1)]
        self.shell_solution=self.solve_shell()
    def clump(self,birth):
        # Coherent angular correlations; each radial cohort normalized to preserve mass.
        # Stellar recoil/lateral momentum exchange are not evolved in this model.
        val=np.zeros((*self.mu.shape,len(birth)))
        for l,m,p in self.k:
            val+=np.cos(l*np.arccos(self.mu)[...,None]+p+np.asarray(birth)/4300)*np.cos(m*self.phi[...,None]+p)*.45
        f=self.A[...,None]*np.exp(val);return f/f.mean(axis=(0,1),keepdims=True)
    def cohort(self,t,rate=False):
        """Conservative per-solid-angle cumulative launched mass, same for force
        balance, free wind and swept shell. Clumps retain their material labels."""
        shape=np.asarray(t).shape;t=np.asarray(t);n=self.mu.size
        if t.shape[0]!=n:raise ValueError('one row per angular sector required')
        j=np.clip(np.searchsorted(self.birth_edges,t,side='right')-1,0,len(self.birth_mid)-1)
        row=np.arange(n).reshape((n,)+(1,)*(t.ndim-1))
        if rate:
            f=self.source_f[row,j]
            f=np.where(t<self.birth_edges[0],self.A.ravel().reshape((n,)+(1,)*(t.ndim-1)),f)
            return np.maximum(0,-self.track.mdotfn(t))*f*self.omega/(4*np.pi)
        dm=self.track.massfn(self.birth_edges[j])-self.track.massfn(np.clip(t,self.birth_edges[0],0))
        result=self.source_C[row,j]+self.source_f[row,j]*dm*self.omega/(4*np.pi)
        before=(self.track.massfn(self.birth_edges[0])-self.track.massfn(t))*self.A.ravel().reshape((n,)+(1,)*(t.ndim-1))*self.omega/(4*np.pi)
        return np.where(t<self.birth_edges[0],before,result)
    def fastwind(self,t):
        s=self.track.state(t);gate=float(smooth((s['T']-8000)/12000))
        mdot=max(0,-float(self.track.mdotfn(t)))*gate
        escape=np.sqrt(2*6.67430e-8*s['M']*MSUN/(s['R']*PC))/1e5
        v=min(2000.,max(150.,1.3*escape))*KM_S
        return mdot,v
    def solve_shell(self):
        # Reduced independent radial sectors share one hot-bubble pressure.
        # Fully radiative swept shell; adiabatic hot bubble, no conductive loss.
        # Conservative shell mass is obtained directly from the track's lost mass.
        shp=self.mu.shape;n=self.mu.size;f=self.A.ravel();vs=self.vslow.ravel();om=self.omega
        R0=np.full(n,4e-5)
        C0=self.cohort(np.zeros(n))
        m0=np.maximum(0,C0-self.cohort(-R0/vs))
        def fun(t,y):
            R=np.maximum(y[:n],1e-8);v=y[n:2*n];E=max(y[2*n],0)
            tr=np.minimum(0,t-R/vs)
            mass=np.maximum(m0,C0-self.cohort(tr))
            sweep=self.cohort(tr,rate=True)*np.maximum(0,v/vs-1)
            vol=np.sum(om*R**3/3);P=2*E/(3*vol)
            acc=(P*om*R**2-(v-vs)*sweep)/np.maximum(mass,1e-15)
            md,vf=self.fastwind(t)
            power=.5*md*vf*vf
            work=P*np.sum(om*R**2*v)
            return np.r_[v,acc,power-work,power,work]
        # E, integrated wind input and integrated work are separate diagnostics.
        sol=solve_ivp(fun,(0,getattr(self,"simulation_end_years",16500)),np.r_[R0,vs,m0.sum()*0,0,0],rtol=2e-5,atol=1e-13,max_step=35,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        self.shell_nfev=sol.nfev
        return sol
    def state(self,t):
        s=self.track.state(t)
        # Log radial cells, entire outer reservoir tracked separately. Inner
        # boundary follows physical stellar radius/dust launch; no fake material.
        rmin=max(1.02*s['R'],1e-7)
        if hasattr(self,"radial_inner_radius_pc"):rmin=max(rmin,self.radial_inner_radius_pc(t))
        rmax=max(.65,12*KM_S*max(t+100000,100000))
        if hasattr(self,"radial_outer_radius_pc"):rmax=max(rmax,self.radial_outer_radius_pc(t))
        edges=np.geomspace(rmin,rmax,self.nr+1);r=np.sqrt(edges[:-1]*edges[1:]);dr=np.diff(edges)
        nr=self.nr;n=self.mu.size
        if t>0:
            sh=self.shell_solution.sol(min(t,getattr(self,"simulation_end_years",16500)));Rsh=sh[:n].reshape(self.mu.shape);Vsh=sh[n:2*n].reshape(self.mu.shape)
        else:Rsh=np.zeros_like(self.mu);Vsh=self.vslow
        # Radial wind already at terminal velocity outside dust formation zone.
        born_lo=np.minimum(0,t-edges[:-1][None,None,:]/self.vslow[...,None])
        born_hi=np.minimum(0,t-edges[1:][None,None,:]/self.vslow[...,None])
        endmass=self.track.massfn(min(0,t))
        cell_lo=np.maximum(edges[:-1][None,None,:],Rsh[...,None])
        b_lo=np.minimum(min(0,t),t-cell_lo/self.vslow[...,None])
        b_hi=np.minimum(min(0,t),t-edges[1:][None,None,:]/self.vslow[...,None])
        mass=(self.cohort(b_lo.reshape(n,nr))-self.cohort(b_hi.reshape(n,nr))).reshape(*self.mu.shape,nr)
        mass=np.maximum(0,mass)
        birthmid=t-r/float(self.vslow.mean())
        f=self.clump(birthmid)
        swept=np.zeros_like(self.mu)
        if t>0:
            tr=np.minimum(0,t-Rsh/self.vslow)
            swept=(self.cohort(np.zeros(n))-self.cohort(tr.ravel())).reshape(self.mu.shape)
            # Finite compressed shell, no post-render contour surfaces. Width
            # and unresolved internal clumping are assumed, not RHD predictions.
            width=np.maximum(.075*Rsh,1e-6)
            cdf=ndtr((edges[None,None,:]-Rsh[...,None])/width[...,None])
            frac=np.diff(cdf,axis=-1);frac/=np.maximum(frac.sum(axis=-1,keepdims=True),1e-30)
            # Keep same conserved swept mass; inherited clump density retained.
            frac*=f;frac/=np.maximum(frac.sum(axis=-1,keepdims=True),1e-30)
            mass+=self.compressed_shell_mass(edges,Rsh,t) if hasattr(self,"compressed_shell_mass") else swept[...,None]*frac
        vol=self.omega*(edges[1:]**3-edges[:-1]**3)/3*PC**3
        rho=mass*MSUN/vol[None,None,:]
        nH=rho*self.track.X/MH
        # Dust in radiative equilibrium; opacity per TOTAL gas mass, not dust mass.
        Tdust=(s['L']/(16*np.pi*5.670374419e-5*(r*PC)**2))**.25
        dust=smooth((1400-Tdust)/300)
        kap=100*dust # assumed cm2/g; includes implicit dust/gas fraction
        tau_rad=np.cumsum(kap[None,None,:]*rho*dr[None,None,:]*PC,axis=-1)
        Q,Lion=self.track.ionizing(t);Q35,_=self.track.ionizing(t,35.12);Q55,_=self.track.ionizing(t,54.94);Q30,_=self.track.ionizing(t,29.6)
        Te=10000.;alpha=2.59e-13
        # Photon-conserving, radial Case-B equilibrium. Diffuse field/on-the-spot,
        # metal competition, time-dependent fronts and thermal balance omitted.
        budget=np.full(self.mu.shape,Q*self.omega/(4*np.pi));x=np.zeros_like(nH)
        qdust=0.;qrec=0.
        for j in range(nr):
            prev=budget.copy();budget*=np.exp(-3*kap[j]*rho[...,j]*dr[j]*PC);qdust+=np.sum(prev-budget)
            need=alpha*1.1*nH[...,j]**2*vol[j]
            used=np.minimum(budget,need);x[...,j]=np.sqrt(np.divide(used,need,out=np.zeros_like(used),where=need>0));budget-=used;qrec+=used.sum()
        ne=1.1*x*nH
        # Attenuated threshold photon rates; common attenuation is an approximate
        # ionic closure, not a multi-frequency ionization/charge-exchange solver.
        atten=np.exp(-3*tau_rad)*np.clip(1-np.cumsum(alpha*1.1*nH*nH*x*x*vol[None,None,:],axis=-1)/(Q*self.omega/(4*np.pi)+1e-100),0,1)
        area=4*np.pi*(r*PC)**2
        g35=Q35/area*1.5e-18*atten;g55=Q55/area*1e-18*atten
        ratio1=g35/(1e-11*ne+1e-99);ratio2=g55/(1e-11*ne+1e-99)
        ofrac=ratio1/(1+ratio1+ratio1*ratio2)
        nfrac=1/(1+Q30/area*1e-18*atten/(3e-12*ne+1e-99))
        ha=3.56e-25*ne*nH*x;hb=1.24e-25*ne*nH*x
        # Effective-collision-strength low-order proxies at fixed 10kK, with
        # collisional de-excitation. Exact atomic data are a subsequent upgrade.
        o3=5.5e-21*ne*nH*self.track.O_H*ofrac*x/(1+ne/7e5)
        n2=3.0e-21*ne*nH*self.track.N_H*nfrac*x/(1+ne/8e4)
        labs=Lion*qrec/max(Q,1e-99)
        # Metal-line emission cannot exceed the available photoelectron heat.
        heat=max(0,labs-13.6*1.602176634e-12*qrec)
        metal=float(np.sum((o3+n2)*vol));mcap=min(1,heat/max(metal,1e-99));o3*=mcap;n2*=mcap
        line=np.stack([ha,hb,o3,n2],axis=-1)
        llines=float(np.sum(line*vol[None,None,:,None]));ecap=min(1,labs/max(llines,1e-99));line*=ecap
        # Assigned narrowband display colors, unit linear luminance. Cyan OIII
        # and amber NII are disclosed colors, not measured temperature or RGB.
        colors=np.array([[1,.012,.004],[.03,.12,1],[.015,.80,1],[1,.20,.008]])
        colors/= (colors@np.array([.2126729,.7151522,.0721750]))[:,None]
        colors*=np.array([1.,.25,3.,1.])[:,None] # disclosed fixed narrowband display gains
        rgb=line@colors/(4*np.pi)/self.track.Yref
        # Single-scattering continuum, isotropic albedo .5, attenuated central
        # illumination; not multiple-scattering transfer or dust IR emission.
        bb=self.track.rgb(s['T'])*self.track.Yref
        Lvis=4*np.pi**2*(s['R']*PC)**2*bb
        scat=kap[None,None,:,None]*rho[...,None]*.5*Lvis/(16*np.pi**2*(r[None,None,:,None]*PC)**2)*np.exp(-tau_rad[...,None])/self.track.Yref
        rgb+=scat
        # Directly absorbed stellar power is re-emitted by equilibrium dust.
        # Only the calculated 1--5 micron blackbody fraction is displayed in a
        # separate false-color IR channel; cold outer dust is not made to glow.
        dtau=kap[None,None,:]*rho*dr[None,None,:]*PC
        tin=tau_rad-dtau
        absorbed_dust=.5*s['L']*self.omega/(4*np.pi)*np.exp(-tin)*(-np.expm1(-dtau))
        opacity=kap[None,None,:]*rho*PC # per pc
        ir_opacity=.12*kap[None,None,:]*rho # per cm
        tau_out=.12*(tau_rad[...,-1:]-tau_rad+.5*dtau)
        thick=.12*tau_rad[...,-1]>.67
        absorbed_sector=absorbed_dust.sum(axis=-1)
        # Grey Eddington diffusion/escape-probability closure for opaque dusty
        # winds. Temperature is set by the absorbed luminosity and optical depth;
        # no new energy source or removal of stellar extinction is introduced.
        T4=.75*absorbed_sector[...,None]/(self.omega*5.670374419e-5*(r*PC)**2)*(tau_out+2/3)
        Tgray=np.maximum(T4,1)**.25
        jbol=ir_opacity*5.670374419e-5*Tgray**4/np.pi
        escape=np.divide(-np.expm1(-tau_out),tau_out,out=np.ones_like(tau_out),where=tau_out>1e-8)
        escaping=4*np.pi*np.sum(jbol*escape*vol[None,None,:],axis=-1)
        jbol*=np.divide(absorbed_sector,escaping,out=np.zeros_like(escaping),where=escaping>0)[...,None]
        # Thin winds retain direct grain heating; thick winds redistribute the
        # absorbed energy by the declared grey closure. Band is 2--40 micron.
        jthin=absorbed_dust/vol[None,None,:]/(4*np.pi)
        jbol=np.where(thick[...,None],jbol,jthin)
        tempIR=np.where(thick[...,None],Tgray,Tdust[None,None,:])
        irfrac=np.interp(np.log(np.maximum(tempIR,2)),np.log(self.ir_temperatures),self.ir_fractions)
        jir=jbol*irfrac/self.track.Yref
        escaped_ir_power=4*np.pi*jir*self.track.Yref*escape*vol[None,None,:]
        Rir=np.zeros_like(self.mu)
        for iu in range(self.nmu):
            for ip in range(self.nphi):
                if thick[iu,ip]:Rir[iu,ip]=np.interp(2/3,tau_out[iu,ip,::-1],r[::-1])
        fields=np.concatenate([rgb,opacity[...,None],jir[...,None],(.12*opacity)[...,None]],axis=-1).astype('f4')
        top=np.minimum(min(t,0),t-edges[-1]/self.vslow.ravel())
        cutoff=np.full(n,min(t,0.)) if t>0 else t-edges[0]/self.vslow.ravel()
        expected_cold=float(np.sum(self.cohort(cutoff)-self.cohort(top)))
        mass_grid=float(mass.sum());lost=float(self.track.massfn(self.track.rel[0])-self.track.massfn(t))
        fastmass=max(0,self.track.M0-float(self.track.massfn(t))) if t>0 else 0.
        report=dict(t_relative_yr=t,T_eff_K=s['T'],L_Lsun=s['L']/LSUN,M_star_Msun=s['M'],R_Rsun=s['R']*PC/6.957e10,mass_grid_Msun=mass_grid,
            total_lost_Msun=lost,outside_or_hot_reservoir_Msun=lost-mass_grid,represented_cold_mass_expected_Msun=expected_cold,wind_shell_mass_relative_error=abs(mass_grid-expected_cold)/max(expected_cold,1e-20),Q_H_s=Q,Q_dust_s=qdust,Q_recomb_s=qrec,Q_escape_s=float(budget.sum()),
            photon_relative_closure=abs(Q-qdust-qrec-budget.sum())/max(Q,1),line_luminosity_erg_s=float(np.sum(line*vol[None,None,:,None])),absorbed_ionizing_luminosity_erg_s=labs,
            metal_heating_cap=mcap,line_energy_cap=ecap,dust_absorbed_direct_luminosity_erg_s=float(absorbed_dust.sum()),dust_NIR_luminosity_erg_s=float(escaped_ir_power.sum()),ionized_mass_Msun=float(np.sum(mass*x)),shell_radius_pc_mean=float(Rsh.mean()),shell_speed_km_s_mean=float(Vsh.mean()/KM_S),
            swept_mass_Msun=float(swept.sum()),fast_reservoir_mass_Msun=fastmass,max_density_cm3=float(nH.max()),radial_tau_V_max=float(tau_rad.max()),
            recombination_time_yr_median=float(np.median(1/(alpha*ne[ne>1e-3]*YR))) if np.any(ne>1e-3) else None)
        report['dust_IR_tau_median']=float(np.median(.12*tau_rad[...,-1]))
        report['dust_IR_photosphere_pc_median']=float(np.median(Rir))
        report['dust_IR_thick_solid_angle_fraction']=float(thick.mean())
        return dict(star=s,r=r,edges=edges,fields=fields,report=report,escaped_ir_power=escaped_ir_power,IR_photosphere=Rir)
