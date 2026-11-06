"""Bounded RT-informed contact reconstruction, NOT nonlinear3D hydrodynamics.

Seeds come from unresolved angular wind loading. The background acceleration
sets the classical local incompressible growth estimate sqrt(A*k*a). Positive
growth is integrated in physical length; nonlinear amplitude/evanescent spatial
extension are explicit reduced-order closures. No independent random texture.
"""
import numpy as np
from scipy.special import sph_harm

class Interface:
    def __init__(self,model):
        self.model=model;t=model.dynamic_t
        R=model.dynamic_r[:,:,0].mean(axis=1).astype(float)
        v=model.dynamic_v[:,:,0].mean(axis=1).astype(float)
        acceleration=np.gradient(v,t)
        m=model.cache['m'];r=model.dynamic_r.astype(float)
        volume=model.omega*(r[:,:,1]**3-r[:,:,0]**3)/3
        # Use the coarse solid angle, not the fine radiation-grid angle.
        coarse_omega=4*np.pi/np.prod(model.coarse_shape)
        cold=(m[:,0][None,:]/(volume*coarse_omega/model.omega)).mean(axis=1)
        hot=np.array([max(0,model.track.M0-model.track.state(float(x))['M']) for x in t])/(4*np.pi*R**3/3)
        atwood=np.clip((cold-hot)/np.maximum(cold+hot,1e-99),0,1)
        ell=[];basis=[];theta=np.arccos(model.mu.ravel());phi=model.phi.ravel()
        for l in range(12,33):
            for mm in range(-l,l+1):
                y=sph_harm(abs(mm),l,phi,theta)
                real=y.real if mm==0 else (np.sqrt(2)*y.real if mm>0 else np.sqrt(2)*y.imag)
                basis.append(real);ell.append(l)
        self.B=np.array(basis,dtype='f4').T;self.ell=np.array(ell)
        fine=model.material_mass[:,:4].sum(axis=-1);fine/=fine.mean()
        coarse=m[:,:4].sum(axis=-1);coarse/=coarse.mean()
        resolved=sum(coarse[idx]*weight for idx,weight in model.angular)
        seed=fine-resolved;seed-=seed.mean()
        # Local energy-driven wind-bubble response R∝density_normalization^(-1/3).
        # This seed rule and its extrapolation to a corrugated contact are assumed.
        self.c0=-(R[0]/3)*(self.B.T@seed)*model.omega
        growth=[]
        for l in range(12,33):
            rate=np.sqrt(np.maximum(acceleration,0)*atwood*np.sqrt(l*(l+1))/R)
            integral=np.r_[0,np.cumsum(.5*(rate[1:]+rate[:-1])*np.diff(t))]
            growth.append(integral)
        self.growth=np.array(growth);self.R=R;self.t=t;self.effective_l=24.
        self.report=dict(mechanism='RT-susceptible accelerated dense/light interface; local incompressible growth estimate',
            mode_degrees=[12,32],mode_count=len(ell),seed='fine AGB mass loading minus interpolated coarse loading; no new random field',
            growth_equation='gamma_l=sqrt(max(a,0)*Atwood*sqrt(l*(l+1))/R)',
            integrated_efolds_l12_l32_at16k=[float(np.interp(16000,t,self.growth[0])),float(np.interp(16000,t,self.growth[-1]))],
            atwood_range=[float(atwood.min()),float(atwood.max())],
            unresolved_closure='smoothly bound total contact displacement at k_eff*xi=0.6 (about0.1wavelength); exponential radial extension; conserve bubble volume',
            scope='prescribed nonlinear cap, no simulated nonlinear instability, mode coupling, mixing, knot formation, or lateral backreaction')
    def offset(self,t,faces):
        amp=np.array([np.exp(np.interp(t,self.t,g)) for g in self.growth])
        coeff=self.c0*amp[self.ell-12]
        displacement=self.B@coeff
        contact=faces[:,0];limit=.6*contact/self.effective_l
        delta=limit*np.tanh(displacement/np.maximum(limit,1e-30))
        # Preserve the background bubble volume, including second-order terms.
        desired=np.sum(contact**3)
        for _ in range(4):delta-=(np.sum((contact+delta)**3)-desired)/(3*np.sum((contact+delta)**2))
        shape=np.exp(-self.effective_l*np.maximum(0,faces-contact[:,None])/contact[:,None])
        return delta[:,None]*shape
