"""Conservative isothermal gas: minmod/Rusanov flux and RK2."""
import numpy as np
from numba import njit
from .dynamics import KM_S, circular_speed
from .spiral import Gravity2 as GridGravity

@njit(cache=True)
def minmod(a,b,c):
    if a>0 and b>0 and c>0:
        return min(a,b,c)
    if a<0 and b<0 and c<0:
        return max(a,b,c)
    return 0.

@njit(cache=True)
def rhs(q, gravity, dx, cs):
    n=q.shape[1]
    w=np.empty((3,n+4,n+4)); sx=np.zeros_like(w); sy=np.zeros_like(w)
    for j in range(n+4):
        jj=min(n-1,max(0,j-2))
        for i in range(n+4):
            ii=min(n-1,max(0,i-2)); d=max(q[0,jj,ii],1e-90)
            w[0,j,i]=d; w[1,j,i]=q[1,jj,ii]/d; w[2,j,i]=q[2,jj,ii]/d
    for k in range(3):
        for j in range(1,n+3):
            for i in range(1,n+3):
                sx[k,j,i]=minmod(2*(w[k,j,i]-w[k,j,i-1]),
                    .5*(w[k,j,i+1]-w[k,j,i-1]),2*(w[k,j,i+1]-w[k,j,i]))
                sy[k,j,i]=minmod(2*(w[k,j,i]-w[k,j-1,i]),
                    .5*(w[k,j+1,i]-w[k,j-1,i]),2*(w[k,j+1,i]-w[k,j,i]))
    fx=np.empty((3,n,n+1)); fy=np.empty((3,n+1,n)); c2=cs*cs
    for j in range(n):
        jj=j+2
        for i in range(n+1):
            ii=i+1
            dl=w[0,jj,ii]+.5*sx[0,jj,ii]; dr=w[0,jj,ii+1]-.5*sx[0,jj,ii+1]
            ul=w[1,jj,ii]+.5*sx[1,jj,ii]; ur=w[1,jj,ii+1]-.5*sx[1,jj,ii+1]
            vl=w[2,jj,ii]+.5*sx[2,jj,ii]; vr=w[2,jj,ii+1]-.5*sx[2,jj,ii+1]
            speed=max(abs(ul),abs(ur))+cs
            fx[0,j,i]=.5*(dl*ul+dr*ur-speed*(dr-dl))
            fx[1,j,i]=.5*(dl*ul*ul+c2*dl+dr*ur*ur+c2*dr-speed*(dr*ur-dl*ul))
            fx[2,j,i]=.5*(dl*ul*vl+dr*ur*vr-speed*(dr*vr-dl*vl))
    for j in range(n+1):
        jj=j+1
        for i in range(n):
            ii=i+2
            dl=w[0,jj,ii]+.5*sy[0,jj,ii]; dr=w[0,jj+1,ii]-.5*sy[0,jj+1,ii]
            ul=w[1,jj,ii]+.5*sy[1,jj,ii]; ur=w[1,jj+1,ii]-.5*sy[1,jj+1,ii]
            vl=w[2,jj,ii]+.5*sy[2,jj,ii]; vr=w[2,jj+1,ii]-.5*sy[2,jj+1,ii]
            speed=max(abs(vl),abs(vr))+cs
            fy[0,j,i]=.5*(dl*vl+dr*vr-speed*(dr-dl))
            fy[1,j,i]=.5*(dl*vl*ul+dr*vr*ur-speed*(dr*ur-dl*ul))
            fy[2,j,i]=.5*(dl*vl*vl+c2*dl+dr*vr*vr+c2*dr-speed*(dr*vr-dl*vl))
    out=np.empty_like(q)
    for j in range(n):
        for i in range(n):
            for k in range(3):
                out[k,j,i]=-(fx[k,j,i+1]-fx[k,j,i]+fy[k,j+1,i]-fy[k,j,i])/dx
            out[1,j,i]+=q[0,j,i]*gravity[0,j,i]
            out[2,j,i]+=q[0,j,i]*gravity[1,j,i]
    # Outward boundary mass/angular-momentum flux, not an interior residual.
    mf=0.; jf=0.; edge=.5*n*dx
    for j in range(n):
        y=(j+.5)*dx-edge
        mf+=dx*(fx[0,j,n]-fx[0,j,0])
        jf+=dx*(edge*fx[2,j,n]-y*fx[1,j,n]
                -(-edge*fx[2,j,0]-y*fx[1,j,0]))
    for i in range(n):
        x=(i+.5)*dx-edge
        mf+=dx*(fy[0,n,i]-fy[0,0,i])
        jf+=dx*(x*fy[2,n,i]-edge*fy[1,n,i]
                -(x*fy[2,0,i]+edge*fy[1,0,i]))
    torque=np.sum(q[0]*gravity[4])*dx*dx
    power=np.sum(q[0]*gravity[3])*dx*dx
    return out,np.array([mf,jf,torque,power])

class GasDisk:
    def __init__(self,n=256,length=36.,forcing=True,gravity_class=GridGravity,cs_kms=7.,sf_model='original',tdep=2000.,clumping=2.):
        self.n=n; self.length=length; self.dx=length/n; self.cs=cs_kms*KM_S
        self.sf_model=sf_model;self.tdep=tdep;self.clumping=clumping
        a=(np.arange(n)+.5)*self.dx-length/2
        self.X,self.Y=np.meshgrid(a,a); R=np.hypot(self.X,self.Y)
        self.gravity=gravity_class(self.X,self.Y); self.forcing=forcing
        sigma=np.exp(-R/4.5)/(1+np.exp((R-15.)/.4))
        sigma*=6e9/(sigma.sum()*self.dx**2)
        self.floor=1e-10*sigma.max()
        vc=circular_speed(R.ravel()).reshape(R.shape)
        # Effective isothermal pressure support in the exponential interior.
        vphi=np.sqrt(np.maximum(vc*vc-self.cs**2*R/4.5,0))
        self.q=np.array([sigma,-sigma*vphi*self.Y/R,sigma*vphi*self.X/R])
        self.initial_mass=self.mass(); self.initial_J=self.angular_momentum()
        self.mass_boundary=0.; self.J_boundary=0.; self.J_external=0.
        self.explicit_work=0.; self.floor_mass=0.; self.formed_mass=0.; self.formed_J=0.
        self.mass_pending=np.zeros_like(sigma)
        self.momentum_pending=np.zeros((2,*sigma.shape))
        self.KE_pending=0.

    def molecular_fraction(self):
        if self.sf_model=='original':return self.q[0]/(self.q[0]+20e6)
        # Solar-metallicity two-phase equilibrium KMT closure, not full KMT+
        # UV/CNM iteration, chemistry, or resolved cloud self-shielding.
        tau=.066*self.clumping*np.maximum(self.q[0]/1e6,1e-30)
        chi=.77*(1+3.1)
        s=np.log(1+.6*chi+.01*chi*chi)/(.6*tau)
        return np.clip(1-.75*s/(1+.25*s),0,1)

    def mass(self):
        return float(self.q[0].sum()*self.dx**2)

    def angular_momentum(self):
        return float(np.sum(self.X*self.q[2]-self.Y*self.q[1])*self.dx**2)

    def energy(self,t):
        q=self.q; rho=q[0]; v2=(q[1]**2+q[2]**2)/rho**2
        phi=self.gravity.at(t,self.forcing)[2]
        return float(np.sum(rho*(.5*v2+self.cs**2*np.log(rho/1e6)+phi))*self.dx**2)

    def dt_limit(self):
        q=self.q
        return min(.125,.38*self.dx/np.max(np.abs(q[1]/q[0])+np.abs(q[2]/q[0])+2*self.cs))

    def step(self,t,dt,form_stars=True):
        g0=self.gravity.at(t,self.forcing); g1=self.gravity.at(t+dt,self.forcing)
        r0,b0=rhs(self.q,g0,self.dx,self.cs)
        q1=self.q+dt*r0
        if q1[0].min()<=0:
            raise RuntimeError('Negative RK intermediate gas density; reduce timestep')
        r1,b1=rhs(q1,g1,self.dx,self.cs)
        self.q=.5*(self.q+q1+dt*r1)
        if self.q[0].min()<=0 or not np.all(np.isfinite(self.q)):
            raise RuntimeError('Nonpositive/nonfinite gas state')
        b=.5*dt*(b0+b1)
        self.mass_boundary+=b[0]; self.J_boundary+=b[1]
        self.J_external+=b[2]; self.explicit_work+=b[3]
        # Explicitly audited floor; never used to conceal a negative state.
        low=self.q[0]<self.floor
        if np.any(low):
            self.floor_mass+=float(np.sum(self.floor-self.q[0,low])*self.dx**2)
            self.q[0,low]=self.floor
        if form_stars:
            mol=self.molecular_fraction()
            fraction=-np.expm1(-dt*mol/self.tdep)
            dm=self.q[0]*fraction*self.dx**2
            dj=float(np.sum((self.X*self.q[2]-self.Y*self.q[1])*fraction)*self.dx**2)
            self.mass_pending+=dm; self.formed_mass+=float(dm.sum()); self.formed_J+=dj
            self.momentum_pending+=self.q[1:]*fraction[None,:,:]*self.dx**2
            self.KE_pending+=float(np.sum(.5*(self.q[1]**2+self.q[2]**2)/self.q[0]*fraction)*self.dx**2)
            self.q*=1-fraction[None,:,:]

    def report(self,t):
        R=np.hypot(self.X,self.Y); q=self.q; rho=q[0]
        vr=(self.X*q[1]+self.Y*q[2])/(np.maximum(R,1e-12)*rho)/KM_S
        return dict(t_Myr=float(t),mass_Msun=self.mass(),formed_Msun=self.formed_mass,
            boundary_out_Msun=self.mass_boundary,floor_added_Msun=self.floor_mass,
            mass_relative_closure=(self.mass()+self.formed_mass+self.mass_boundary
                                   -self.floor_mass-self.initial_mass)/self.initial_mass,
            J_gas=self.angular_momentum(),J_formed=self.formed_J,J_external=self.J_external,
            J_boundary_out=self.J_boundary,
            J_relative_residual=(self.angular_momentum()+self.formed_J+self.J_boundary
                                -self.J_external-self.initial_J)/abs(self.initial_J),
            barotropic_energy=self.energy(t),explicit_potential_work=self.explicit_work,
            energy_scope='isothermal/barotropic diagnostic; shocks and closure dissipate energy; no total thermal/radiative energy solution',
            surface_density_Msun_pc2_p50_p95_p99=np.quantile(rho[R<13]/1e6,[.5,.95,.99]).tolist(),
            radial_velocity_km_s_mass_rms=float(np.sqrt(np.sum(rho*vr**2)/rho.sum())),
            effective_sound_speed_km_s=self.cs/KM_S,gas_cell_pc=1000*self.dx,
            sfr_Msun_yr=float(np.sum(rho*self.molecular_fraction())*self.dx**2/self.tdep/1e6))
