"""Inexpensive Lagrangian radial columns; not a 3D hydrodynamic simulation.

Fixed cold masses, pressure-gradient momentum equation, shared adiabatic wind
bubble and compressive von Neumann viscosity. Photoheated gas is isothermal at
10kK; neutral gas100K. The thermal thermostat is explicit, not total-energy RHD.
Grey Case-B H equilibrium supplies the pressure ionization closure. The more
detailed line-emission calculation is separate and audited against this closure.
"""
import numpy as np
from numba import njit
from .evolution import KM_S,MSUN,PC,MH,YR,KB

@njit(cache=True)
def evolve(r,v,m,omega,times,wind_time,wind_power,qh,star_mass,E0,progress_limit=200000):
    n,nf=r.shape;nc=nf-1
    rr=np.empty((len(times),n,nf),np.float32);vv=np.empty_like(rr)
    ee=np.zeros((len(times),5))
    P=np.zeros((n,nc));a=np.zeros_like(r);x=np.zeros_like(m)
    t=times[0];E=E0;inputE=0.;workE=0.;steps=0;dtmin=1e99;maxv=0.
    cs_neutral2=(.8*KM_S)**2;cs_ion2=(11.7*KM_S)**2
    Gcode=6.67430e-8*MSUN*YR*YR/PC**3
    for out in range(len(times)):
        while t<times[out]-1e-10:
            Q=np.interp(t,wind_time,qh)*.88 # He/metal photon competition approximation, pressure only
            power=np.interp(t,wind_time,wind_power)
            massstar=np.interp(t,wind_time,star_mass)
            vol=0.
            for k in range(n):vol+=omega*r[k,0]**3/3
            Pb=2*E/(3*max(vol,1e-30));dt=25.
            for k in range(n):
                budget=Q*omega/(4*np.pi)
                for j in range(nc):
                    dV=omega*(r[k,j+1]**3-r[k,j]**3)/3
                    rho=m[k,j]/max(dV,1e-40)
                    nH=rho*MSUN/PC**3*.7004434805295725/MH
                    need=2.59e-13*1.1*nH*nH*dV*PC**3
                    used=min(budget,need);xx=np.sqrt(used/max(need,1e-100));budget-=used
                    x[k,j]=xx
                    cs2=cs_neutral2+(cs_ion2-cs_neutral2)*xx
                    dv=v[k,j+1]-v[k,j]
                    # Shock viscosity only on compression, no arbitrary stirring.
                    visc=0. if dv>=0 else rho*(2.*dv*dv+.5*np.sqrt(cs2)*(-dv))
                    P[k,j]=rho*cs2+visc
                    dx=r[k,j+1]-r[k,j]
                    if dx<=0:raise ValueError('Lagrangian material cells crossed')
                    dt=min(dt,.18*dx/(np.sqrt(cs2)+2*abs(dv)+1e-30))
                for j in range(nf):
                    left=Pb if j==0 else P[k,j-1]
                    right=0. if j==nc else P[k,j]
                    mf=.5*m[k,0] if j==0 else (.5*m[k,-1] if j==nc else .5*(m[k,j-1]+m[k,j]))
                    a[k,j]=omega*r[k,j]**2*(left-right)/max(mf,1e-30)-Gcode*massstar/max(r[k,j]**2,1e-30)
                    # Resolved gravity is small; constrain early interface acceleration.
                    dt=min(dt,.15*np.sqrt(max(r[k,min(j+1,nc)]-r[k,max(j-1,0)],1e-20)/(abs(a[k,j])+1e-40)))
            dt=min(dt,times[out]-t)
            dtmin=min(dtmin,dt)
            for k in range(n):
                for j in range(nf):
                    r[k,j]+=v[k,j]*dt+.5*a[k,j]*dt*dt
                    v[k,j]+=a[k,j]*dt
                    maxv=max(maxv,abs(v[k,j])/KM_S)
            vnew=0.
            for k in range(n):vnew+=omega*r[k,0]**3/3
            # Exact adiabatic volume scaling of the one-zone bubble, trapezoid source.
            Enew=(E+.5*power*dt)*(vol/max(vnew,1e-40))**(2/3)+.5*power*dt
            inputE+=power*dt;workE+=E+power*dt-Enew;E=Enew
            t+=dt;steps+=1
            if steps>progress_limit:raise ValueError('radial-column cost guard exceeded')
        rr[out]=r;vv[out]=v
        ee[out]=np.array([E,inputE,workE,float(steps),dtmin])
    return rr,vv,ee,maxv

def initial_state(model,ncells=160,start=100.):
    """Conservative terminal-wind cohorts, inherited100yr startup shell state.

    The tiny startup rim is a mass-preserving ordered compression, not a resolved
    photospheric wind launch. This closure ends once pressure evolution starts.
    """
    n=model.mu.size;vs=model.vslow.ravel()
    ages=np.r_[0.,np.geomspace(30.,300000.,ncells)]
    born=np.broadcast_to(-ages,(n,len(ages)))
    cum=model.cohort(born);m=np.maximum(0,-np.diff(cum,axis=1))
    sh=model.shell_solution.sol(start);R=sh[:n];V=sh[n:2*n]
    radii=vs[:,None]*(start+ages[None,:])
    swept_age=np.maximum(0,R/vs-start)
    # Preserve radial ordering and every cohort; no overlapping Gaussian kernels.
    u=ages[None,:]/np.maximum(swept_age[:,None],1e-10)
    radii=np.where(u<1,R[:,None]*(.90+.10*u),radii)
    radii[:,0]=.90*R
    velocities=np.broadcast_to(vs[:,None],radii.shape).copy()
    mix=np.exp(-((radii-R[:,None])/(.12*R[:,None]))**2)
    velocities+=(V[:,None]-vs[:,None])*mix
    return radii,velocities,m,float(sh[2*n])

@njit(cache=True)
def evolve_time(r,v,m,omega,times,wind_time,wind_power,qh,star_mass,E0,ions,progress_limit=3000000):
    """Sharp-front Case-B time-dependent H, finite photon inventory.

    Ions is a per-column count, NOT a fraction of an arbitrary voxel. Inside the
    front H is ionized; the partial front cell has a volume filling fraction.
    Backward-Euler recombination is solved exactly piecewise in material mass.
    This remains direct-ray/on-the-spot, grey, no lateral transport or full RHD.
    """
    n,nf=r.shape;nc=nf-1
    rr=np.empty((len(times),n,nf),np.float32);vv=np.empty_like(rr)
    ii=np.empty((len(times),n));ee=np.zeros((len(times),8))
    P=np.zeros((n,nc));a=np.zeros_like(r);Qcol=np.zeros(n)
    NH=m*MSUN*.7004434805295725/MH
    t=times[0];E=E0;inputE=0.;workE=0.;steps=0;dtmin=1e99;maxv=0.
    photons=0.;recombined=0.;escaped=0.
    c0=(.8*KM_S)**2;c1=(11.7*KM_S)**2
    Gcode=6.67430e-8*MSUN*YR*YR/PC**3
    for out in range(len(times)):
        while t<times[out]-1e-10:
            Q=np.interp(t,wind_time,qh)*.88
            power=np.interp(t,wind_time,wind_power);massstar=np.interp(t,wind_time,star_mass)
            vol=0.
            for k in range(n):vol+=omega*r[k,0]**3/3
            Pb=2*E/(3*max(vol,1e-30));dt=25.
            for k in range(n):
                tau=0.
                for j in range(nc):
                    dV=omega*(r[k,j+1]**3-r[k,j]**3)/3*PC**3
                    tau+=100*m[k,j]*MSUN/max(dV,1e-50)*(r[k,j+1]-r[k,j])*PC
                Qray=Q*omega/(4*np.pi)*np.exp(-min(700,3*tau));Qcol[k]=Qray
                leftions=ions[k];foundfront=False;recprefix=0.
                for j in range(nc):
                    f=min(1.,max(0.,leftions/max(NH[k,j],1e-100)));leftions-=f*NH[k,j]
                    dV=omega*(r[k,j+1]**3-r[k,j]**3)/3
                    Rfull=2.59e-13*1.1*(NH[k,j]/max(dV*PC**3,1e-50))*NH[k,j]
                    if f<1 and not foundfront:
                        # Limit heating/front movement within the first neutral
                        # material cell; prevents a large pressure jump in one step.
                        dt=min(dt,.05*NH[k,j]/max(abs(Qray-recprefix-f*Rfull)*YR,1e-90));foundfront=True
                    recprefix+=f*Rfull
                    rho=m[k,j]/max(dV,1e-40);cs2=c0+(c1-c0)*f
                    dv=v[k,j+1]-v[k,j]
                    visc=0. if dv>=0 else rho*(2*dv*dv+.5*np.sqrt(cs2)*(-dv))
                    P[k,j]=rho*cs2+visc;dx=r[k,j+1]-r[k,j]
                    if dx<=0:raise ValueError('time-dependent material cells crossed')
                    dt=min(dt,.18*dx/(np.sqrt(cs2)+2*abs(dv)+1e-30))
                for j in range(nf):
                    left=Pb if j==0 else P[k,j-1];right=0. if j==nc else P[k,j]
                    mf=.5*m[k,0] if j==0 else (.5*m[k,-1] if j==nc else .5*(m[k,j-1]+m[k,j]))
                    a[k,j]=omega*r[k,j]**2*(left-right)/max(mf,1e-30)-Gcode*massstar/max(r[k,j]**2,1e-30)
                    dt=min(dt,.15*np.sqrt(max(r[k,min(j+1,nc)]-r[k,max(j-1,0)],1e-20)/(abs(a[k,j])+1e-40)))
            dt=min(dt,times[out]-t);dtmin=min(dtmin,dt)
            for k in range(n):
                for j in range(nf):
                    r[k,j]+=v[k,j]*dt+.5*a[k,j]*dt*dt;v[k,j]+=a[k,j]*dt
                    maxv=max(maxv,abs(v[k,j])/KM_S)
            vnew=0.
            for k in range(n):vnew+=omega*r[k,0]**3/3
            Enew=(E+.5*power*dt)*(vol/max(vnew,1e-40))**(2/3)+.5*power*dt
            inputE+=power*dt;workE+=E+power*dt-Enew;E=Enew
            for k in range(n):
                supply=ions[k]+Qcol[k]*dt*YR;newions=0.;rec=0.
                for j in range(nc):
                    dV=omega*(r[k,j+1]**3-r[k,j]**3)/3*PC**3
                    density=NH[k,j]/max(dV,1e-50)
                    R=2.59e-13*1.1*density*NH[k,j]
                    demand=NH[k,j]+dt*YR*R
                    f=min(1.,max(0.,supply/max(demand,1e-100)))
                    supply-=f*demand;newions+=f*NH[k,j];rec+=f*dt*YR*R
                ions[k]=newions;recombined+=rec;escaped+=max(0,supply)
            photons+=Qcol.sum()*dt*YR;t+=dt;steps+=1
            if steps>progress_limit:raise ValueError('time-dependent cost guard exceeded')
        rr[out]=r;vv[out]=v;ii[out]=ions
        ee[out]=np.array([E,inputE,workE,float(steps),dtmin,photons,recombined,escaped])
    return rr,vv,ee,maxv,ii
