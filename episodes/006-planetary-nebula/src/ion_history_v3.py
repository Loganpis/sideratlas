"""Finite photon inventory on passively advected material, independent of camera.

Grey dusty screen +12% He photon allowance; sharp H front, Case-B on-the-spot,
10kK ionized material. This is a time-dependent reduced transfer approximation,
not a multi-element Cloudy calculation or diffuse3D transfer.
"""
from pathlib import Path
import numpy as np
from numba import njit
from .evolution import PC,MSUN,MH,YR

@njit(cache=True)
def advance(faces,NH,ions,omega,dt,q,he_fraction):
    n,nf=faces.shape;nc=nf-1;dust=0.;helium=0.;rec=0.;escape=0.
    emitted=q*dt*YR
    for k in range(n):
        tau=0.
        for j in range(nc):
            vol=omega*(faces[k,j+1]**3-faces[k,j]**3)/3*PC**3
            rho=NH[k,j]*MH/.7004434805295725/max(vol,1e-99)
            tau+=100*rho*(faces[k,j+1]-faces[k,j])*PC
        qray=q*omega/(4*np.pi);trans=np.exp(-min(700,3*tau))
        dust+=qray*(1-trans)*dt*YR;helium+=qray*trans*he_fraction*dt*YR
        supply=ions[k]+qray*trans*(1-he_fraction)*dt*YR;ni=0.
        for j in range(nc):
            vol=omega*(faces[k,j+1]**3-faces[k,j]**3)/3*PC**3
            R=2.59e-13*1.1*NH[k,j]**2/max(vol,1e-99)
            demand=NH[k,j]+dt*YR*R
            f=min(1.,max(0.,supply/max(demand,1e-99)))
            supply-=f*demand;ni+=f*NH[k,j];rec+=f*dt*YR*R
        ions[k]=ni;escape+=max(0,supply)
    return np.array([emitted,dust,helium,rec,escape])

def calculate(model,step_factor=1.):
    end=float(model.dynamic_t[-1]);ions=np.zeros(model.mu.size)
    NH=model.material_mass*MSUN*model.track.X/MH
    times=np.unique(np.r_[np.arange(100,min(1000,end),2*step_factor),np.arange(1000,min(10000,end),10*step_factor),np.arange(10000,end,50*step_factor),end])
    times=times[(times>=100)&(times<=end)];history=np.zeros((len(times),len(ions)));budget=np.zeros((len(times),5));cum=np.zeros(5)
    for j in range(1,len(times)):
        mid=.5*(times[j]+times[j-1]);faces,_=model.material(mid,velocities=False);q=model.track.ionizing(mid)[0]
        he=min(.12,model.track.ionizing(mid,24.59)[0]/max(q,1e-99))
        cum+=advance(faces,NH,ions,model.omega,times[j]-times[j-1],q,he)
        history[j]=ions;budget[j]=cum
    return times,history,budget

def filling(NH,ion_counts):
    remaining=ion_counts[:,None]-np.c_[np.zeros(len(NH)),np.cumsum(NH,axis=-1)[:,:-1]]
    return np.clip(remaining/np.maximum(NH,1e-100),0,1)
