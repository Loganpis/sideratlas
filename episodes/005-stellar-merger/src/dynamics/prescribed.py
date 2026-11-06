"""Small gravitational parcel model and energy-limited illustrative scenario.

There is no self-consistent binary inspiral, fluid pressure solver, resolved
shock or disruption calculation here. The existing MESA profiles provide finite
stellar gravity and the initial photospheric parameters. Every imposed closure
is recorded with the resulting cache and budgets.
"""

from pathlib import Path

import argparse,hashlib,json,time,sys

import numpy as np

from scipy.integrate import cumulative_trapezoid

from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

from src.stellar_structure import read_mesa_profile,G,MSUN,RSUN

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def smooth(x):x=np.clip(x,0,1);return x*x*(3-2*x)

class Scenario:
    def __init__(self,config):
        self.config=config;v=ROOT/'data/stellar-structures'
        self.donor=json.loads((v/'donor-fine-structure.json').read_text());self.companion=json.loads((v/'companion-fine-structure.json').read_text())
        self.M=self.donor['mass_msun'];self.m=self.companion['mass_msun'];self.R=self.donor['radius_rsun'];self.r2=self.companion['radius_rsun']
        profile=ROOT/'data/stellar-structures/profiles/donor-fine.profile.data';_,p=read_mesa_profile(profile)
        self.r=np.r_[0.,10**p['logR']];self.mass=np.r_[0.,p['mass']]
        self.grid=np.geomspace(max(self.r[1],1e-6),self.R,6000);mg=self.enclosed(self.grid)
        integ=cumulative_trapezoid(mg/self.grid**2,self.grid,initial=0.)
        self.phi=-self.M/self.R-(integ[-1]-integ)
        self.sources=dict(donor_profile=str(profile.relative_to(ROOT)),donor_profile_sha256=sha(profile),companion_profile_sha256=self.companion['sha256'])
        self.time_unit=np.sqrt(RSUN**3/(G*MSUN));self.energy_unit=G*MSUN**2/RSUN
    def enclosed(self,r):return np.interp(np.asarray(r)**3,self.r**3,self.mass,left=0,right=self.M)
    def potential(self,r):
        r=np.asarray(r);inside=np.interp(r,self.grid,self.phi,left=self.phi[0],right=self.phi[-1]);return np.where(r>self.R,-self.M/np.maximum(r,1e-9),inside)
    def state(self,t):
        c=self.config;u=np.clip(t/c['inspiral_end_code'],0,1)
        a=c['final_separation_rsun']+(c['initial_separation_rsun']-c['final_separation_rsun'])*(1-u**3)**1.5
        disruption=float(smooth((u-.84)/.16));alive=self.m*(1-disruption)
        radius=self.R*(1+.65*smooth((t-c['inspiral_end_code'])/130.))
        return a,disruption,alive,radius
    def omega(self,t):
        a,_,alive,_=self.state(t);return np.sqrt((float(self.enclosed(a))+alive)/a**3)
    def centers(self,t,phase):
        a,_,alive,_=self.state(t);e=np.array([np.cos(phase),np.sin(phase),0.]);xd=-a*alive/(self.M+self.m)*e;return xd,xd+a*e
    def drive(self,t):
        a,_,_,_=self.state(t);a0=self.config['initial_separation_rsun']
        energy=lambda x:self.m*(float(self.potential(x))+.5*float(self.enclosed(x))/x)
        return max(0.,energy(a0)-energy(a))
    def acceleration(self,x,t,phase,gas_mass):
        xd,xc=self.centers(t,phase);_,disruption,alive,_=self.state(t)
        q=x-xd;r=np.maximum(np.linalg.norm(q,axis=1),1e-8);mg=self.enclosed(r)*(1-gas_mass/self.M)+self.m*disruption*(r/np.maximum(r,.3))**3
        a=-mg[:,None]*q/r[:,None]**3
        q=x-xc;r=np.maximum(np.linalg.norm(q,axis=1),1e-8);a-=alive*q/np.maximum(r,self.r2)[:,None]**3
        return a
