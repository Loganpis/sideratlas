"""Inertial, source-resolved reduced transport; NOT a hydrodynamic solver.

Prescribed orbit, pressure-work launch closure and local dissipative encounters
are explicit. No gas path is rotated with the stars after launch.
"""
import numpy as np
from scipy.optimize import brentq
from src.dynamics.prescribed import Scenario,smooth

def sphere_phi(r,m,h):
    r=np.maximum(r,1e-12)
    return np.where(r>=h,-m/r,-m*(3-(r/h)**2)/(2*h))

class TransportScenario(Scenario):
    def __init__(self,config):
        super().__init__(config);self.gas_mass=0.;self.secondary_accreted=0.
    def components(self,t):
        _,dis,_,_=super().state(t)
        secondary=(self.m+self.secondary_accreted)*(1-dis)
        profile=self.M-self.gas_mass-self.secondary_accreted
        core=(self.m+self.secondary_accreted)*dis
        return profile,core,secondary
    def state(self,t):
        a,dis,_,R=super().state(t);return a,dis,self.components(t)[2],R
    def centers(self,t,phase):
        a,_,m2,_=self.state(t);total=self.M+self.m-self.gas_mass
        e=np.array([np.cos(phase),np.sin(phase),0.]);return -a*m2/total*e,a*(1-m2/total)*e
    def omega(self,t):
        a,_,m2,_=self.state(t);mp,mc,_=self.components(t)
        inner=self.enclosed(a)*mp/self.M+mc*min(1.,(a/.3)**3)
        return np.sqrt((inner+m2)/a**3)
    def center_velocities(self,t,phase):
        h=1e-4;om=self.omega(t)
        before=self.centers(t-h,phase-om*h);after=self.centers(t+h,phase+om*h)
        return tuple((b-a)/(2*h) for a,b in zip(before,after))
    def potential_field(self,x,t,phase):
        p,c,m2=self.components(t);x1,x2=self.centers(t,phase)
        r=np.linalg.norm(x-x1,axis=-1);q=np.linalg.norm(x-x2,axis=-1)
        return p/self.M*self.potential(r)+sphere_phi(r,c,.3)+sphere_phi(q,m2,max(self.r2*(m2/self.m)**(1/3),1e-5))
    def force(self,x,t,phase):
        p,c,m2=self.components(t);x1,x2=self.centers(t,phase)
        dr=x-x1;r=np.maximum(np.linalg.norm(dr,axis=-1),1e-10)
        gravity=(p/self.M*self.enclosed(r)/r**3+c/np.maximum(r,.3)**3)
        out=-gravity[:,None]*dr;dr=x-x2;r=np.maximum(np.linalg.norm(dr,axis=-1),1e-10)
        h=max(self.r2*(m2/self.m)**(1/3),1e-5)
        return out-m2*dr/np.maximum(r,h)[:,None]**3
    def lagrange(self,t,phase):
        a,_,m2,_=self.state(t)
        if a<1.05*self.R or m2<.1*self.m:return None
        x1,x2=self.centers(t,0.);om=self.omega(t)
        def f(z):return float(self.force(np.array([[z,0.,0.]]),t,0.)[0,0]+om*om*z)
        return np.array([brentq(f,x1[0]+.01*a,x2[0]-.01*a),brentq(f,x2[0]+.01*a,x2[0]+2*a)])

def surface_radius(direction,R,tensor,omega,mass,engulf):
    flatten=min(.15,.35*omega*omega*R**3/mass)*engulf
    vals,basis=np.linalg.eigh(np.eye(3)+tensor+np.diag([flatten/2,flatten/2,-flatten]));axes=R*vals/np.prod(vals)**(1/3)
    return 1/np.sqrt(np.sum(((direction@basis)/axes)**2,axis=-1))

def collision_impulse(x,v,u,ids,mass,radius,restitution,births,time,event_log=None):
    """Disjoint approaching pairs; central impulses conserve vector J.

    A coarse parcel encounter/shock closure, not resolved fluid shocks. Same
    injection packets and subsonic/co-moving neighbors do not collide.
    """
    from scipy.spatial import cKDTree
    if len(ids)<2:return 0.,0,0.
    pairs=cKDTree(x[ids]).query_pairs(radius,output_type='ndarray')
    if not len(pairs):return 0.,0,0.
    ii=ids[pairs[:,0]];jj=ids[pairs[:,1]];dr=x[jj]-x[ii];r=np.linalg.norm(dr,axis=1);n=dr/np.maximum(r[:,None],1e-12);vn=np.sum((v[jj]-v[ii])*n,axis=1)
    good=(vn<-.035)&(abs(births[ii]-births[jj])>3.)&(r>1e-8)
    ii=ii[good];jj=jj[good];n=n[good];vn=vn[good];r=r[good]
    used=set();loss=0.;events=0;Jerror=0.
    for k in np.argsort(r):
        i=int(ii[k]);j=int(jj[k])
        if i in used or j in used:continue
        impulse=.5*(1+restitution)*vn[k]*n[k];before_v=v[[i,j]].copy();before_u=u[[i,j]].copy()
        before=np.cross(x[i],v[i])+np.cross(x[j],v[j]);v[i]+=impulse;v[j]-=impulse
        delta=.25*mass*(1-restitution**2)*vn[k]**2
        u[i]+=delta/(2*mass);u[j]+=delta/(2*mass);loss+=delta;events+=1;used.update([i,j])
        if event_log is not None:event_log.append(dict(t=float(time),ids=[i,j],position=x[[i,j]].tolist(),velocity_before=before_v.tolist(),velocity_after=v[[i,j]].tolist(),internal_before=before_u.tolist(),internal_after=u[[i,j]].tolist(),kinetic_loss=float(delta)))
        Jerror=max(Jerror,float(np.linalg.norm(np.cross(x[i],v[i])+np.cross(x[j],v[j])-before)*mass))
    return loss,events,Jerror
