"""Event-aware inertial trajectories and a non-connected Eulerian mass grid."""
from pathlib import Path
import json,numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d
from src.gas_transport import TransportScenario
from src.rendering.photosphere_base import cic
from src.stellar_structure import RSUN

class MotionData:
    def enable_force_quintic(self,path):
        """Optional display interpolation; original integrated arrays unchanged."""
        a=np.load(path);self.playback_acc={k:a[k] for k in a.files}
    def __init__(self,source,data=None,report=None):
        self.source=Path(source);self.report=report if report is not None else json.loads((self.source/'reduced-model-report.json').read_text())
        if data is None:
            a=np.load(self.source/'reduced-merger-cache.npz');self.data={k:a[k] for k in a.files}
        else:self.data=data
        self.times=self.data['t'];self.mass=float(self.data['parcel_mass_msun']);self.n=self.data['x'].shape[1]
        self.gas_history=np.array([q['gas_mass'] for q in self.report['diagnostics']])
        self.clock=PchipInterpolator([0,8,16,24,32],[0,140,210,300,390]);self.phase=PchipInterpolator(self.times,self.data['phase']);self.scene=TransportScenario(self.report['config'])
        self.birth=self.data['birth_times'];self.death=self.data['death_times'];self.initial_x=np.zeros((self.n,3));self.initial_v=self.initial_x.copy();self.initial_u=np.zeros(self.n);self.accepted=np.zeros(self.n,bool)
        self.launches=json.loads((self.source/'exact-launch-records.json').read_text())
        for q in self.launches:
            j=q['id'];self.initial_x[j]=q['position'];self.initial_v[j]=q['velocity'];self.initial_u[j]=q['internal_energy'];self.accepted[j]=True
        self.events=json.loads((self.source/'encounter-events.json').read_text());self.event_time=np.array([q['t'] for q in self.events])
        self.actual_J0=self.report['diagnostics'][0]['orbital_J']+self.report['diagnostics'][0]['envelope_J']+self.report['diagnostics'][0]['unresolved_stellar_J_reservoir'][2]
    def state(self,t):
        self.scene.gas_mass=float(np.interp(t,self.times,self.gas_history))
        self.scene.secondary_accreted=float(np.interp(t,self.times,self.data['secondary_accreted']))
        return self.scene
    def at(self,t):
        a=max(0,min(len(self.times)-2,np.searchsorted(self.times,t)-1));b=a+1;ta=self.times[a];tb=self.times[b]
        active=self.accepted&(self.birth<=t)&(~np.isfinite(self.death)|(self.death>t))
        lo=np.full(self.n,ta);hi=np.full(self.n,tb);x0=self.data['x'][a].astype('f8').copy();x1=self.data['x'][b].astype('f8').copy();v0=self.data['v'][a].astype('f8').copy();v1=self.data['v'][b].astype('f8').copy();u0=self.data['specific_internal_energy'][a].astype('f8').copy();u1=self.data['specific_internal_energy'][b].astype('f8').copy()
        quintic=hasattr(self,'playback_acc')
        if quintic:a0=self.playback_acc['snapshots'][a].copy();a1=self.playback_acc['snapshots'][b].copy()
        born=(self.birth>ta)&self.accepted;lo[born]=self.birth[born];x0[born]=self.initial_x[born];v0[born]=self.initial_v[born];u0[born]=self.initial_u[born]
        died=np.isfinite(self.death)&(self.death<tb);hi[died]=self.death[died]
        if quintic:a0[born]=self.playback_acc['birth'][born];a1[died]=self.playback_acc['death'][died]
        # Collision positions are continuous; endpoint velocities have genuine
        # impulses. Never interpolate through a discontinuity with one cubic.
        # Include impulses exactly on saved snapshot boundaries. The left
        # interval must end at the PRE-impulse velocity; its cached endpoint
        # otherwise contains the POST-impulse state and smears the jump over
        # the entire preceding interval. The right side starts post-impulse.
        for k in np.flatnonzero((self.event_time>=ta)&(self.event_time<=tb)):
            q=self.events[k];et=q['t']
            for m,j in enumerate(q['ids']):
                if et<=t and et>=lo[j]:lo[j]=et;x0[j]=q['position'][m];v0[j]=q['velocity_after'][m];u0[j]=q['internal_after'][m]
                elif et>t and et<=hi[j]:hi[j]=et;x1[j]=q['position'][m];v1[j]=q['velocity_before'][m];u1[j]=q['internal_before'][m]
                else:continue
                if quintic:
                    if et<=t:a0[j]=self.playback_acc['events'][k,m]
                    else:a1[j]=self.playback_acc['events'][k,m]
        h=np.maximum(hi-lo,1e-12);f=np.clip((t-lo)/h,0,1);f2=f*f;f3=f2*f
        x=(2*f3-3*f2+1)[:,None]*x0+((f3-2*f2+f)*h)[:,None]*v0+(-2*f3+3*f2)[:,None]*x1+((f3-f2)*h)[:,None]*v1
        v=((6*f2-6*f)/h)[:,None]*x0+(3*f2-4*f+1)[:,None]*v0+((-6*f2+6*f)/h)[:,None]*x1+(3*f2-2*f)[:,None]*v1
        if quintic:
            H=h[:,None];F=f[:,None];D=x1-x0-v0*H-.5*a0*H*H;V=(v1-v0)*H-a0*H*H;A=(a1-a0)*H*H;c3=10*D-4*V+.5*A;c4=-15*D+7*V-A;c5=6*D-3*V+.5*A
            x=x0+v0*H*F+.5*a0*H*H*F**2+c3*F**3+c4*F**4+c5*F**5;v=(v0*H+a0*H*H*F+3*c3*F**2+4*c4*F**3+5*c5*F**4)/H
        u=np.exp((1-f)*np.log(np.maximum(u0,1e-20))+f*np.log(np.maximum(u1,1e-20)))
        s=self.state(t);phase=float(self.phase(t));xd,xc=s.centers(t,phase);vd=np.array([np.interp(t,self.times,self.data['donor_velocity'][:,j]) for j in range(3)])
        rel=x-xd;r=np.linalg.norm(rel,axis=1);vr=np.sum((v-vd)*rel,axis=1)/np.maximum(r,1e-12);energy=.5*np.sum(v*v,axis=1)+s.potential_field(x,t,phase)
        fc=np.clip((t-ta)/(tb-ta),0,1);rmax=np.maximum(r,(1-fc)*self.data['maximum_radius'][a]+fc*self.data['maximum_radius'][b])
        fallback=active&(vr<-.002)&(rmax>1.1*self.data['initial_radii'])
        classes=np.full(self.n,3,np.int8);classes[active&(energy<=0)]=1;classes[active&(energy>0)&(vr>0)]=0;classes[fallback]=2
        unitv=(RSUN/100)/self.report['physical_time_unit_seconds'];T=np.maximum(0.,u*(2/3)*1.3*1.67262192369e-27*unitv**2/1.380649e-23)
        return dict(t=t,x=x,v=v,u=u,T=T,live=active,radial_velocity=vr,energy=energy,classes=classes,phase=phase,donor=xd,companion=xc,donor_velocity=vd,radius=s.state(t)[3],m2=s.state(t)[2])
    def camera_curve(self):
        ts=np.linspace(0,390,161);width=[]
        for t in ts:
            q=self.at(float(t));x=q['x'][q['live']]
            inc=.82-.12*np.clip((t-215)/175,0,1);az=-np.pi/2+.35*(np.interp(t,[0,140,210,300,390],[0,.25,.5,.75,1]))
            eye=np.array([np.sin(inc)*np.cos(az),np.sin(inc)*np.sin(az),np.cos(inc)]);right=np.array([-np.sin(az),np.cos(az),0]);up=np.cross(eye,right)
            needed=np.maximum(abs(x@right),np.where(x@up>=0,(x@up)/(1.28*16/9),-(x@up)/(.72*16/9))) if len(x) else np.zeros(1)
            width.append(max(6.85,float(np.quantile(needed,.95))*1.15))
        width=np.maximum.accumulate(width);width=gaussian_filter1d(width,2,mode='nearest')
        return PchipInterpolator(ts,width)

def mass_grid(points,masses,temperatures,halfwidth,gw=96,nz=80,center_y=.64,eye=None,right=None,up=None):
    """CIC + fixed one-cell anti-alias basis, no material connectors or ages.

    Every contribution is local to its actual inertial position. The only
    broadening is the declared Eulerian resolution response, not a gas packet.
    """
    from scipy.ndimage import gaussian_filter
    if eye is None:eye=np.array([0.,0.,1.]);right=np.array([1.,0.,0.]);up=np.array([0.,1.,0.])
    camera=np.column_stack([right,up,eye]);p=points@camera;gh=round(gw*16/9);pad=4;dx=2*halfwidth/gw
    Z=max(halfwidth,float(np.quantile(abs(p[:,2]),.995))+dx) if len(p) else halfwidth
    dz=2*Z/nz;origin=np.array([-halfwidth-pad*dx,-(1-center_y)*gh*dx-pad*dx,-Z-pad*dz]);spacing=np.array([dx,dx,dz]);shape=(nz+2*pad,gh+2*pad,gw+2*pad)
    m,h=cic(p,masses,temperatures,shape,origin,spacing)
    # A constant sub-voxel reconstruction response. It never follows velocity,
    # birth ID, stellar rotation or parcel age.
    m=gaussian_filter(m,.65,mode='constant',truncate=4);h=gaussian_filter(h,.65,mode='constant',truncate=4)
    zc=origin[2]+(np.arange(shape[0])+.5)*dz;zs=np.full(shape[0],dz)
    return m,h,origin,spacing,zc,zs
