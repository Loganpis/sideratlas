"""Pressure-work surface sources + inertial gravity + dissipative encounters.

This intentionally remains a prescribed-orbit reduced model. Reports include
external-potential work and mass-source/sink energies, not a false closed binary.
"""
from pathlib import Path
import sys,json,time,argparse,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import ndtri
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.gas_transport import TransportScenario,surface_radius,collision_impulse
from src.dynamics.prescribed import smooth
from src.stellar_structure import G,MSUN,RSUN

def build(config,out):
    begin=time.monotonic();s=TransportScenario(config);n=config['gas_parcels'];mass=config['gas_mass_msun']/n
    step=config['integration_step_code'];end=config['duration_code'];gamma=5/3;unitv=(RSUN/100)/s.time_unit
    photo_c2=1.380649e-23*s.donor['effective_temperature_k']/(1.3*1.67262192369e-27)/unitv**2;photo_u=photo_c2/(gamma-1)
    q=s.M/s.m;f=.49*q**(2/3)/(.6*q**(2/3)+np.log1p(q**(1/3)))
    guide=np.linspace(0,end,10001);separation=np.array([s.state(t)[0] for t in guide]);start=float(np.interp(s.R/f,separation[::-1],guide[::-1]))
    R=s.R;bury=smooth((1.10*R-separation)/(.65*R));rate=(.12+bury**2)*smooth((guide-start)/15)*(1-smooth((guide-295)/20))
    cumulative=np.r_[0,np.cumsum((rate[1:]+rate[:-1])/2*np.diff(guide))];cumulative/=cumulative[-1]
    # Stratified chronological packets. Each packet samples a finite photosphere
    # and velocity distribution; none is an arm attached to a rotating object.
    packet=8;births=np.interp(((np.arange(n)//packet)+.5)/(n/packet),cumulative,guide)
    z=ndtri(np.clip(qmc.Sobol(6,scramble=True,seed=5005).random_base2(int(np.ceil(np.log2(n))))[:n],1e-5,1-1e-5))
    drive=np.array([s.drive(t) for t in guide]);drive_rate=np.maximum(0,np.gradient(drive,guide));mdot=rate/config.get('mass_rate_normalization',1.)
    mdot*=config['gas_mass_msun']/np.trapz(rate,guide)
    targets=np.unique(np.r_[np.linspace(0,end,config['cache_frames']),PchipInterpolator([0,8,16,24,32],[0,140,210,300,390])((np.arange(384)+.5)/12),[300,340,390]])
    x=np.zeros((n,3));v=np.zeros_like(x);u=np.zeros(n);live=np.zeros(n,bool);launched=np.zeros(n,bool);rmax=np.zeros(n);r0=np.ones(n);death=np.full(n,np.nan);sink_kind=np.zeros(n,np.int8)
    om0=s.omega(0);I0=s.donor['axial_inertia_g_cm2']/(MSUN*RSUN**2);J0=s.M*s.m/(s.M+s.m)*config['initial_separation_rsun']**2*om0+I0*om0
    phase=0.;pattern=0.;envelope=om0;t=0.;tensor=np.zeros((3,3));heat=0.;escaped=0.;launch_work=0.;external_work=0.;mechanical_sources=0.;mechanical_sinks=0.;internal_sources=0.;internal_sinks=0.;gas_radiated=0.;expansion_work=0.;shock_loss=0.;collision_count=0;collision_J_error=0.;sinkJ=np.zeros(3)
    snapshots=[];logs=[];launch_log=[];encounter_log=[];ci=0;next_collision=.25;tracked=np.linspace(0,n-1,48,dtype=int);track=[];next_track=0.;lastdrive=0.;snapshot_index=0
    def state_energy(ids,tt,ph):return mass*(.5*np.sum(v[ids]**2,axis=1)+s.potential_field(x[ids],tt,ph))
    def record():
        nonlocal snapshots
        a,dis,m2,R=s.state(t);xd,xc=s.centers(t,phase);vd,vc=s.center_velocities(t,phase);ids=np.flatnonzero(live);r=np.linalg.norm(x-xd,axis=1);vr=np.sum((v-vd)*(x-xd),axis=1)/np.maximum(r,1e-12)
        temperature=np.maximum(300.,u*(gamma-1)*1.3*1.67262192369e-27*unitv**2/1.380649e-23)
        energy=np.zeros(n);energy[ids]=.5*np.sum(v[ids]**2,axis=1)+s.potential_field(x[ids],t,phase)
        # A reconstruction footprint, NOT age-broadened material geometry.
        width=np.full(n,config['reconstruction_sigma_Rsun'])
        E=float(state_energy(ids,t,phase).sum());Eint=float(mass*u[ids].sum())
        Jgas=mass*np.sum(np.cross(x[ids],v[ids]),axis=0);M1=s.M+s.m-s.gas_mass-m2;Jorb=M1*m2/(M1+m2)*a*a*s.omega(t);I=I0*max(.01,s.components(t)[0]/s.M)*(R/s.R)**2
        residualJ=np.array([0.,0.,J0-Jorb-I*envelope])-Jgas
        mass_error=M1+m2+s.gas_mass-s.M-s.m
        mechanical_error=E-(mechanical_sources+external_work-mechanical_sinks-shock_loss)
        total_error=E+Eint+gas_radiated+expansion_work-(mechanical_sources+internal_sources+external_work-mechanical_sinks-internal_sinks)
        snap=dict(t=t,phase=phase,launch_pattern_phase=pattern,envelope_omega=envelope,separation=a,donor_center=xd,companion_center=xc,donor_velocity=vd,companion_velocity=vc,donor_radius=R,companion_mass=m2,primary_mass=M1,secondary_accreted=s.secondary_accreted,
          x=x.copy().astype('f4'),v=v.copy().astype('f4'),live=live.copy(),temperature=temperature.astype('f4'),specific_internal_energy=u.astype('f4'),width=width.astype('f4'),radial_velocity=vr.astype('f4'),specific_energy=energy.astype('f4'),maximum_radius=rmax.astype('f4'),
          gas_luminosity_code=.5*heat/config['thermal_diffusion_time_code']+mass*u[ids].sum()/config['gas_cooling_time_code'],remnant_luminosity_code=.5*heat/config['thermal_diffusion_time_code'],heat_code=heat,available_orbital_work_code=s.drive(t),disruption=dis,tidal_tensor=tensor.copy())
        snapshots.append(snap);logs.append(dict(t=t,gas_mass=s.gas_mass,primary_mass=M1,secondary_mass=m2,mass_error=mass_error,mechanical_energy=E,internal_energy=Eint,external_potential_work=external_work,mechanical_sources=mechanical_sources,mechanical_sinks=mechanical_sinks,internal_sources=internal_sources,internal_sinks=internal_sinks,shock_dissipation=shock_loss,gas_radiated=gas_radiated,expansion_work=expansion_work,mechanical_residual=mechanical_error,total_gas_residual=total_error,gas_J=Jgas.tolist(),orbital_J=Jorb,envelope_J=I*envelope,unresolved_stellar_J_reservoir=residualJ.tolist(),source_pressure_work=launch_work,source_work_allowance=config['ejecta_launch_energy_fraction_ceiling']*s.drive(t),collision_count=collision_count,heat=heat,escaped_heat=escaped))
    while t<end+1e-9:
        if t>=next_track-1e-8:
            track.append((t,x[tracked].copy().astype('f4'),v[tracked].copy().astype('f4'),live[tracked].copy()));next_track+=.1
        oldids=np.flatnonzero(live);phi_before=s.potential_field(x[oldids],t,phase)
        newborn=np.flatnonzero((births<=t+1e-8)&~launched)
        for j in newborn:
            launched[j]=True;a,dis,m2,R=s.state(t);xd,xc=s.centers(t,phase);vd,_=s.center_velocities(t,phase);buried=float(smooth((1.10*R-a)/(.65*R)));om=s.omega(t)
            u_shock=min(.8*(s.M+s.m-s.gas_mass)/R,.10*float(np.interp(t,guide,drive_rate))/max(float(np.interp(t,guide,mdot)),1e-9))
            fan=min(1.15,np.sqrt(gamma*(gamma-1)*u_shock)*config['envelope_response_time_code']/R)
            az=(phase if buried<.05 else pattern)-buried*np.arctan(om*config['envelope_response_time_code'])
            az+=z[j,0]*(.025+buried*fan);lat=np.clip(z[j,1]*(.025+buried*.65*fan),-1.1,1.1)
            direction=np.array([np.cos(az)*np.cos(lat),np.sin(az)*np.cos(lat),np.sin(lat)]);tangent=np.array([-np.sin(az),np.cos(az),0.]);vertical=np.cross(direction,tangent)
            engulf=float(smooth((1.05*R-a)/(.65*R)));rr=float(surface_radius(direction[None],R,tensor,envelope,s.M+dis*s.m,engulf)[0]);pos=xd+rr*1.012*direction
            base=vd+envelope*np.cross(np.array([0.,0.,1.]),pos-xd)
            family=0 if buried<.05 else 1
            if family==0:
                lp=s.lagrange(t,phase);radial=np.sqrt(photo_c2)
                if lp is not None:
                    saddle=lp[0]*np.array([np.cos(phase),np.sin(phase),0.]);psi1=s.potential_field(saddle[None],t,phase)[0]-.5*om*om*np.dot(saddle,saddle);psis=s.potential_field(pos[None],t,phase)[0]-.5*om*om*np.dot(pos,pos)
                    radial=np.sqrt(photo_c2+2*max(0.,psi1-psis))
                kick=radial*direction+np.sqrt(photo_c2)*(.45*z[j,2]*tangent+.45*z[j,3]*vertical);extra_u=0.
            else:
                sigma=np.sqrt(.15*u_shock+photo_c2)
                radial=np.sqrt(max(0.,u_shock))+.55*sigma*z[j,2]
                kick=max(np.sqrt(photo_c2),radial)*direction+sigma*(z[j,3]*tangent+z[j,4]*vertical);extra_u=.25*u_shock
            allowance=max(0.,config['ejecta_launch_energy_fraction_ceiling']*s.drive(t)-launch_work)
            def cost(alpha):return mass*(max(0.,.5*(np.dot(base+alpha*kick,base+alpha*kick)-np.dot(base,base)))+alpha*alpha*extra_u)
            alpha=1.
            if cost(1)>allowance:
                lo=0.;hi=1.
                for _ in range(30):
                    mid=(lo+hi)/2
                    if cost(mid)>allowance:hi=mid
                    else:lo=mid
                alpha=lo
            if alpha<.1:death[j]=t;sink_kind[j]=6;continue
            live[j]=True;x[j]=pos;v[j]=base+alpha*kick;u[j]=photo_u+alpha*alpha*extra_u;r0[j]=np.linalg.norm(pos-xd);rmax[j]=r0[j];s.gas_mass+=mass;launch_work+=cost(alpha)
            launch_log.append(dict(id=int(j),t=t,family='surface_Roche_pressure_nozzle' if family==0 else 'buried_wake_surface_fan',position=pos.tolist(),velocity=v[j].tolist(),internal_energy=float(u[j]),primary_velocity=vd.tolist(),inherited_velocity=base.tolist(),pressure_velocity=(alpha*kick).tolist(),pressure_scale=alpha,pressure_specific_energy=u_shock,launch_radius=r0[j],surface_radius=rr,latitude=lat,fan_sigma=fan*buried,launch_J=np.cross(pos,v[j]).tolist()))
        newids=np.flatnonzero(live)
        if len(newborn):
            if len(oldids):external_work+=float(mass*np.sum(s.potential_field(x[oldids],t,phase)-phi_before))
            accepted=newborn[live[newborn]];mechanical_sources+=float(state_energy(accepted,t,phase).sum());internal_sources+=float(mass*u[accepted].sum())
        if snapshot_index<len(targets) and t>=targets[snapshot_index]-1e-8:
            record();snapshot_index+=1
        if t>=end-1e-9:break
        # Land on exact cache and birth times; no gas interpolation from an
        # unborn placeholder and no missing sub-frame source motion.
        pending_birth=births[~launched];next_birth=float(pending_birth[0]) if len(pending_birth) else end
        dt=min(step,end-t,float(targets[snapshot_index])-t if snapshot_index<len(targets) else step,max(1e-9,next_birth-t) if next_birth>t+1e-9 else step)
        ids=np.flatnonzero(live);om=s.omega(t);phase_next=phase+om*dt;om_next=s.omega(t+dt);phase_next=phase+.5*(om+om_next)*dt
        Jgas=mass*np.sum(np.cross(x[ids],v[ids]),axis=0) if len(ids) else np.zeros(3);a,dis,m2,R=s.state(t);M1=s.M+s.m-s.gas_mass-m2;Jorb=M1*m2/(M1+m2)*a*a*om;I=I0*max(.01,s.components(t)[0]/s.M)*(R/s.R)**2
        spin_available=J0-Jorb-Jgas[2];target=min((I0*om0+config['orbital_J_to_envelope_fraction']*max(0.,J0-I0*om0-Jorb))/I,max(0.,spin_available)/I,config['launch_breakup_fraction']*np.sqrt((M1)/R**3))
        new_envelope=min(envelope+(1-np.exp(-dt/config['envelope_response_time_code']))*(target-envelope),max(0.,spin_available)/I)
        b=float(smooth((1.10*R-a)/(.65*R)));source_rate=(1-b)*om+b*.5*(envelope+new_envelope);pattern+=source_rate*dt;envelope=new_envelope
        e=np.array([np.cos(phase),np.sin(phase),0.]);eps=min(.16,.75*(m2/max(M1,1e-8))*(R/max(a,R))**3)*(1-dis);tidal_target=eps*(3*np.outer(e,e)-np.eye(3));tensor+=(1-np.exp(-dt/config['envelope_response_time_code']))*(tidal_target-tensor)
        d=s.drive(t+dt);heat+=config['thermal_deposition_fraction']*(d-lastdrive);lastdrive=d;cool=heat*(-np.expm1(-dt/config['thermal_diffusion_time_code']));heat-=cool;escaped+=cool
        if len(ids):
            before=x[ids].copy();v[ids]+=.5*dt*s.force(before,t,phase);x[ids]+=dt*v[ids];v[ids]+=.5*dt*s.force(x[ids],t+dt,phase_next)
            midpoint=.5*(before+x[ids]);external_work+=float(mass*np.sum(s.potential_field(midpoint,t+dt,phase_next)-s.potential_field(midpoint,t,phase)))
            xd,xc=s.centers(t+dt,phase_next);vd,_=s.center_velocities(t+dt,phase_next);r=np.linalg.norm(x[ids]-xd,axis=1);newmax=np.maximum(rmax[ids],r);lossu=u[ids]*(1-(np.maximum(newmax,r0[ids])/np.maximum(rmax[ids],r0[ids]))**(-4/3));u[ids]-=lossu;expansion_work+=float(mass*lossu.sum());rmax[ids]=newmax
            lossu=u[ids]*(-np.expm1(-dt/config['gas_cooling_time_code']));u[ids]-=lossu;gas_radiated+=float(mass*lossu.sum())
            if t+dt>=next_collision-1e-8 and config['collision_restitution']<1:
                lost,count,error=collision_impulse(x,v,u,ids,mass,config['collision_radius_Rsun'],config['collision_restitution'],births,t+dt,encounter_log);shock_loss+=lost;collision_count+=count;collision_J_error=max(collision_J_error,error)
                interval=config.get('collision_interval_code',.1);next_collision=(np.floor((t+dt)/interval)+1)*interval
            a,dis,m2,R=s.state(t+dt);direction=(x[ids]-xd)/np.maximum(r[:,None],1e-12);engulf=float(smooth((1.05*R-a)/(.65*R)));surface=surface_radius(direction,R,tensor,envelope,s.M+dis*s.m,engulf)
            primary_hit=r<surface*.997;secondary_hit=(np.linalg.norm(x[ids]-xc,axis=1)<s.r2*(max(m2,0)/s.m)**(1/3))&~primary_hit
            removed=ids[primary_hit|secondary_hit]
            if len(removed):
                survivor=ids[~(primary_hit|secondary_hit)];phi0=s.potential_field(x[survivor],t+dt,phase_next)
                mechanical_sinks+=float(state_energy(removed,t+dt,phase_next).sum());internal_sinks+=float(mass*u[removed].sum());sinkJ+=mass*np.sum(np.cross(x[removed],v[removed]),axis=0)
                vr=np.sum((v[ids]-vd)*direction,axis=1)
                for local,j in enumerate(ids):
                    if primary_hit[local]:sink_kind[j]=2 if vr[local]<0 and rmax[j]>1.1*r0[j] else (3 if vr[local]<0 else 4)
                    elif secondary_hit[local]:sink_kind[j]=5
                s.secondary_accreted+=float(secondary_hit.sum()*mass/max(1-dis,1e-8));s.gas_mass-=len(removed)*mass;live[removed]=False;death[removed]=t+dt
                if len(survivor):external_work+=float(mass*np.sum(s.potential_field(x[survivor],t+dt,phase_next)-phi0))
        t+=dt;phase=phase_next
    out.mkdir(parents=True,exist_ok=True);stack=lambda key:np.array([f[key] for f in snapshots]);keys=list(snapshots[0])
    np.savez_compressed(out/'reduced-merger-cache.npz',**{k:stack(k) for k in keys},parcel_mass_msun=mass,birth_times=births,death_times=death,sink_kind=sink_kind,initial_radii=r0)
    np.savez_compressed(out/'tracked-parcels.npz',time=np.array([q[0] for q in track]),ids=tracked,x=np.array([q[1] for q in track]),v=np.array([q[2] for q in track]),live=np.array([q[3] for q in track]))
    denom=max(abs(mechanical_sources)+abs(external_work)+abs(mechanical_sinks)+shock_loss,1e-10);max_error=max(abs(q['mechanical_residual']) for q in logs)
    result=dict(scope=__doc__,config=config,stellar_inputs=s.sources,donor=dict(mass=s.M,radius=s.R,Teff=s.donor['effective_temperature_k']),companion=dict(mass=s.m,radius=s.r2,Teff=s.companion['effective_temperature_k']),physical_time_unit_seconds=s.time_unit,physical_duration_days=end*s.time_unit/86400,elapsed_seconds=time.monotonic()-begin,diagnostics=logs,launch_counts={'total':len(launch_log),'Roche_nozzle':sum(q['family']=='surface_Roche_pressure_nozzle' for q in launch_log),'wake_fan':sum(q['family']=='buried_wake_surface_fan' for q in launch_log)},sink_counts={str(i):int((sink_kind==i).sum()) for i in range(7)},collision_count=collision_count,collision_J_absolute_error=collision_J_error,relative_mechanical_balance_error=max_error/denom,pressure_work=launch_work,launch_allowance=config['ejecta_launch_energy_fraction_ceiling']*s.drive(end),max_mass_error=max(abs(q['mass_error']) for q in logs),start_Roche_contact_code=start,gas_remaining_msun=s.gas_mass,model_is_hydrodynamic=False,classification='Instantaneous inertial energy and primary-relative radial motion; fallback requires previous outward excursion. Not guaranteed asymptotic escape in a time-dependent potential.',limitations=['Prescribed separation/disruption and pressure-work source; no calculated instability or merger','Collision closure is disjoint central restitution impulses, not resolved shocks or gas pressure','External potential work and stellar momentum/angular-momentum reservoirs are reported, not self-consistent orbital backreaction','Kernel size is a reconstruction footprint; it is not an expanding physical packet or a material stream connector'])
    convert=lambda v:v.item() if isinstance(v,np.generic) else str(v)
    (out/'reduced-model-report.json').write_text(json.dumps(result,indent=2,default=convert)+'\n');(out/'exact-launch-records.json').write_text(json.dumps(launch_log,indent=2,default=convert)+'\n')
    (out/'encounter-events.json').write_text(json.dumps(encounter_log,indent=2,default=convert)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['diagnostics','config','stellar_inputs']},indent=2,default=convert),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if (a.output/'reduced-merger-cache.npz').exists():raise RuntimeError('Do not overwrite completed transport result')
    build(json.loads(a.config.read_text()),a.output)
if __name__=='__main__':main()
