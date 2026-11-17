"""Selected calculation and recoverable checkpoints; explicit external gravity.

The original physical loop is retained. Only filesystem layout and the local
preparation/review-file prerequisites are replaced by the selected configuration.
Running this command is expensive and was not repeated for publication.
"""
import argparse,fcntl,hashlib,json,os,time
from pathlib import Path
import numpy as np
from . import ROOT,CONFIG
from .dynamics import BASE,evaluate,orbit_step,initial_cohorts
from .gas import GasDisk
from .starformation import sample_sink
from .spiral import field2,Gravity2,WAVES2
GAS_BUDGETS=['mass_boundary','floor_mass','formed_mass','formed_J','J_boundary','J_external','explicit_work','initial_mass','initial_J']

def atomic_json(path,obj):
    path=Path(path); tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n'); tmp.replace(path)

def run(args):
    out=args.output.resolve();parent=out.parent
    parent.mkdir(parents=True,exist_ok=True)
    lock=(parent/'science.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('Another science calculation is active')
    lock.seek(0);lock.truncate();lock.write(json.dumps(dict(pid=os.getpid(),model='trailing22')));lock.flush()
    if not np.array_equal(WAVES2,np.array(CONFIG['model']['spiral_parameters'])):
        raise RuntimeError('Selected spiral parameters changed')
    if out.exists() and not args.resume:raise RuntimeError('Calculation exists; preserve it, do not duplicate')
    if args.resume and (out/'run-report.json').exists():raise RuntimeError('Calculation already finished; reuse its results')
    out.mkdir(exist_ok=True);(out/'snapshots').mkdir(exist_ok=True)
    waves=WAVES2.copy();gas=GasDisk(args.grid,gravity_class=lambda X,Y:Gravity2(X,Y,waves),cs_kms=4.,sf_model='kmt',tdep=1500.,clumping=2.)
    capacity=600000;yx=np.zeros((capacity,3));yv=yx.copy();ya=np.zeros((capacity,6));ym=np.zeros(capacity);yb=np.zeros(capacity)
    code=['src/dynamics.py','src/gas.py','src/starformation.py','src/spiral.py','src/simulate.py','config.json']
    hashes={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in code}
    if args.resume:
        with np.load(out/'checkpoint.npz') as cp:
            meta=json.loads(str(cp['metadata']));x=cp['x'].copy();v=cp['v'].copy();a=cp['a'].copy();mass=cp['mass'].copy();cohort=cp['cohort'].copy()
            ny=len(cp['young_mass']);yx[:ny]=cp['young_x'];yv[:ny]=cp['young_v'];ya[:ny]=cp['young_a'];ym[:ny]=cp['young_mass'];yb[:ny]=cp['young_birth']
            gas.q=cp['gas_q'].copy();gas.mass_pending=cp['gas_pending'].copy();gas.momentum_pending=cp['gas_momentum_pending'].copy();gas.KE_pending=float(cp['gas_KE_pending'])
        if meta['source_sha256']!=hashes:raise RuntimeError('Checkpoint source differs; audit before any continuation')
        for k,value in meta['gas'].items():setattr(gas,k,value)
        t=meta['t'];W=meta['W'];J=meta['J'];wy=meta['wy'];jy=meta['jy'];Ey_birth=meta['Ey_birth'];Jy_birth=meta['Jy_birth'];E0=meta['E0'];J0=meta['J0']
        spawn=meta['spawn'];snap_at=meta['snap_at'];gas_at=meta['gas_at'];diag_at=meta['diag_at'];records=meta['records'];gas_records=meta['gas_records'];budgets=meta['budgets'];birth_controls=meta['birth_controls']
        rng=np.random.default_rng();rng.bit_generator.state=meta['rng'];elapsed=meta['wall_seconds'];segment=meta['segment']+1
    else:
        if 200000%args.stars:raise ValueError('Use an exact divisor of the accepted200k initial sample')
        x0,v0,c0,m0=initial_cohorts(200000);stride=200000//args.stars
        ids=np.r_[np.arange(0,120000,stride),np.arange(120000,200000,stride)]
        x=x0[ids].copy();v=v0[ids].copy();cohort=c0[ids];mass=m0[ids]*stride;a=evaluate(x,0.,BASE,waves,field2)
        np.savez(out/'cohorts.npz',id=cohort,mass=mass,original_ids=ids)
        E0=float(np.sum(mass*(.5*np.sum(v*v,axis=1)+a[:,3])));J0=float(np.sum(mass*np.cross(x,v)[:,2]))
        W=J=wy=jy=Ey_birth=Jy_birth=0.;ny=0;t=0.;spawn=.5;snap_at=560.;gas_at=560.;diag_at=0.
        records=[];gas_records=[];budgets=[];birth_controls=[];rng=np.random.default_rng(27007);elapsed=0.;segment=0
        config=dict(model='trailing22',stars=len(x),gas_grid=args.grid,gas_box_kpc=36.,
            record_start_Myr=560.,record_step_Myr=2.,gas_record_step_Myr=.5,run_end_Myr=720.,playback_Myr_per_second=5.,
            base_parameters=BASE.tolist(),spiral_parameters=waves.tolist(),source_sha256=hashes,
            initial_sampling='Deterministic every-Nth subset of exact accepted200k realization; masses scaled conservatively; same cohort positions/velocities, not shifted onto spirals',
            gas_closure=dict(effective_cs_kms=4.,sf='KMT two-phase equilibrium approximation',clumping_factor=2.,molecular_depletion_Myr=1500.),
            scope='Calculated orbit/gas response to prescribed scalar potential; no live self-gravity, MHD, feedback, or claim of spontaneously generated spiral structure')
        atomic_json(out/'config.json',config)
    start=time.monotonic();last_cp=start
    while t<720-1e-8:
        dt=min(.125,gas.dt_limit(),spawn-t,720-t)
        if snap_at>t+1e-8:dt=min(dt,snap_at-t)
        gas.step(t,dt);w,j=orbit_step(x,v,a,dt,t,BASE,waves,field2,mass);W+=w;J+=j
        if ny:
            w,j=orbit_step(yx[:ny],yv[:ny],ya[:ny],dt,t,BASE,waves,field2,ym[:ny]);wy+=w;jy+=j
        t+=dt
        if t>=spawn-1e-8:
            pos,vel,weights,control=sample_sink(gas,rng);count=len(weights)
            if ny+count>capacity:raise RuntimeError('Association capacity reached; preserve checkpoint before a targeted capacity review')
            yx[ny:ny+count]=pos;yv[ny:ny+count]=vel;ym[ny:ny+count]=weights;yb[ny:ny+count]=t
            ya[ny:ny+count]=evaluate(pos,t,BASE,waves,field2)
            Ey_birth+=float(np.sum(weights*(.5*np.sum(vel*vel,axis=1)+ya[ny:ny+count,3])))
            Jy_birth+=float(np.sum(weights*np.cross(pos,vel)[:,2]));ny+=count;birth_controls.append(control)
            gas.mass_pending[:]=0;gas.momentum_pending[:]=0;gas.KE_pending=0.;spawn+=.5
        if t>=snap_at-1e-8:
            name=f'snap-{len(records):03d}-s{segment}.npz'
            np.savez(out/'snapshots'/name,t=t,old_x=x.astype('f4'),old_v=v.astype('f4'),young_x=yx[:ny].astype('f4'),young_v=yv[:ny].astype('f4'),gas_v=(gas.q[1:]/gas.q[0]).astype('f4'))
            records.append(dict(t_Myr=float(t),file=name,young_count=ny));snap_at+=2.
        if t>=gas_at-1e-8:
            name=f'gas-{len(gas_records):03d}-s{segment}.npy';np.save(out/'snapshots'/name,gas.q[0].astype('f4'))
            gas_records.append(dict(t_Myr=float(t),file=name));gas_at+=.5
        if t>=diag_at-1e-8 or t>=720-1e-8:
            E=float(np.sum(mass*(.5*np.sum(v*v,axis=1)+a[:,3])));L=float(np.sum(mass*np.cross(x,v)[:,2]))
            ey=float(np.sum(ym[:ny]*(.5*np.sum(yv[:ny]*yv[:ny],axis=1)+ya[:ny,3])));ly=float(np.sum(ym[:ny]*np.cross(yx[:ny],yv[:ny])[:,2]))
            report=gas.report(t);report.update(stellar_energy_relative_residual=(E-E0-W)/abs(E0),stellar_J_relative_residual=(L-J0-J)/abs(J0),
                young_energy_relative_residual=(ey-Ey_birth-wy)/max(abs(Ey_birth),1.),young_J_relative_residual=(ly-Jy_birth-jy)/max(abs(Jy_birth),1.),
                old_stellar_J=L,young_stellar_J=ly,old_external_torque=J,young_external_torque=jy,
                young_associations=ny,wall_seconds=elapsed+time.monotonic()-start,
                cohort_z_rms_kpc=[float(np.sqrt(np.mean(x[cohort==k,2]**2))) for k in (0,1)])
            budgets.append(report);atomic_json(out/'progress.json',report);diag_at+=10.
            print(f'{args.model} t={t:7.3f}/720Myr wall={report["wall_seconds"]:.1f}s gasJ={report["J_relative_residual"]:.3e} sources={ny}',flush=True)
        if time.monotonic()-last_cp>60 or t>=720-1e-8:
            meta=dict(t=t,W=W,J=J,wy=wy,jy=jy,Ey_birth=Ey_birth,Jy_birth=Jy_birth,E0=E0,J0=J0,spawn=spawn,snap_at=snap_at,gas_at=gas_at,diag_at=diag_at,
                rng=rng.bit_generator.state,gas={k:getattr(gas,k) for k in GAS_BUDGETS},records=records,gas_records=gas_records,budgets=budgets,birth_controls=birth_controls,
                wall_seconds=elapsed+time.monotonic()-start,segment=segment,source_sha256=hashes)
            np.savez(out/'checkpoint.tmp.npz',metadata=json.dumps(meta),x=x,v=v,a=a,cohort=cohort,mass=mass,
                young_x=yx[:ny],young_v=yv[:ny],young_a=ya[:ny],young_mass=ym[:ny],young_birth=yb[:ny],
                gas_q=gas.q,gas_pending=gas.mass_pending,gas_momentum_pending=gas.momentum_pending,gas_KE_pending=gas.KE_pending)
            (out/'checkpoint.tmp.npz').replace(out/'checkpoint.npz');last_cp=time.monotonic()
    np.savez(out/'associations.npz',birth=yb[:ny],mass=ym[:ny])
    atomic_json(out/'snapshots.json',records);atomic_json(out/'gas-snapshots.json',gas_records);atomic_json(out/'diagnostics.json',budgets);atomic_json(out/'birth-controls.json',birth_controls)
    final=budgets[-1];checks=dict(gas_mass=abs(final['mass_relative_closure'])<1e-8,gas_J=abs(final['J_relative_residual'])<.02,
        stellar_work=max(abs(r['stellar_energy_relative_residual']) for r in budgets)<1e-4,stellar_torque=max(abs(r['stellar_J_relative_residual']) for r in budgets)<1e-4,
        young_work=max(abs(r['young_energy_relative_residual']) for r in budgets)<1e-4,young_torque=max(abs(r['young_J_relative_residual']) for r in budgets)<1e-4,
        birth_source=max(max(abs(r[k]) for k in ['relative_mass_error','relative_J_error','relative_KE_error','vector_P_relative_residual']) for r in birth_controls)<1e-8,
        frames=len(records)==81 and len(gas_records)==321)
    checks={k:bool(v) for k,v in checks.items()};result=dict(completed=True,passed=all(checks.values()),checks=checks,wall_seconds=elapsed+time.monotonic()-start,final=final,
        scope='Numerical stability/accounting only. Gas is isothermal; imposed fields exchange work and torque with tracers. This is not live-disk stability, self-consistency or visual acceptance.')
    atomic_json(out/'run-report.json',result);print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:raise RuntimeError('Response failed declared numerical gates; preserve results and do not render as an accepted candidate')
    fcntl.flock(lock,fcntl.LOCK_UN);lock.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description='Regenerate only the selected100k/512trailing22calculation')
    p.add_argument('--output',type=Path,default=ROOT/'work/calculation')
    p.add_argument('--resume',action='store_true')
    args=p.parse_args();args.model='trailing22';args.stars=100000;args.grid=512
    run(args)
