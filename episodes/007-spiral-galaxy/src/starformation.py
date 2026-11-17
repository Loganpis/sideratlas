"""Nonnegative conservative quadrature of calculated gas sinks."""
import numpy as np

def positive_calibration(A,b):
    """Small active-set minimum-change quadrature, no negative source mass."""
    n=A.shape[1]; target=np.full(n,1/n); active=np.ones(n,bool)
    for _ in range(n-A.shape[0]):
        aa=A[:,active]; w0=target[active]
        delta=np.linalg.lstsq(aa@aa.T,b-aa@w0,rcond=1e-13)[0]
        w=w0+aa.T@delta
        if np.min(w)>=-1e-12:
            out=np.zeros(n); out[active]=np.maximum(0,w)
            residual=float(np.max(abs(A@out-b)))
            if residual<1e-9:
                return out,residual
            break
        bad=np.flatnonzero(active)[np.argmin(w)]; active[bad]=False
    return None,None

def sample_sink(gas,rng):
    dm=gas.mass_pending; mom=gas.momentum_pending
    M=float(dm.sum())
    symmetry=float(np.max(abs(dm-dm[::-1,::-1]))/max(float(dm.max()),1e-90))
    J=float(np.sum(gas.X*mom[1]-gas.Y*mom[0])); E=gas.KE_pending
    jref=max(abs(J/M),1e-8); eref=max(E/M,1e-8)
    P=mom.sum(axis=(1,2)); rs=8.; vs=.2
    centroid=np.array([np.sum(gas.X*dm),np.sum(gas.Y*dm)])/M
    b=np.array([1.,centroid[0]/rs,centroid[1]/rs,P[0]/M/vs,P[1]/M/vs,J/M/jref,E/M/eref])
    cdf=np.cumsum(dm.ravel()); cdf/=cdf[-1]
    n=max(128,int(np.ceil(M/5000)))
    for trial in range(4):
        quant=(np.arange(n)+rng.random(n))/n
        index=np.searchsorted(cdf,quant); jj,ii=np.unravel_index(index,dm.shape)
        px=gas.X[jj,ii]+rng.uniform(-.5,.5,n)*gas.dx
        py=gas.Y[jj,ii]+rng.uniform(-.5,.5,n)*gas.dx
        pos=np.c_[px,py,np.zeros(n)]
        # Calculated, time-averaged sink velocity of the selected hydro cell.
        vel=np.c_[mom[0,jj,ii]/dm[jj,ii],mom[1,jj,ii]/dm[jj,ii],np.zeros(n)]
        ji=np.cross(pos,vel)[:,2]; ei=.5*np.sum(vel*vel,axis=1)
        A=np.vstack([np.ones(n),px/rs,py/rs,vel[:,0]/vs,vel[:,1]/vs,ji/jref,ei/eref])
        w,residual=positive_calibration(A,b)
        if w is not None:
            keep=w>1e-12; pos=pos[keep]; vel=vel[keep]; w=w[keep]
            positions=pos; velocities=vel; weights=w*M
            jp=np.sum(np.cross(positions,velocities)*weights[:,None],axis=0)
            ep=float(.5*np.sum(weights*np.sum(velocities*velocities,axis=1)))
            return positions,velocities,weights,dict(
                cell_sink_point_symmetry_relative_error=symmetry,
                relative_mass_error=float((weights.sum()-M)/M),
                relative_J_error=float(abs(jp[2]-J)/max(abs(J),1e-90)),
                relative_KE_error=float(abs(ep-E)/max(E,1e-90)),
                vector_P_relative_residual=float(np.linalg.norm(np.sum(weights[:,None]*velocities,axis=0)[:2]-P)/(M*vs)),
                max_birth_velocity_correction_km_s=0.,
                weight_ratio_max_to_uniform=float(w.max()*n),trials=trial+1,
                scope='Nonnegative source quadrature matching mass, centroid, momentum, Jz, KE with exact local cell velocities; not an inferred cluster initial mass function')
        n*=2
    raise RuntimeError('Nonnegative birth quadrature failed; preserve evidence rather than change local velocities')
