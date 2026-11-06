"""Audit native MESA profiles before SPH mapping; cgs throughout.

Shell integrals are independent of MESA's reported scalar energy/inertia.
They do not validate the later three-dimensional particle realization.
"""
from pathlib import Path
import hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid

G = 6.67430e-8
MSUN = 1.9884098706980504e33
RSUN = 6.957e10

def read_mesa_profile(path):
    lines = Path(path).read_text().splitlines()
    header = dict(zip(lines[1].split(), lines[2].split()))
    names = lines[5].split()
    data = np.loadtxt(lines[6:])
    if data.ndim != 2 or data.shape[1] != len(names):
        raise ValueError('Unsupported MESA profile layout')
    columns = dict(zip(names, data.T))
    order = np.argsort(columns['mass'])
    columns = {k: v[order] for k, v in columns.items()}
    for k, v in header.items():
        try: header[k] = float(v.replace('D', 'e'))
        except ValueError: pass
    return header, columns

def audit_profile(path, expected_mass):
    header, c = read_mesa_profile(path)
    m = c['mass'] * MSUN
    r = 10**c['logR'] * RSUN
    rho = 10**c['logRho']
    p = 10**c['logP']
    temperature = 10**c['logT']
    inner_m = np.r_[0., m[:-1]]
    inner_r = np.r_[0., r[:-1]]
    dm = c.get('dm', m-inner_m)
    mmid = (m + inner_m) / 2
    rmid = c.get('rmid', (r+inner_r)/2/RSUN) * RSUN
    if np.any(np.diff(m) <= 0) or np.any(np.diff(r) <= 0):
        raise ValueError('Nonmonotonic mass or radius')
    # MESA uses logarithmic cell densities and face radii. Comparing their
    # volume-integrated mass tests the geometric/density mapping directly.
    volume = 4*np.pi/3 * (r-inner_r)*(r*r + r*inner_r + inner_r*inner_r)
    mass_from_density = float(np.sum(rho*volume))
    potential_energy = -float(np.sum(G*mmid/rmid*dm))
    pressure_virial = float(3*np.sum(p/rho*dm))
    surface_term = float(4*np.pi*r[-1]**3*p[-1])
    virial_residual = (pressure_virial + potential_energy - surface_term)/abs(potential_energy)
    # Exact constant-density shell integral for the axial moment of inertia.
    r2_mean = 3/5 * (r**5-inner_r**5)/(r**3-inner_r**3)
    inertia = float(2/3*np.sum(r2_mean*dm))
    x, y, z = (c[k] for k in ['x_mass_fraction_H','y_mass_fraction_He','z_mass_fraction_metals'])
    # Compare pressure at the same mass-centred locations, excluding the
    # unresolved centre and outer atmosphere from this derivative diagnostic.
    gravity_gradient = G*mmid/(4*np.pi*rmid**4)
    pressure_integral = -cumulative_trapezoid(gravity_gradient, mmid, initial=0.)
    fitted_pressure = p[-1] + pressure_integral - pressure_integral[-1]
    interior = (mmid/m[-1] > .001) & (mmid/m[-1] < .999)
    pressure_error = abs(fitted_pressure-p)/p
    core = x < .01
    # Only the contiguous central depleted region is a helium-core measure.
    first_hydrogen = np.flatnonzero(~core)
    last_core = int(first_hydrogen[0]-1) if len(first_hydrogen) else len(x)-1
    core_summary = None if last_core < 0 else dict(
        definition='Outer face of contiguous central cells with hydrogen mass fraction <0.01',
        mass_msun=float(m[last_core]/MSUN),radius_rsun=float(r[last_core]/RSUN),
        cells=last_core+1)
    internal_energy = float(np.sum(c['energy']*dm)) if 'energy' in c else None
    envelope_summary = None
    if last_core >= 0:
        envelope = np.arange(len(m)) > last_core
        env_w = -float(np.sum(G*mmid[envelope]/rmid[envelope]*dm[envelope]))
        env_u = float(np.sum(c['energy'][envelope]*dm[envelope])) if 'energy' in c else None
        envelope_summary = dict(mass_msun=float(dm[envelope].sum()/MSUN),
            gravitational_energy_erg=env_w,internal_energy_erg=env_u,
            binding_energy_erg=None if env_u is None else -(env_w+env_u),
            axial_inertia_g_cm2=float(2/3*np.sum(r2_mean[envelope]*dm[envelope])))
    checks = dict(
        finite_fields=bool(all(np.isfinite(v).all() for v in [m,r,rho,p,temperature,x,y,z])),
        positive_thermodynamics=bool(np.all(rho>0) and np.all(p>0) and np.all(temperature>0)),
        mass_target=bool(abs(m[-1]/MSUN/expected_mass-1)<1e-5),
        composition_closure=bool(np.max(abs(x+y+z-1))<1e-8),
        density_volume_mass=bool(abs(mass_from_density/m[-1]-1)<.02),
        virial=bool(abs(virial_residual)<.02),
        integrated_pressure_interior=bool(np.median(pressure_error[interior])<.05))
    return dict(
        scope='One-dimensional stellar-structure consistency only; SPH equilibrium and binary not validated',
        file=str(path),sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        cells=len(m),header=header,mass_msun=float(m[-1]/MSUN),radius_rsun=float(r[-1]/RSUN),
        age_years=header.get('star_age'),effective_temperature_k=header.get('Teff'),
        central_density_cgs=float(rho[0]),central_temperature_k=float(temperature[0]),
        central_hydrogen=float(x[0]),core=core_summary,envelope=envelope_summary,
        density_integrated_mass_msun=mass_from_density/MSUN,
        gravitational_energy_erg=potential_energy,internal_energy_erg=internal_energy,
        binding_energy_including_MESA_internal_energy_erg=None if internal_energy is None else -(potential_energy+internal_energy),
        energy_reference_note='Uses native MESA EOS internal energy reference, including its ionization terms; not a newly calculated radiation release',
        axial_inertia_g_cm2=inertia,gyration_coefficient=inertia/(m[-1]*r[-1]**2),
        virial_residual=virial_residual,
        median_interior_integrated_pressure_error=float(np.median(pressure_error[interior])),
        percentile95_interior_integrated_pressure_error=float(np.percentile(pressure_error[interior],95)),
        max_composition_closure_error=float(np.max(abs(x+y+z-1))),checks=checks,
        passed_1d_checks=all(checks.values()))

def compare_profiles(coarse_path, fine_path):
    _, a = read_mesa_profile(coarse_path)
    _, b = read_mesa_profile(fine_path)
    q = np.linspace(.001,.999,999)
    result = {}
    for name in ['logR','logRho','logP','logT']:
        av = np.interp(q,a['mass']/a['mass'][-1],a[name])
        bv = np.interp(q,b['mass']/b['mass'][-1],b[name])
        error = abs(10**(av-bv)-1)
        result[name] = dict(median_relative_difference=float(np.median(error)),
                            percentile95_relative_difference=float(np.percentile(error,95)))
    ac=audit_profile(coarse_path,float(a['mass'][-1]))
    bc=audit_profile(fine_path,float(b['mass'][-1]))
    scalars={}
    for key in ['radius_rsun','age_years','effective_temperature_k','central_density_cgs',
                'axial_inertia_g_cm2','gravitational_energy_erg','internal_energy_erg',
                'binding_energy_including_MESA_internal_energy_erg']:
        scalars[key]=dict(base=ac[key],fine=bc[key],fine_over_base_minus_one=float(bc[key]/ac[key]-1))
    for group,keys in [('core',['mass_msun','radius_rsun']),('envelope',['mass_msun','binding_energy_erg'])]:
        if ac[group] is None or bc[group] is None:continue
        for key in keys:
            av,bv=ac[group][key],bc[group][key]
            scalars[f'{group}.{key}']=dict(base=av,fine=bv,fine_over_base_minus_one=float(bv/av-1))
    return dict(scope='Mass-coordinate comparison of two 1D evolutionary models, not SPH resolution convergence',
                endpoint_note='Actual first-radius-crossing endpoints and corresponding coeval companion ages; no forced radius rescaling',
                fields=result,scalars=scalars)
