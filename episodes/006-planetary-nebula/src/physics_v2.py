"""History-dependent, nonperiodic launch inhomogeneities for the No.006 prototype.

This subclasses the existing conservative wind/sector model. Spherical patches
are PRESCRIBED finite-duration launch conditions, not RHD turbulence or calculated
filament instabilities. The same allocation drives source mass, swept mass and
shell forces. No density/emissivity texture is added after the physical state.
"""
from __future__ import annotations

import numpy as np
from .physics import Nebula as OriginalNebula


class Nebula(OriginalNebula):
    """API-compatible replacement for the original ``physics.Nebula`` class.

    The inherited stellar evolution, mass history, opacity, line radiation,
    central-wind prescription and shared-pressure shell equations are retained.
    Only the angular/time distribution of the conserved launched mass changes.
    """

    simulation_end_years = 45000
    patch_seed = 6006

    def __init__(self, nmu=32, nphi=64, nr=256):
        super().__init__(nmu=nmu, nphi=nphi, nr=nr)

    def radial_outer_radius_pc(self, t):
        """Enclose every physical swept sector plus eight Gaussian shell widths.

        The original domain followed the 12 km/s wind only. Faster swept material
        can exceed it at late epochs; normalizing a truncated shell then returned
        exterior mass into the grid. This expands the REAL integration domain,
        populated only by the original conservative wind history.
        """
        if t <= 0:
            return 0.
        n = self.mu.size
        radii = self.shell_solution.sol(min(float(t), self.simulation_end_years))[:n]
        return float(np.max(radii) * (1 + 8*.075))

    def _prepare_patches(self):
        if hasattr(self, "patch_sites"):
            return
        rng = np.random.default_rng(self.patch_seed)
        # Each finite event is a coherent source patch. Broad episodes describe
        # long-lived global wind inhomogeneity; smaller patches describe assumed
        # unresolved variability, not freshly solved convection.
        specs = [
            # count, earliest birth yr, width rad, lifetime sigma yr, log amplitude
            (64, -300000., (.38, .80), (5000., 30000.), .35),
            (160, -60000., (.15, .40), (500., 3500.), .50),
            (256, -20000., (.07, .22), (150., 1100.), .55),
        ]
        sites = []; widths = []; centers = []; lifetimes = []; amplitudes = []
        for count, earliest, wrange, trange, amp in specs:
            xyz = rng.normal(size=(count, 3)); xyz /= np.linalg.norm(xyz, axis=1)[:, None]
            sites.append(xyz)
            widths.append(rng.uniform(*wrange, size=count))
            # Irregular nonperiodic event centers; no longitude harmonics and no
            # reset when the camera completes a turn.
            centers.append(rng.uniform(earliest, 0., size=count))
            lifetimes.append(np.exp(rng.uniform(np.log(trange[0]), np.log(trange[1]), size=count)))
            amplitudes.append(np.clip(rng.normal(0., amp, size=count), -1.2, 1.2))
        self.patch_sites = np.concatenate(sites)
        self.patch_widths = np.concatenate(widths)
        self.patch_centers = np.concatenate(centers)
        self.patch_lifetimes = np.concatenate(lifetimes)
        self.patch_amplitudes = np.concatenate(amplitudes)
        st = np.sqrt(np.maximum(0., 1 - self.mu.ravel()**2))
        normals = np.column_stack([st*np.cos(self.phi.ravel()), st*np.sin(self.phi.ravel()), self.mu.ravel()])
        self.patch_normals = normals
        self.patch_basis = self.spherical_basis(normals)

    def spherical_basis(self, normals):
        """Gaussian covariance on chord distance, smooth at both poles/seam."""
        dot = np.asarray(normals) @ self.patch_sites.T
        return np.exp((np.clip(dot, -1., 1.) - 1.) / self.patch_widths[None, :]**2)

    def patch_coefficients(self, birth):
        birth = np.atleast_1d(np.asarray(birth, dtype=float))
        u = (birth[None, :] - self.patch_centers[:, None]) / self.patch_lifetimes[:, None]
        # Compact C-infinity pulse: Gaussian near its center and exactly zero
        # beyond three lifetime parameters, without on/off discontinuities.
        q = u*u; window = np.zeros_like(q)
        use = q < 9.
        window[use] = np.exp(-.5*q[use]/(1 - q[use]/9))
        return self.patch_amplitudes[:, None] * window

    def clump(self, birth):
        self._prepare_patches()
        birth = np.atleast_1d(np.asarray(birth, dtype=float))
        logf = self.patch_basis @ self.patch_coefficients(birth)
        # No time-dependent renormalization of total stellar mass loss: normalize
        # the angular allocation separately for every original mass cohort.
        f = self.A.ravel()[:, None] * np.exp(logf)
        f /= f.mean(axis=0, keepdims=True)
        return f.reshape(*self.mu.shape, len(birth))

    def radial_inner_radius_pc(self,t):
        if t<=0:return 0.
        radii=self.shell_solution.sol(min(t,self.simulation_end_years))[:self.mu.size]
        # Cold shell begins more than 8.6 Gaussian widths above this boundary.
        # The absent hot-bubble optical emission was already outside the model;
        # no physical photosphere or unresolved stellar gas is discarded.
        return .35*float(radii.min())

    def compressed_shell_mass(self,edges,Rsh,t):
        """Preserve material launch history through finite shell compression.

        A prescribed ordered compression map, not simulated turbulent mixing.
        Old material lies outward, recently swept material inward. Its density
        modulation is inherited from actual source cohorts, not post-render noise.
        """
        from scipy.special import ndtr
        n=self.mu.size;rad=Rsh.reshape(n,1);width=.075*rad
        cdf=ndtr((edges[None,:]-rad)/np.maximum(width,1e-12))
        u=(cdf-cdf[:,:1])/np.maximum(cdf[:,-1:]-cdf[:,:1],1e-30)
        oldest=np.minimum(0,t-Rsh.ravel()/self.vslow.ravel())
        mass_old=self.track.massfn(oldest);m0=self.track.M0
        target=m0+(mass_old[:,None]-m0)*u
        # Invert the SAME monotone mass history. Endpoint labels are exact;
        # internal birth sampling only sets the declared compression ordering.
        masses=self.track.massfn(self.birth_edges)
        born=np.interp(target,masses[::-1],self.birth_edges[::-1])
        born=np.maximum(oldest[:,None],np.minimum(0.,born))
        born[:,0]=0.;born[:,-1]=oldest
        cumulative=self.cohort(born)
        packed=np.maximum(0,-np.diff(cumulative,axis=-1))
        return packed.reshape(*self.mu.shape,len(edges)-1)
