"""Serve the exact selected atmospheric samples without a raw-grid download.

The original Optics implementation consumes these same 65 interpolated flux
values. The full ResearchAtmosphere implementation remains in mode_model.py for
independent rebuilding from provider inputs; this adapter changes no radiation,
limb-darkening, mode, or optical equation.
"""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class ResearchAtmosphere:
    def __init__(self):
        self.record = json.loads((ROOT / 'data/final-optics-input.json').read_text())
        self.cmf = np.loadtxt(ROOT / 'data/CIE_xyz_1931_2deg.csv', delimiter=',')

    def flux_spectrum(self, temperature, logg, wavelength_nm):
        r = self.record
        wavelength = np.asarray(wavelength_nm, dtype=float)
        if (abs(float(temperature) - r['temperature_K']) > 1e-9
                or abs(float(logg) - r['logg_cgs']) > 1e-12
                or wavelength.shape != (65,)
                or not np.allclose(wavelength, r['wavelength_nm'], rtol=0, atol=1e-10)):
            raise ValueError('This cache serves only the accepted model/PSF grid; rebuild from the documented provider inputs for another model')
        return np.asarray(r['flux_spectrum'], dtype=float).copy()
