"""Explicit-coalescence 32-second preview; unchanged reduced-order trajectories.

Mass and trajectories are retained; the grey emitting-layer closure is revised.
True grey optical depths and
energy-limited thermal sources are computed first. Power-law compressed
optical depth and neutral gas-radiance compression are DISPLAY choices needed
to see through intrinsically opaque gas; this is not a literal optical image.
"""

from pathlib import Path

import argparse,hashlib,json,sys,time

import numpy as np

from scipy.interpolate import PchipInterpolator

from scipy.ndimage import gaussian_filter,gaussian_filter1d

from scipy.signal import fftconvolve

from scipy.spatial import cKDTree

from scipy.integrate import cumulative_trapezoid

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

from PIL import Image,ImageDraw,ImageFont

from src.dynamics.prescribed import Scenario

from src.radiation.spectrum import make_lut,color

from src.video_export import encoder

from src.stellar_structure import G,MSUN,RSUN

from src.stellar_render_reuse import sample_six_charts,grade_photosphere,ReusedOptics,YROW

SIGMA=5.670374419e-8

def smooth(x):x=np.clip(x,0,1);return x*x*(3-2*x)

def cic(points,weight,temperature,shape,origin,spacing):
    """Conservative cloud-in-cell deposition; clipping is reported separately."""
    coords=(points-origin)/spacing-.5;lower=np.floor(coords).astype(int);frac=coords-lower
    mass=np.zeros(shape,np.float32);heat=np.zeros(shape,np.float32);flat=mass.ravel();hf=heat.ravel()
    for iz in [0,1]:
        for iy in [0,1]:
            for ix in [0,1]:
                inds=lower+np.array([ix,iy,iz]);valid=np.all((inds>=0)&(inds<np.array(shape[::-1])),axis=1)
                w=weight*np.prod(np.where(np.array([ix,iy,iz]),frac,1-frac),axis=1)
                index=(inds[valid,2]*shape[1]+inds[valid,1])*shape[2]+inds[valid,0]
                np.add.at(flat,index,w[valid]);np.add.at(hf,index,w[valid]*temperature[valid])
    return mass,heat

def fast_gaussian(field,sigmas):
    """Same finite Gaussian as ndimage, with FFT for broad separable kernels.

    Float64 transforms keep roundoff below the physical-opacity skin. No new
    physical texture, changed kernel, support cutoff or simulation is added.
    """
    result=field.astype(np.float64)
    for axis,sigma in enumerate(sigmas):
        radius=int(8*sigma+.5)
        if sigma<=3:
            result=gaussian_filter1d(result,sigma,axis=axis,mode='constant',truncate=8.)
        else:
            grid=np.arange(-radius,radius+1);kernel=np.exp(-.5*(grid/sigma)**2);kernel/=kernel.sum()
            shape=[1]*3;shape[axis]=len(kernel)
            result=fftconvolve(result,kernel.reshape(shape),mode='same',axes=(axis,))
    return np.maximum(result,0)

class StoryPreview:
    def __init__(self,source,width=540,base_width=256):
        self.source=source;self.report=json.loads((source/'reduced-model-report.json').read_text());self.config=self.report['config'];self.scene=Scenario(self.config)
        self.data=np.load(source/'reduced-merger-cache.npz');self.times=self.data['t'];self.phase=PchipInterpolator(self.times,self.data['phase'])
        self.width=width;self.height=round(width*16/9);self.w=base_width;self.h=round(base_width*16/9)
        self.pixel_ys,self.pixel_xs=np.mgrid[:self.h,:self.w];self.lut=make_lut();self.metrics=[]
        self.star_ys,self.star_xs=np.mgrid[:2*self.h,:2*self.w]
        self.clock=PchipInterpolator([0,.25,.5,.75,1],[0,140,210,300,390])
        self.surface_report=json.loads((ROOT/'cache/surface-statistics/surface-atlas-report.json').read_text())
        self.spin=PchipInterpolator(self.times,np.r_[0.,cumulative_trapezoid(self.data['envelope_omega'],self.times)])
        self.optics=ReusedOptics((self.h,self.w),exposure=1.05)
        self.grade=dict(relative_intensity_knots=[0,.45,.9,1.5,4],linear_rgb_gains_at_knots=[[1.20,1,.64],[1.17,1,.68],[1.13,1,.78],[1.06,1,.90],[1,1,1]])
        self.code_lum=G*MSUN**2/RSUN/1e7/self.report['physical_time_unit_seconds']
        self.base0=color(self.report['donor']['Teff'],self.lut)
        ids=np.arange(self.config['gas_parcels']);birth=300*(.25+.75*((ids+.5)/self.config['gas_parcels'])**(1/3))
        before=np.array([self.scene.state(t)[0]>1.05*self.scene.R for t in birth])
        self.channel=np.where(before,np.where(ids%5==0,0,1),2)
        self.births=birth;self.birth_phase=self.phase(birth)
        if 'launch_pattern_phase' in self.data.files:
            launch_phase=PchipInterpolator(self.times,self.data['launch_pattern_phase'])(birth)
            self.birth_phase=np.where(before,self.birth_phase,launch_phase)
        self.data={key:self.data[key] for key in self.data.files}
        work=self.data['available_orbital_work_code']
        self.deposition=np.maximum(0,np.gradient(work,self.times))
        # Delayed low-order equilibrium-tide response. The envelope's retained
        # acoustic response is preparation-independent and supplies the lag.
        tau=self.config['envelope_response_time_code'];tensor=np.zeros((3,3));history=[]
        for i,tt in enumerate(self.times):
            a,dis,alive,R=self.scene.state(float(tt));phi=float(self.phase(tt));e=np.array([np.cos(phi),np.sin(phi),0.])
            eps=min(.16,.75*(self.scene.m/self.scene.M)*(self.scene.R/max(a,self.scene.R))**3)*(1-dis)
            target=eps*(3*np.outer(e,e)-np.eye(3))
            if i:tensor+=(1-np.exp(-(tt-self.times[i-1])/tau))*(target-tensor)
            history.append(tensor.copy())
        self.tidal_response=np.array(history)
        ei,pi=np.nonzero(self.data['live'][:-1]&~self.data['live'][1:])
        self.impact_times=self.times[ei+1]
        positions=self.data['x'][ei+1,pi]-self.data['donor_center'][ei+1]
        self.impact_directions=positions/np.maximum(np.linalg.norm(positions,axis=1,keepdims=True),1e-12)
        donor_velocity=np.gradient(self.data['donor_center'],self.times,axis=0)
        spin_velocity=self.data['envelope_omega'][ei,None]*np.cross(np.array([0.,0.,1.]),positions)
        relative=self.data['v'][ei,pi]-donor_velocity[ei]-spin_velocity
        self.impact_energy=.5*float(self.data['parcel_mass_msun'])*np.sum(relative*relative,axis=1)
