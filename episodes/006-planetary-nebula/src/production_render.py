"""Accepted v4 scenes, native area sampling and a short display shutter.

The original frame method produces all physical states, camera, exposures and
HDR source terms. Its unused intermediate optical images are suppressed. The
two physical-time samples are averaged in detector-linear HDR before the actual
original optical operators and display response are applied once.
"""
from pathlib import Path
import json,time
import numpy as np
from . import render_v3,render_v4
from .metal_projector import MetalProjector
from .production_numerics import stable_early_state
from .evolution import smooth,PC,RSUN
ROOT=Path(__file__).resolve().parents[1]

class ProductionRenderer(render_v4.Renderer):
    def __init__(self,spatial_samples=4,temporal_samples=2):
        assert json.loads((ROOT/'validation/production-metal-qualification.json').read_text())['passed']
        super().__init__(1080,1920)
        assert type(self.model).__module__=='src.physics_v3'
        physical_state=self.model.state
        self.model.state=lambda t:stable_early_state(self.model,physical_state,t)
        self.gpu=MetalProjector(spatial_samples)
        self.previous_projectors=(render_v3.project_material,render_v4.project_material)
        render_v3.project_material=render_v4.project_material=self.gpu
        self.solar_process=self.solar_optics.process
        self.compact_linear=self.compact_optics.linear_process
        self.sharp_linear=self.sharp_optics.linear_process
        self.zero=np.zeros((1920,1080,3))
        self.temporal_samples=temporal_samples
        self.suppress_unused_optics()

    def suppress_unused_optics(self):
        self.solar_optics.process=lambda _:self.zero
        self.compact_optics.linear_process=lambda _:self.zero
        self.sharp_optics.linear_process=lambda _:self.zero

    def legacy_reference(self,t):
        """Original v4 result for a targeted equivalence check, same GPU kernel."""
        self.solar_optics.process=self.solar_process
        self.compact_optics.linear_process=self.compact_linear
        self.sharp_optics.linear_process=self.sharp_linear
        try:return render_v4.Renderer.frame(self,t,30)
        finally:self.suppress_unused_optics()

    def old_display(self,hdr,t):
        weight=float(smooth((t-3)/13))
        if weight<1:a=self.solar_process(hdr)
        if weight>0:
            linear=self.compact_linear(hdr)
            b=linear/(1+np.max(linear,axis=-1)[...,None])
        return a if weight==0 else (b if weight==1 else (1-weight)*a+weight*b)

    def production_frame(self,index):
        begin=time.monotonic();t=(index+.5)/30
        if self.temporal_samples==1:offsets=[0.]
        else:
            assert self.temporal_samples==2
            half=1/(4*30*np.sqrt(3));offsets=[-half,half]
        old=None;new=None;records=[];sample_weight=1/len(offsets)
        calls_start=len(self.gpu.calls)
        for offset in offsets:
            sample_t=t+offset
            render_v4.Renderer.frame(self,sample_t,30)
            record=self.records[-1];a=record['accounting']
            assert a['wind_shell_mass_relative_error']<1e-10 and a['photon_relative_closure']<1e-10,(sample_t,a['wind_shell_mass_relative_error'],a['photon_relative_closure'])
            if 'photon_integral_relative_closure' in a:assert a['photon_integral_relative_closure']<1e-10
            assert a['line_luminosity_erg_s']<=a['absorbed_ionizing_luminosity_erg_s']*(1+1e-10)+1e-50
            for key in ['H1_fraction_range','O2_fraction_range','N1_fraction_range']:
                if key in a:assert min(a[key])>=0 and max(a[key])<=1
            if t<9:
                exposed=self.last_layers['hdr']*self.last_layers['exposure']
                if old is None:old=sample_weight*exposed
                else:old+=sample_weight*exposed
            if t>8:
                exposed=self.last_layers['v4_detector_input']
                if new is None:new=sample_weight*exposed
                else:new+=sample_weight*exposed
            records.append(record);self.records=[]
        if t<9:old_display=self.old_display(old,t)
        if t>8:
            linear=self.sharp_linear(new);maximum=np.max(linear,axis=-1)
            exponent=1+.35*float(smooth((t-14)/5))
            power=maximum**exponent;mapped=power/(1+power)
            new_display=linear*np.divide(mapped,maximum,out=np.zeros_like(mapped),where=maximum>0)[...,None]
            blend=float(smooth((t-8)/1))
            display=new_display if t>=9 else old_display*(1-blend)+new_display*blend
        else:display=old_display
        display=display*(1-self.alpha[...,None])+self.alpha[...,None]*.9
        srgb=np.where(display<=.0031308,12.92*display,1.055*np.maximum(0,display)**(1/2.4)-.055)
        image=np.uint16(np.clip(srgb,0,1)*65535+.5)
        physical=float(self.track.clock(t));star=self.track.state(physical)
        ledger=dict(index=index,playback_s=t,physical_age_yr=physical+self.track.t0,post_AGB_reference_yr=physical,
            star_center=dict(M=star['M'],R_Rsun=star['R']*PC/RSUN,Teff_K=star['T'],L_erg_s=star['L']),
            samples=records,render_s=time.monotonic()-begin,ray_calls=self.gpu.calls[calls_start:],
            temporal_samples=len(offsets),shutter_degrees=180 if len(offsets)>1 else 0,
            spatial_samples_unresolved=self.gpu.spatial_samples,
            physical_model='accepted v3; rejected v4/v5 disabled',imaging='accepted v4; normalized two-time HDR integration before original optics/tone')
        return image,ledger

    def close(self):
        render_v3.project_material,render_v4.project_material=self.previous_projectors
        self.gpu.close()
