"""Temperature-ordered gas colors anchored to the launch photosphere.

No gas temperature, energy, source power, mass or trajectory is changed.
The Planckian display-temperature compression and No.004 grade are explicit
visualization operations, not an optical gas spectrum prediction.
"""
import numpy as np
from src.detail.consistent_absorption import ConsistentAbsorptionRenderer
from src.detail.fields import Y
from src.stellar_render_reuse import grade_photosphere
from src.detail.microflow import ROOT

def display_temperature(T,T0,hot_exponent=.7):
 T=np.maximum(np.asarray(T),30.);T0=np.asarray(T0);r=T/T0;upper=np.maximum(6500.,T0+200.)
 return np.where(r>=1,T0+(upper-T0)*(-np.expm1(-hot_exponent*np.log(np.maximum(r,1)))),np.maximum(800.,T0*np.minimum(r,1)**.7))

def reference_spectral_color(renderer,T,T0=None):
 Td=display_temperature(T,renderer.report['donor']['Teff'] if T0 is None else T0,renderer.micro.settings.get('gas_hot_compression_exponent',.7));raw=np.stack([np.interp(np.log(Td),np.log(renderer.gas_lut[0]),renderer.gas_lut[1][:,j]) for j in range(3)],-1);rgb=raw;return rgb/np.maximum(rgb@Y,1e-100)[...,None]

class LaunchColorRenderer(ConsistentAbsorptionRenderer):
 def __init__(self,*args,**kwargs):
  super().__init__(*args,**kwargs);a=np.load(ROOT/self.micro.settings['launch_color_reference']);self.launch_reference={k:a[k] for k in a.files};self.anchor_correction=self.launch_reference['chromaticity']/reference_spectral_color(self,self.launch_reference['temperature_K'],self.launch_reference['temperature_K'])
 def gas_chromaticity(self,T,ids,t,q):
  c=reference_spectral_color(self,T,self.launch_reference['temperature_K'][ids])*self.anchor_correction[ids];return c/np.maximum(c@Y,1e-100)[:,None]
 def gas_display_mapping(self,rgb):
  # Same luminance/contrast map as the old palette; its 75% brightness-hue
  # paint is removed. Color now comes from thermal state + source reference.
  y=np.maximum(rgb@Y,0);ref=float(np.interp(np.log(self.report['donor']['Teff']),np.log(self.gas_lut[0]),self.gas_lut[1]@Y));cfg=dict(self.grade);cfg['linear_rgb_gains_at_knots']=(np.asarray(cfg['linear_rgb_gains_at_knots'])**self.micro.settings.get('gas_grade_warmth_exponent',1.)).tolist();graded=grade_photosphere(rgb,ref,cfg);shown=.20*(y/.20)**1.3;return graded/np.maximum(y,1e-100)[...,None]*shown[...,None]
 def render_linear(self,p):
  im,mask=super().render_linear(p);self.records[-1].update(warm_grade='Source-local physical spectrum; monotone compressed temperature; No.004 luminance-preserving warm grade; not measured temperature',gas_palette='No.004 local physical launch-spectrum anchor; compressed temperature; shared stellar warm grade and luminance-only contrast',temperature_data_unchanged=True,source_cooling_power_unchanged=True,old_brightness_hue_palette_removed=True);return im,mask
