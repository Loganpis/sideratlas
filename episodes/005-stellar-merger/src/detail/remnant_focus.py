"""Remnant-focused camera and launch-relative thermal false colors.

Only the external camera and chromaticity change. The cached orbit, parcel
positions, velocities, thermal histories, stellar radii, opacity, mass and
scalar cooling sources remain unchanged. Colors are NOT absolute thermometry:
they show current temperature relative to a fixed per-parcel launch reference.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from src.detail.launch_color import LaunchColorRenderer
from src.detail.fields import Y


def thermal_display_temperature(T,photosphere,birth):
    """Monotone in T for each fixed source; no time-dependent renormalization.

    Cold photospheric launches retain the old T/Tphotosphere mapping. Hot
    launches instead show fractional cooling from max(Tbirth,Tphotosphere).
    This exposes genuine cooling hidden by the old 6000-K display saturation.
    The same color need not mean the same absolute T for different parents.
    """
    T=np.maximum(np.asarray(T),30.);P=np.asarray(photosphere);B=np.maximum(np.asarray(birth),P)
    r=T/B;weight=np.clip(np.log(B/P)/np.log(20.),0,1);weight=weight*weight*(3-2*weight)
    exponent=.7*(1-weight)+.16*weight
    cool=np.maximum(1000.,P*np.minimum(r,1.)**exponent)
    upper=np.maximum(7000.,P+300.)
    hot=P+(upper-P)*(-np.expm1(-.35*np.log(np.maximum(r,1.))))
    return np.where(r<=1,cool,hot)


class RemnantFocusRenderer(LaunchColorRenderer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        # Preserve early framing exactly, then stop trying to contain every
        # distant ejecta particle. Fixed scale after 24 s exposes expansion.
        self.focus_knots=np.array([0.,8.,14.,16.,18.,20.,22.,24.,32.])
        early=[super(RemnantFocusRenderer,self).camera_geometry(p)[0] for p in self.focus_knots[:4]]
        self.focus_width=PchipInterpolator(self.focus_knots,np.log(early+[11.7,22.,42.,55.,55.]))
        self.palette_mode=self.micro.settings.get('thermal_palette','launch_relative')
        self.focus_camera=self.micro.settings.get('remnant_focused_camera',True)

    def camera_geometry(self,playback,t=None):
        L,inc,az,camera=super().camera_geometry(playback,t)
        if self.focus_camera and playback>16.:
            target=float(np.exp(self.focus_width(np.clip(playback,0,32))))
            u=np.clip((playback-16.)/2.,0,1);blend=u**3*(10+u*(-15+6*u))
            L=(1-blend)*L+blend*target
        return L,inc,az,camera

    def gas_chromaticity(self,T,ids,t,q):
        if self.palette_mode=='legacy':return super().gas_chromaticity(T,ids,t,q)
        ref=self.launch_reference
        Td=thermal_display_temperature(T,ref['temperature_K'][ids],ref['birth_gas_temperature_K'][ids])
        raw=np.stack([np.interp(np.log(Td),np.log(self.gas_lut[0]),self.gas_lut[1][:,c]) for c in range(3)],-1)
        c=raw*self.anchor_correction[ids]
        return c/np.maximum(c@Y,1e-100)[:,None]

    def render_linear(self,p):
        im,mask=super().render_linear(p);a=self.records[-1]
        a.update(remnant_focused_camera=self.focus_camera,camera_independent_of_dynamics=True,
                 gas_palette_mode=self.palette_mode,absolute_temperature_color_scale=False if self.palette_mode!='legacy' else None,
                 warm_grade='Fixed launch-relative thermal false colors; actual T and bolometric cooling unchanged; not absolute-temperature thermometry',
                 no_new_mass_forces_or_heating=True)
        return im,mask
