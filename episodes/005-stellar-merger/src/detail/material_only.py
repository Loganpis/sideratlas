"""Material-following cooling display; real stellar extinction is retained.

The separate radial radiation cache is audited but omitted from this display.
Omitted energy is not reassigned to gas, stars, or new heating.
"""
import numpy as np
from src.detail.particle_cooling import ParticleRenderer

class MaterialOnlyRenderer(ParticleRenderer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        if self.micro.settings.get("material_geometry_support"):
            from src.detail.material_geometry import MaterialGeometryFlow
            self.micro=MaterialGeometryFlow(self.control_micro,self.motion)
    def emitting_layer(self,t,q,right,up,eye,L,power):
        zero=np.zeros((self.H,self.W,3),np.float32);self.reprocessed_reference=zero.copy();_,w=np.polynomial.laguerre.laggauss(3);mu=np.maximum(0,self.radiation['directions']@eye);view=float(np.sum(power*4*mu*(.4+.6*mu)/.8))
        self.emitting_metrics=dict(allocated_physical_luminosity_W=float(power.sum()),view_isotropic_equivalent_source_W=view,view_power_in_native_footprints_W=0.,last_interaction_weight_sum=float(w.sum()),reprocessing_source_samples=0,samples_per_cell_node=0,continuous_emission_displayed=False,surface_samples_are_radiation_not_material=True,persistent_source_labels=True,bounded_model_paths=True,path_segments=0,path_pixel_step=self.path_pixel_step,radiation_display_omitted=True,omitted_physical_luminosity_W=float(power.sum()))
        return zero,np.zeros((self.H,self.W),np.float32),np.zeros(len(power)*3)
    def render_linear(self,p):
        image,mask=super().render_linear(p);self.records[-1]['visual_support_kind']='local material-neighborhood covariance' if self.micro.settings.get('material_geometry_support') else 'original velocity-aligned covariance';self.records[-1]['fine_footprints']['evaluated_but_not_displayed']=True;self.records[-1].update(material_only_display=True,physical_stellar_absorption_retained=True,omitted_reprocessed_luminosity_W=self.records[-1]['rendered_photosphere_luminosity_W'],displayed_reprocessed_luminosity_W=0.)
        return image,mask
