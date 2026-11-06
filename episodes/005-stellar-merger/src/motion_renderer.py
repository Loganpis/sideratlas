"""Unchanged detailed photospheres over corrected, unlinked gas transport."""
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from src.quality.renderer import QualityRenderer
from src.motion_view import MotionData,mass_grid
from src.transport_volume import adaptive_grid,cooling_rays,absorbed_fraction,reprocessed_upper
from src.stellar_structure import MSUN,RSUN

class CachedScene:
    def __init__(self,motion):self.motion=motion;self.base=motion.scene;self.M=self.base.M;self.m=self.base.m;self.R=self.base.R;self.r2=self.base.r2
    def state(self,t):return self.motion.state(t).state(t)
    def centers(self,t,phase):return self.motion.state(t).centers(t,phase)

class MotionRenderer(QualityRenderer):
    def __init__(self,source,width=540,gas_width=128):
        super().__init__(source,width,gas_width)
        self.motion=MotionData(source,self.data,self.report);self.scene=CachedScene(self.motion);self.clock=self.motion.clock;self.camera_width=self.motion.camera_curve();self.tidal_response=self.data['tidal_tensor']
        # Actual primary-return events only. Accretor impacts and engulfment
        # by surface expansion must not be relabelled as fallback heating.
        ids=np.flatnonzero(np.isin(self.data['sink_kind'],[2,3]));self.impact_times=self.data['death_times'][ids]
        q=self.data['x'][-1,ids].astype('f8');vel=self.data['v'][-1,ids].astype('f8');center=np.column_stack([np.interp(self.impact_times,self.times,self.data['donor_center'][:,j]) for j in range(3)]);vd=np.column_stack([np.interp(self.impact_times,self.times,self.data['donor_velocity'][:,j]) for j in range(3)]);pos=q-center
        self.impact_directions=pos/np.maximum(np.linalg.norm(pos,axis=1,keepdims=True),1e-12);om=np.interp(self.impact_times,self.times,self.data['envelope_omega']);relative=vel-vd-om[:,None]*np.cross(np.array([0,0,1.]),pos);self.impact_energy=.5*self.motion.mass*np.sum(relative*relative,axis=1)
        self.records=[]
        # Temperature chromaticity extends into the hot-gas regime; the
        # spectral shape never sets the bolometric luminosity of a mass cell.
        from src.radiation.spectrum import H,C,K,XYZ_RGB
        cmf=np.loadtxt(Path(__file__).resolve().parents[1]/'data/spectral/CIE_xyz_1931_2deg.csv',delimiter=',')
        wave=cmf[:,0]*1e-9;temps=np.geomspace(300,1e7,1000);B=2*H*C*C/wave[None,:]**5/np.expm1(H*C/(wave[None,:]*K*temps[:,None]));xyz=np.trapz(B[:,:,None]*cmf[None,:,1:4],wave,axis=1);rgb=np.maximum(0,xyz@XYZ_RGB.T)
        B0=2*H*C*C/wave**5/np.expm1(H*C/(wave*K*self.report['donor']['Teff']));base=np.trapz(B0*cmf[:,2],wave);self.gas_lut=(temps,rgb/base)
    def render_linear(self,playback):
        f=np.clip(playback/32,0,1);t=float(self.clock(playback));q=self.motion.at(t);phase=q['phase'];s=self.scene;R=q['radius']
        # Finish pulling back BEFORE the final interval. Fixed position,
        # orientation and scale then make genuine expansion directly visible.
        from src.dynamics.prescribed import smooth
        L=6.85+(float(self.camera_width(390))-6.85)*float(smooth((t-210)/90))
        inc=.82-.12*np.clip((min(t,300)-215)/175,0,1);az=-np.pi/2+.35*min(f,.75);eye=np.array([np.sin(inc)*np.cos(az),np.sin(inc)*np.sin(az),np.cos(inc)]);right=np.array([-np.sin(az),np.cos(az),0]);up=np.cross(eye,right);center=np.zeros(3)
        lum=float(np.interp(t,self.times,self.data['gas_luminosity_code']+self.data['remnant_luminosity_code']))*self.code_lum
        star,depth,mask,Lgas=self.stars_native(t,phase,right,up,eye,L,center,lum);active=q['live'];n=int(active.sum());survival=1.;mass_in_volume=0.;upper=0.;absorbed=0.;reprocess_upper=0.;stellar_physical_luminosity=0.
        if n:
            mass,heat,origin,spacing,zc,zs,support_info=adaptive_grid(q['x'][active],np.full(n,self.motion.mass),q['T'][active],L,self.gas_width,96,eye=eye,right=right,up=up,reference_center=q['donor'],reference_radius=R)
            pixel=spacing[0];kappa=self.config['grey_opacity_cm2_g']*.1*(MSUN/1000)/(pixel*RSUN/100)**2
            cv=1.5*1.380649e-23/(1.3*1.67262192369e-27);cooling=self.config['gas_cooling_time_code']*self.report['physical_time_unit_seconds'];raw_upper=float(heat.sum(dtype='f8')*(MSUN/1000)*cv/cooling);dilution=min(1.,Lgas/max(raw_upper,1e-30));upper=raw_upper*dilution;mass_in_volume=float(mass.sum(dtype='f8'))
            source_center=q['donor']@np.column_stack([right,up,eye]);_,Tmean,_,_=self.stellar_state(t,lum);stellar_physical_luminosity=4*np.pi*(R*RSUN/100)**2*5.670374419e-8*Tmean**4
            absorbed=absorbed_fraction(mass,origin,spacing,source_center,R,kappa);rep_scale=0.
            result,surviving=cooling_rays(mass,heat,star,depth,origin,spacing,2*L/self.W,L,.64*self.H,zc,zs,kappa,np.log(self.gas_lut[0]),self.gas_lut[1],self.report['donor']['Teff'],cooling,self.config['grey_opacity_cm2_g']*.1,cv*dilution,.001,source_center,R,stellar_physical_luminosity,rep_scale,RSUN/100,0.)
            survival=float(surviving.sum(dtype='f8')/max(star.sum(dtype='f8'),1e-30))
            # No.004's outside-disk halo must follow the actually visible
            # photosphere, not the geometric silhouette of a buried star.
            visible=np.divide(surviving,star.sum(axis=-1),out=np.zeros_like(surviving),where=star.sum(axis=-1)>1e-20)
            mask*=np.clip(visible,0,1)
        else:result=star
        xx=q['x']@right;yy=q['x']@up;inframe=(abs(xx)<L)&(yy>-.72*16/9*L)&(yy<1.28*16/9*L);out=active&(q['classes']==0);back=active&(q['classes']==2)
        gas=n*self.motion.mass;m2=q['m2'];M1=self.scene.M+self.scene.m-gas-m2
        self.records.append(dict(playback_s=float(playback),t_code=t,primary_mass_msun=M1,companion_mass_msun=m2,gas_mass_msun=gas,mass_identity_error=abs(M1+m2+gas-1.68),camera_halfwidth_Rsun=L,camera_azimuth=az,camera_inclination=inc,remnant_diameter_fraction=R/L,gas_centers_in_frame_fraction=float(np.sum(active&inframe)/max(n,1)),outward_candidates=int(out.sum()),fallback=int(back.sum()),outward_candidates_in_frustum_fraction=float(np.sum(out&inframe)/max(out.sum(),1)),fallback_in_frustum_fraction=float(np.sum(back&inframe)/max(back.sum(),1)),gas_source_upper_W=float(upper),assigned_gas_power_W=float(Lgas),mass_in_volume_msun=mass_in_volume,stellar_display_flux_survival=survival,stellar_luminosity_W=stellar_physical_luminosity,stellar_absorbed_fraction=absorbed,reprocessed_stellar_upper_W=float(reprocess_upper),reconstruction='Current-position compact isotropic density estimate, no material links, velocity tubes, age-growth or density textures',gas_display='Transparent cooling-emissivity projection, 0.001 fixed radiance normalization and 0.06 asinh luminance map; true stellar extinction retained; reprocessing withheld; not optical photometry'))
        return result,mask
