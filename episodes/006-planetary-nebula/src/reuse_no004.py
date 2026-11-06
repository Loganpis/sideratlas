"""Actual No.004 code/data adapter; original project is read-only.

The original sampler, MURaM mapper, golden grade and optional C IV projector
are imported from their source files. Evolved-star morphology is explicitly a
pressure-scale-height-rescaled solar-radiation surrogate, NOT giant RHD.
"""
from pathlib import Path
from contextlib import contextmanager
import importlib.util,sys,os,json,copy,hashlib
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent/'004-main-sequence-star'
YROW=np.array([.2126729,.7151522,.0721750])

class No004Adapter:
    def __init__(self,spatial_average=2):
        if spatial_average not in (1,2,3):raise ValueError('Use a conservative integer source average: 1, 2 or 3')
        self.spatial_average=spatial_average;self.modules={};self.calls={};self.atmospheres={}
        self.cfg=json.loads((SOURCE/'golden_orbit/config-final.json').read_text())
        os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'cache/no004-numba'))
        with self._aliases():
            for name,rel in [('continuum','src/continuum.py'),('geometry','src/geometry.py'),('solar_surface','src/solar_surface.py'),('active_region','src/active_region.py'),('surface','photosphere/surface.py'),('optics','src/optics.py'),('video_export','src/video_export.py'),('stellar_renderer','golden/stellar_renderer.py'),('atmosphere','golden/atmosphere.py'),('selected_render','golden_orbit/render.py')]:
                spec=importlib.util.spec_from_file_location('_sideratlas_no004_'+name,SOURCE/rel)
                mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;sys.modules[name]=mod
                spec.loader.exec_module(mod);self.modules[name]=mod
            self.geometry=self.modules['geometry'];self.surface_module=self.modules['surface']
            self.grade=self.modules['selected_render'].grade_photosphere
            self.YROW=self.modules['solar_surface'].YROW.copy()
            files=sorted(p for p in (SOURCE/'data/continuum').glob('*-17-2-solar.npz') if p.with_suffix('.json').exists())
            first=np.load(files[0]);n=spatial_average
            shape=first['rgb'].shape;assert shape[2]%n==shape[3]%n==0
            self.fields=np.empty((len(files),shape[0],shape[1],shape[2]//n,shape[3]//n,3),'f4')
            self.times=np.empty(len(files))
            for i,p in enumerate(files):
                with np.load(p) as d:
                    a=d['rgb'];self.fields[i]=a.reshape(shape[0],shape[1],shape[2]//n,n,shape[3]//n,n,3).mean(axis=(3,5));self.times[i]=float(d['time_s'])
            self.mu=first['mu'];self.mean_y=(self.fields[0]@self.YROW).mean(axis=(1,2,3));self.normalization=float(self.mean_y[-1])
            self.fields/=self.normalization;self.mean_y/=self.normalization
            self.redshift=(1-2*6.67430e-11*1.98847e30/(6.957e8*299792458.**2))**.5
            quiet=SimpleNamespace(mu=self.mu,mean_mu_y=self.mean_y)
            self.active=self.surface_module.CoherentActiveRegion(self.cfg,quiet)
        self.file_paths=files
    @contextmanager
    def _aliases(self):
        names=['continuum','geometry','solar_surface','active_region','surface','optics','video_export','stellar_renderer','atmosphere','selected_render']
        old={k:sys.modules.get(k) for k in names};paths=sys.path.copy();dont=sys.dont_write_bytecode
        sys.dont_write_bytecode=True
        sys.path[:0]=[str(SOURCE/'golden'),str(SOURCE/'photosphere'),str(SOURCE),str(SOURCE/'src')]
        try:
            for k,v in self.modules.items():sys.modules[k]=v
            yield
        finally:
            sys.path[:]=paths;sys.dont_write_bytecode=dont
            for k,v in old.items():
                if v is None:sys.modules.pop(k,None)
                else:sys.modules[k]=v
    def _quiet(self,normals,mu,phase,observer,magnification=10.,contrast=3.25,delays=None):
        n=np.asarray(normals,float);mu=np.asarray(mu,float);o=np.array(observer,dtype=float,copy=True);o/=np.linalg.norm(o)
        delay=np.zeros(len(n)) if delays is None else np.asarray(delays,float)
        pt=self.cfg['source_start_s']+20*float(phase)
        if pt-delay.max(initial=0)<self.times[0] or pt+300>self.times[-1]:raise ValueError('No004 source clock outside actual cached data')
        self.calls['photosphere/surface.py:sample_patches']=self.calls.get('photosphere/surface.py:sample_patches',0)+1
        return self.surface_module.sample_patches(n,o,mu,delay,self.fields,self.times,pt,float(magnification),self.mu,self.mean_y,float(contrast))*self.redshift**4
    def source_solar(self,normals,mu,phase,observer,grade=True,active=True,delays=None):
        """Original quiet + coherent MURaM RGB, normalized as No.004.

        normals and observer are world/star-frame vectors; mu is the physical
        emission cosine. This does not render an image, enlarge radii or add gas.
        """
        rgb=self._quiet(normals,mu,phase,observer)
        if active:
            a=self.active;delay=np.zeros(len(mu)) if delays is None else np.asarray(delays,float)
            t=a.start+20*float(phase)
            if not a.times[0]<=t<=a.times[-1]:raise ValueError('No004 MURaM data coverage exceeded')
            patch,weight=self.surface_module.active_rectangle(np.asarray(normals,float),delay,t,a.center,a.east,a.north,a.intensity,a.large_scale,a.times,np.asarray(mu,float),self.mu,self.mean_y,a.colors,a.tempgrid,a.cfg['spatial_magnification'],a.quiet,a.cfg['contrast_exponent'])
            rgb=rgb*(1-weight[:,None])+patch*weight[:,None]
            self.calls['photosphere/surface.py:active_rectangle']=self.calls.get('photosphere/surface.py:active_rectangle',0)+1
        if grade:
            ref=np.interp(mu,self.mu,self.mean_y)*self.redshift**4
            rgb=self.grade(rgb,ref,self.cfg['color_grade'])
            self.calls['golden_orbit/render.py:grade_photosphere']=self.calls.get('golden_orbit/render.py:grade_photosphere',0)+1
        return rgb
    def morphology(self,normals,mu,phase,observer,radius_rsun,teff,mass_msun,contrast=1.):
        """Actual sampled quiet radiation divided by its original solar limb.

        Only dimensionless spatial morphology is reused for evolved stars.
        Caller must provide evolved spectrum, limb response and luminosity.
        H_p~T/g scale, with solar display magnification removed as cells grow.
        """
        hp_ratio=(float(teff)/5772)*float(radius_rsun)**2/float(mass_msun)
        solar_display_mag=10**(1-np.clip(np.log10(max(radius_rsun,1))/1.5,0,1))
        mag=solar_display_mag*hp_ratio/float(radius_rsun)
        rgb=self._quiet(normals,mu,phase,observer,magnification=mag,contrast=contrast)
        limb=np.interp(mu,self.mu,self.mean_y)*self.redshift**4
        ratio=np.divide(rgb@self.YROW,limb,out=np.ones_like(limb),where=limb>0)
        return ratio,dict(source='No004 actual quiet radiation + original sample_patches',H_p_ratio_to_solar=hp_ratio,turnover_time_ratio_to_solar=np.sqrt(float(teff)/5772)*float(radius_rsun)**2/float(mass_msun),effective_sampler_magnification=mag,solar_display_magnification=solar_display_mag,claim='Rescaled source morphology only; not giant radiation hydrodynamics')
    def atmosphere_image(self,phase,width,height,center_fraction,radius_fraction,camera=None,depth_samples=None):
        """Original C IV 2.5D projector, applicable only to the solar stage."""
        cfg=copy.deepcopy(self.cfg);cfg.update(width=width,height=height,center_fraction=list(center_fraction),radius_fraction_of_width=float(radius_fraction))
        if camera is not None:cfg['camera']=copy.deepcopy(camera)
        # One shared physical emissivity sequence. Changing projected radius or
        # camera must not allocate another 108MB source or turn patches to face
        # the viewer. Source geometry remains anchored to its initial solar frame.
        if not self.atmospheres:
            with self._aliases():self.atmospheres['source']=self.modules['atmosphere'].Atmosphere(cfg)
        atmosphere=self.atmospheres['source'];atmosphere.cfg=cfg;atmosphere.a=cfg['atmosphere']
        self.calls['golden/atmosphere.py:Atmosphere.image']=self.calls.get('golden/atmosphere.py:Atmosphere.image',0)+1
        return atmosphere.image(float(phase),depth_samples=depth_samples)
    def optics(self,shape,cfg):
        """Instantiate the original solar spectrum-weighted No004 Optics."""
        with self._aliases():return self.modules['optics'].Optics(shape,cfg)
    def provenance(self):
        sources=['golden_orbit/render.py','golden/stellar_renderer.py','photosphere/surface.py','src/solar_surface.py','src/active_region.py','src/geometry.py','src/optics.py','golden/atmosphere.py','src/continuum.py']
        return dict(original_root=str(SOURCE),calls=self.calls.copy(),source_average=self.spatial_average,fields_shape=list(self.fields.shape),source_time_range_s=[float(self.times[0]),float(self.times[-1])],actual_source_files=[dict(path=str(SOURCE/s),sha256=hashlib.sha256((SOURCE/s).read_bytes()).hexdigest()) for s in sources],scope='Actual original source functions and MHD-derived maps; averaging/morphology adaptation explicitly disclosed; no original files changed')
