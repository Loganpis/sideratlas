"""Execute applicable original No.003 code with new MIST hot-star inputs.

No WDEC structure, ZZ Ceti mode, cellular model, DA spectral cache or temperature
palette is used. Original callable bodies are AST-loaded from byte-verified
published source, avoiding imports that initialize the old pulsation model.
"""
from pathlib import Path
from functools import lru_cache
import ast,hashlib,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
from scipy.fft import rfft2,irfft2,rfftfreq,fftfreq
from .evolution import Track,PC,MSUN,YR

ROOT=Path(__file__).resolve().parents[1]
VENDOR=ROOT/'vendor/no003'
PROVENANCE=json.loads((VENDOR/'PROVENANCE.json').read_text())
CALL_COUNTS={}
LOADED={}
H=6.62607015e-34; C=299792458.;K=1.380649e-23;G=6.67430e-11

for record in PROVENANCE['files']:
    path=VENDOR/record['path']
    if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
        raise ValueError(f'Modified original No.003 source: {record["path"]}')

def _load(relative,names,env):
    path=VENDOR/relative;source=path.read_text();parsed=ast.parse(source,str(path))
    selected=[n for n in parsed.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    if {n.name for n in selected}!=set(names):raise ValueError('Missing original symbols')
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),env)
    for node in selected:
        segment=ast.get_source_segment(source,node)
        LOADED[f'{relative}:{node.name}']=dict(path=str(path),line=node.lineno,
            original_commit=PROVENANCE['commit'],body_sha256=hashlib.sha256(segment.encode()).hexdigest())
    return {n:env[n] for n in names}

_env=dict(np=np,leggauss=leggauss,RegularGridInterpolator=RegularGridInterpolator,
          Path=Path,json=json,hashlib=hashlib,ROOT=VENDOR,C=C)
_original=_load('src/white_dwarf.py',['ray_map','limb_law','limb_integral','DAAtmosphere'],_env)
_geometry=_load('src/visualization.py',['circle_rectangle_area','orientation','camera_direction','linear_to_srgb','srgb_to_linear'],dict(np=np))
circle_rectangle_area=_geometry['circle_rectangle_area'];orientation=_geometry['orientation']
limb_law=_original['limb_law'];limb_integral=_original['limb_integral']

class HotGreyAtmosphere(_original['DAAtmosphere']):
    """NEW Planck flux + grey angular input, not an extrapolated DA atmosphere.

    Spectral angular normalization, g^5 wavelength redshift, CIE integration and
    matrix conversion execute the ORIGINAL DAAtmosphere.colour_radiance method.
    Pure-H non-LTE line blanketing and a post-AGB atmosphere are NOT solved.
    """
    def __init__(self,temperature,limb_u=.6):
        self.temperature=float(temperature);self.u=float(limb_u)
        self.cmf=np.loadtxt(ROOT/'data/spectral/CIE_xyz_1931_2deg.csv',delimiter=',')
    def flux_spectrum(self,temperature,logg,wavelength_nm):
        w=np.asarray(wavelength_nm)*1e-9
        return np.pi*2*H*C*C/w**5/np.expm1(H*C/(w*K*temperature))*1e-9
    def coefficients(self,temperature,logg,wavelength_nm):
        out=np.zeros((len(wavelength_nm),4));out[:,1]=self.u
        return out


def parameters(star):
    mass=float(star['M'])*MSUN/1000;radius=float(star['R'])*PC/100
    u=2*G*mass/(radius*C*C);g=np.sqrt(1-u)
    return dict(mass_kg=mass,radius_m=radius,compactness=u,redshift_factor=g,
                logg_cgs=float(np.log10(100*G*mass/radius**2)),temperature_k=float(star['T']))

@lru_cache(maxsize=64)
def _radiance_lut(T,logg,g,limb_u):
    atmo=HotGreyAtmosphere(T,limb_u)
    mu=np.linspace(0,1,129)**2
    CALL_COUNTS['original_DAAtmosphere.colour_radiance']=CALL_COUNTS.get('original_DAAtmosphere.colour_radiance',0)+1
    rgb,xyz=atmo.colour_radiance(T,logg,mu,redshift_factor=g,step_nm=.5)
    # Original method returns photometric 683*integral in SI. No.006 uses
    # CIE-weighted cgs energy radiance without the 683 lm/W factor.
    return mu,rgb*(1000/683),xyz*(1000/683)

@lru_cache(maxsize=64)
def _rays(u):
    nodes=np.linspace(0,1,257)
    CALL_COUNTS['original_ray_map']=CALL_COUNTS.get('original_ray_map',0)+1
    cp,delay=_original['ray_map'](nodes,u,quadrature=96)
    return nodes,cp,delay

def direct_ray_samples(mu,star):
    p=parameters(star);nodes,cp,dt=_rays(p['compactness'])
    return np.interp(mu,nodes,cp),np.interp(mu,nodes,dt)*p['radius_m']/C

def compact_surface(normals,mu,star,y_reference=None,return_diagnostics=False):
    """Return RGB in the No.006 Track.Yref normalization (default).

    star: M in Msun, R in parsecs, T in K, L in erg/s. Radius is the supplied
    evolutionary radius, never the No.003 cold mass-radius relation. `normals`
    are used only as geometry here: no imported pulsations or texture is applied.
    """
    normals=np.asarray(normals);mu=np.clip(np.asarray(mu),0,1)
    if normals.shape!=mu.shape+(3,):raise ValueError('Matching normals and emission mu required')
    p=parameters(star);nodes,rgb,xyz=_radiance_lut(p['temperature_k'],p['logg_cgs'],p['redshift_factor'],.6)
    result=np.stack([np.interp(mu,nodes,rgb[:,i]) for i in range(3)],axis=-1)
    cp,delay=direct_ray_samples(mu,star)
    if y_reference is None:y_reference=Track().Yref
    result=result/float(y_reference)
    if return_diagnostics:return result,dict(**p,cos_psi=cp,relative_delay_seconds=delay,
        atmosphere='New Planck/grey input; original No.003 spectral integration; NOT hot DA/non-LTE atmosphere')
    return result

# The original Optics constructor is preserved. Its two model factories are
# explicitly replaced, so it never instantiates InteriorModes or ResearchAtmosphere.
def _optics_class(star):
    # A scalar temperature supplies unredshifted virtual-PSF weights only.
    # No mass/radius or DA-atmosphere inference is made from that scalar.
    p=parameters(star) if isinstance(star,dict) else dict(temperature_k=float(star),logg_cgs=0.,redshift_factor=1.)
    class HotStar:
        temperature_k=p['temperature_k'];logg_cgs=p['logg_cgs'];redshift_factor=p['redshift_factor']
    class StarHolder:
        star=HotStar()
    env=dict(np=np,rfft2=rfft2,irfft2=irfft2,rfftfreq=rfftfreq,fftfreq=fftfreq,
             InteriorModes=StarHolder,ResearchAtmosphere=lambda:HotGreyAtmosphere(p['temperature_k']))
    funcs=_load('luminous/src/luminous_render.py',['airy_otf','Optics'],env)
    matrix=np.array([[3.24096994,-1.53738318,-.49861076],[-.96924364,1.8759675,.04155506],[.05563008,-.20397696,1.05697151]])
    env['YROW']=np.linalg.inv(matrix)[1]
    return funcs['Optics'],env['YROW']

class CompactOptics:
    """Original No.003 normalized Airy/jitter/Moffat optics, no mode coloring.

    The original constructor/process execute directly. Detector sampling imposes
    a minimum 2.5-pixel Airy first zero and .4-pixel jitter below native width,
    limiting sampled PSF ringing for an unresolved source; the original code retains its nonnegative FFT treatment. At1080,
    the original first zero3.8 and jitter.6 are retained. This is a virtual PSF,
    not an angularly calibrated telescope and not an enlarged stellar radius.
    """
    def __init__(self,shape,star,exposure=1.,scatter_fraction=.12,preview_floor_px=2.5):
        cls,self.YROW=_optics_class(star)
        h,w=shape;s=w/1080
        self.cfg=dict(padding_at_1080=384,airy_first_zero_px_at_550nm_1080=max(3.8,preview_floor_px/s),
            jitter_sigma_px_at_1080=max(.6,.4/s),moffat_alpha_px_at_1080=85,moffat_beta=2.5,
            scatter_fraction=scatter_fraction,false_color_strength=0.,temperature_scale_K=160.,
            asinh_softening_K=4.,display_exposure=exposure)
        CALL_COUNTS['original_Optics.__init__']=CALL_COUNTS.get('original_Optics.__init__',0)+1
        self.original=cls(shape,self.cfg,temperature=star['T'] if isinstance(star,dict) else float(star),psf_nodes=33)
        self.palette=np.ones((2,3))
        if preview_floor_px<2.5:self.pixel_integrated_psf(star)
    def pixel_integrated_psf(self,star):
        """Positive detector-integrated Airy kernel for under-Nyquist previews.

        Sampling a continuous MTF only up to detector Nyquist discarded power
        above Nyquist and introduced a tiny axial ringing cross around a bright
        unresolved source. Integrate a physical nonnegative PSF over each pixel,
        then transform that sampled kernel. No gas geometry or absorption changes.
        """
        from scipy.special import j1
        from scipy.ndimage import gaussian_filter
        o=self.original;hh,ww=o.shape
        yy=np.fft.fftfreq(hh)*hh;xx=np.fft.fftfreq(ww)*ww
        Y,X=np.meshgrid(yy,xx,indexing='ij');radius=np.hypot(X,Y)
        wave=o.wave;T=star['T'] if isinstance(star,dict) else float(star)
        at=HotGreyAtmosphere(T);weight=at.flux_spectrum(T,8.,wave)*np.interp(wave,at.cmf[:,0],at.cmf[:,2])
        weight[[0,-1]]*=.5;weight/=weight.sum()
        first=self.cfg['airy_first_zero_px_at_550nm_1080']*o.scale
        kernel=np.zeros_like(radius);near=radius<24
        def airy(rr,zero):
            z=3.8317059702075125*rr/zero
            a=np.divide(2*j1(z),z,out=np.ones_like(z),where=z!=0)
            # Unit integral of the continuous Airy PSF; chromatic width changes
            # its central intensity, not its integrated spectral weight.
            return a*a*3.8317059702075125**2/(4*np.pi*zero**2)
        offsets=(np.arange(8)+.5)/8-.5
        nx=X[near];ny=Y[near]
        for lam,wgt in zip(wave,weight):
            zero=first*lam/550
            kernel+=wgt*airy(radius,zero)
            central=np.zeros(len(nx))
            for dy in offsets:
                for dx in offsets:central+=airy(np.hypot(nx+dx,ny+dy),zero)/64
            kernel[near]+=wgt*(central-airy(radius[near],zero))
        kernel=gaussian_filter(kernel,self.cfg['jitter_sigma_px_at_1080']*o.scale,mode='wrap')
        kernel=np.maximum(kernel,0);kernel/=kernel.sum()
        o.core_otf=rfft2(kernel,workers=2)
        alpha=self.cfg['moffat_alpha_px_at_1080']*o.scale
        wing=(1+(radius/alpha)**2)**(-self.cfg['moffat_beta']);wing/=wing.sum()
        o.wing_otf=rfft2(wing,workers=2)
        self.sampling_check=dict(method='8x8 pixel-integrated Airy inside24pixels, center-sampled asymptotic wings; nonnegative normalized Moffat',core_sum=float(kernel.sum()),core_min=float(kernel.min()),wing_sum=float(wing.sum()),preview_floor_px=first,scope='new detector-sampling adaptation; original No003 convolve methods retained')
    def process(self,rgb):
        rgb=np.asarray(rgb)
        if rgb.shape[:2]!=(self.original.h,self.original.w):raise ValueError('Shape mismatch')
        if not np.any(rgb>0):return np.zeros_like(rgb)
        source=np.concatenate((rgb,(rgb@self.YROW)[...,None],np.zeros(rgb.shape[:2]+(1,))),axis=-1)
        CALL_COUNTS['original_Optics.process']=CALL_COUNTS.get('original_Optics.process',0)+1
        return self.original.process(source,self.palette,0.)
    def linear_process(self,rgb):
        """Original No003 FFT/PSF methods, before its old channel-wise display map.

        Nebular line composites need a hue-preserving display response. Both
        optical components still redistribute the same input; no new halo power.
        """
        rgb=np.asarray(rgb);o=self.original;eta=self.cfg['scatter_fraction']
        core=np.stack([o.convolve(rgb[...,i],o.core_otf) for i in range(3)],axis=-1)
        halo=np.stack([o.convolve(rgb[...,i],o.wing_otf) for i in range(3)],axis=-1)
        hdr=(1-eta)*np.maximum(core,0)+eta*np.maximum(halo,0)
        original=float(np.sum(rgb@self.YROW));total=float(np.sum(hdr@self.YROW))
        if original>0:o.full_flux_error.append(abs(total/original-1))
        p=o.pad;crop=hdr[p:p+o.h,p:p+o.w]
        if total>0:o.crop_fractions.append(float(np.sum(crop@self.YROW)/total))
        CALL_COUNTS['original_Optics.convolve']=CALL_COUNTS.get('original_Optics.convolve',0)+6
        return crop
    @property
    def flux_errors(self):return self.original.full_flux_error
    @property
    def crop_fractions(self):return self.original.crop_fractions

def reuse_evidence():
    return dict(repository=PROVENANCE['repository'],commit=PROVENANCE['commit'],
                source_functions=LOADED,executed_call_counts=dict(CALL_COUNTS),
                old_physics_not_loaded=['InteriorModes','WDEC structure','DA flux interpolation/cache','ZZ Ceti modes','cellular_convection','temperature false-color overlay'])


# Actual original packed ray samples and tiled spherical projector. Only the
# unused mode tuple is replaced by empty input; its temperature method is never
# invoked. Source temperatures always come from the supplied evolutionary star.
_tile_env=dict(np=np,C=C,MODES=(),circle_rectangle_area=circle_rectangle_area)
_tiles=_load('src/visualization.py',['Samples','Tile'],_tile_env)

class CompactAdapter:
    """Public integration API for the original No.003 hot-star components."""
    def __init__(self,cmf_path=None,y_reference=None):
        expected=ROOT/'data/spectral/CIE_xyz_1931_2deg.csv'
        if cmf_path is not None:
            actual=Path(cmf_path)
            if hashlib.sha256(actual.read_bytes()).digest()!=hashlib.sha256(expected.read_bytes()).digest():
                raise ValueError('A changed CMF needs a separately fingerprinted calibration')
        self.y_reference=Track().Yref if y_reference is None else float(y_reference)
    def surface(self,mu_array,star):
        mu=np.asarray(mu_array);normals=np.stack((np.sqrt(np.maximum(0,1-mu*mu)),np.zeros_like(mu),mu),axis=-1)
        return compact_surface(normals,mu,star,y_reference=self.y_reference)
    def render_disk(self,width,height,pixel_radius,star,center=None):
        """HDR RGB disk from ORIGINAL Samples, Tile and Tile.project bodies.

        `pixel_radius` is the projected physical-radius scale; No.003's geometry
        applies R_app/R=1/g. This is an explicit camera scale, not a changed star.
        Eight-by-eight limb sampling and exact circle/pixel coverage are retained.
        """
        p=parameters(star)
        class StarGeometry:
            radius_m=p['radius_m'];redshift_factor=p['redshift_factor']
        center=(width*.5,height*.5) if center is None else center
        cfg=dict(width=int(width),height=int(height),supersample=1,limb_supersample=8,
                 center_fraction=[center[0]/width,center[1]/height],radius_fraction_of_width=pixel_radius/width)
        nodes,cp,dt=_rays(p['compactness']);ray_grid=(nodes,cp,dt)
        result=np.zeros((height,width,3),'f4')
        radius=pixel_radius/p['redshift_factor']
        for row0 in range(0,height,48):
            row1=min(height,row0+48)
            if row1<center[1]-radius-1 or row0>center[1]+radius+1:continue
            tile=_tiles['Tile'](row0,row1,cfg,StarGeometry(),ray_grid)
            regular=self.surface(tile.regular.mu,star);edge=self.surface(tile.edge_samples.mu,star)
            # Uniform-source tangent limit in extremely thin limb slivers;
            # avoids the historical fallback round-off at circle intersection.
            fallback=self.surface(np.zeros(int(tile.fallback.sum())),star)
            result[row0:row1]=tile.project(regular,edge,fallback)
            CALL_COUNTS['original_Samples.__init__']=CALL_COUNTS.get('original_Samples.__init__',0)+3
            CALL_COUNTS['original_Tile.__init__']=CALL_COUNTS.get('original_Tile.__init__',0)+1
            CALL_COUNTS['original_Tile.project']=CALL_COUNTS.get('original_Tile.project',0)+1
        return result
