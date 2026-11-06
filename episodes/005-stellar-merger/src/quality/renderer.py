"""Native-resolution photospheres and 3D gas rays over retained trajectories.

Surface morphology is a rescaled, dimensionless No.004 solar-atmosphere
surrogate. Illustrative activity and statistical gas clumping are NOT calculated
magnetic fields or a new hydro simulation. MESA spectrum is recalculated.
"""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.special import roots_laguerre
from scipy.ndimage import gaussian_filter
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from src.rendering.photosphere_base import StoryPreview,smooth,SIGMA
from src.radiation.spectrum import color
from src.stellar_structure import G,MSUN,RSUN
from src.stellar_render_reuse import ReusedOptics,grade_photosphere,YROW
from src.quality.kernels import surface_sample,deposit_packets,raycast,bolometric_upper

class QualityOptics(ReusedOptics):
    """No.004 outside-disk halo and smooth common-channel shoulder."""
    def __init__(self,shape):
        super().__init__(shape,exposure=1.6)
    def process_masked(self,rgb,mask):
        eta=.035;core=np.stack([self.convolve(rgb[...,c],self.core) for c in range(3)],-1)
        halo=np.stack([self.convolve(rgb[...,c],self.wing) for c in range(3)],-1)
        physical=(1-eta)*core+eta*halo
        self.flux_errors.append(float(abs(physical.sum(dtype=np.float64)/rgb.sum(dtype=np.float64)-1)))
        p=self.pad;self.crop_fractions.append(float(physical[p:p+self.h,p:p+self.w].sum()/physical.sum()))
        outside=1-gaussian_filter(np.pad(mask.astype('f4'),p),max(1.,self.w/180))
        hdr=np.maximum((physical+eta*2.*halo*outside[...,None])[p:p+self.h,p:p+self.w],0)
        a=1.6*hdr
        return a/(1+a.max(axis=-1))[...,None]

class QualityRenderer(StoryPreview):
    def __init__(self,source,width=1080,gas_width=256):
        super().__init__(source,width,gas_width)
        self.gas_width=gas_width;self.H=round(width*16/9);self.W=width
        q=ROOT/'data/quality-surface';a=np.load(q/'morphology.npz');self.morph=a['z'];self.mt=a['time'];self.mmu=a['mu']
        self.meso=np.load(q/'mesoscales.npz')['z'];a=np.load(q/'illustrative-active-region.npz');self.active=a['ratio'];self.active_time=a['time']
        self.meta=json.loads((q/'provenance.json').read_text())
        # Explicit unresolved-density closure: a single seeded Gaussian field
        # with a band-limited Matérn-like spectrum, frozen in material space.
        rng=np.random.default_rng(5005);n=48;k=np.fft.fftfreq(n)*n
        kz,ky,kx=np.meshgrid(k,k,k,indexing='ij');k2=kx*kx+ky*ky+kz*kz
        filt=(1+k2/4.)**(-11/12)*np.exp(-k2/(2*18**2));filt[0,0,0]=0
        field=np.fft.ifftn(np.fft.fftn(rng.normal(size=(n,n,n)))*filt).real
        self.clump_field=((field-field.mean())/field.std()).astype('f4')
        self.spot_center=np.array([.15,-.87,.47]);self.spot_center/=np.linalg.norm(self.spot_center)
        self.spot_east=np.cross(np.array([0.,0.,1.]),self.spot_center);self.spot_east/=np.linalg.norm(self.spot_east);self.spot_north=np.cross(self.spot_center,self.spot_east)
        self.optics=QualityOptics((self.H,self.W));self.records=[]
        # Formal LTE grey transfer replaces the Eddington-Barbier shortcut.
        self.atmT=np.geomspace(1800,25000,256);self.atmMu=np.linspace(0,1,65)
        u,w=roots_laguerre(12);radiation=np.zeros((256,65,3))
        for node,weight in zip(u,w):radiation+=weight*color(self.atmT[:,None]*(.75*(self.atmMu[None,:]*node+2/3))**.25,self.lut)
        self.atmRGB=radiation
        q,w=np.polynomial.legendre.leggauss(80);q=(q+1)/2;w=w/2
        f=(1-q)**6*(1+6*q+35*q*q/3)
        self.support=float(np.sqrt(3*np.sum(w*q*q*f)/np.sum(w*q**4*f)))
        self.title_core,self.title_outline=self.make_title()
        self.grade=dict(relative_intensity_knots=[0,.25,.6,1.,2.,5.],linear_rgb_gains_at_knots=[[1.30,1,.62],[1.25,1,.66],[1.14,1,.76],[1.07,1,.86],[1.015,1,.97],[1,1,1]])
    def make_title(self):
        from PIL import Image,ImageDraw,ImageFont
        scale=self.W/1080;ss=3;core=Image.new('L',(ss*self.W,ss*self.H));outline=Image.new('L',core.size)
        d=ImageDraw.Draw(core);do=ImageDraw.Draw(outline)
        for x,y,label,size,font,anchor in [(.08*self.W,.16*self.H,'@sideratlas',35,'lmroman10-regular.otf','lt'),(.92*self.W,.16*self.H,'NO.005',35,'lmroman10-regular.otf','rt'),(.5*self.W,.205*self.H,'THE STELLAR MERGER',63,'lmroman12-bold.otf','mt')]:
            ft=ImageFont.truetype(str(ROOT/'fonts'/font),round(size*scale*ss));point=(x*ss,y*ss)
            d.text(point,label,font=ft,anchor=anchor,fill=255)
            do.text(point,label,font=ft,anchor=anchor,fill=255,stroke_width=max(1,round(1.2*scale*ss)))
        return np.asarray(core.resize((self.W,self.H),Image.Resampling.LANCZOS),'f4')/255,np.asarray(outline.resize((self.W,self.H),Image.Resampling.LANCZOS),'f4')/255
    def atmosphere(self,T,mu):
        x=(np.log(np.maximum(T,1800))-np.log(self.atmT[0]))/np.log(self.atmT[1]/self.atmT[0]);i=np.clip(x.astype(int),0,len(self.atmT)-2);f=np.clip(x-i,0,1)
        y=np.clip(mu,0,1)*64;j=np.minimum(y.astype(int),63);g=y-j
        return (1-f[:,None])*((1-g[:,None])*self.atmRGB[i,j]+g[:,None]*self.atmRGB[i,j+1])+f[:,None]*((1-g[:,None])*self.atmRGB[i+1,j]+g[:,None]*self.atmRGB[i+1,j+1])
    def thermal_pattern(self,t,normals):
        # Positive full angular kernels retain localized heating; unlike the
        # former dipole reduction they do not average away every impact patch.
        end=min(t,300.);sel=(self.times<=end)&(self.times>end-35.)
        q=self.times[sel];buried=np.clip((1.05*self.scene.R-self.data['separation'][sel])/(.65*self.scene.R),0,1)
        age=end-q;tau=self.config['envelope_response_time_code']
        w=self.deposition[sel]*buried*np.exp(-age/tau)*smooth(age/(.4*tau))
        orbital=np.ones(len(normals));energy_orb=0.
        if w.sum()>1e-20:
            phase=self.data['phase'][sel]+self.spin(t)-self.spin(end);direction=np.column_stack([np.cos(phase),np.sin(phase),np.zeros(len(phase))])
            kernel=np.exp(6*(normals@direction.T))@(w/w.sum())*(6/np.sinh(6))
            fraction=.70*np.exp(-max(0.,t-300)/self.config['thermal_diffusion_time_code']);orbital=1-fraction+fraction*kernel
            energy_orb=float(np.sum(w*np.gradient(self.times)[sel])*self.config['thermal_deposition_fraction'])
        sel=(self.impact_times<t)&(self.impact_times>t-4*tau)
        if not np.any(sel):return orbital
        et=self.impact_times[sel];age=t-et;w=self.impact_energy[sel]*np.exp(-age/tau)*smooth(age/(.4*tau))
        if w.sum()<1e-20:return orbital
        # Bin nearby impact directions before full-kernel summation to bound
        # the work; selected modes are a source reconstruction, not new physics.
        d=self.impact_directions[sel].copy();angle=self.spin(t)-self.spin(et);c=np.cos(angle);s=np.sin(angle);xx=d[:,0].copy();yy=d[:,1].copy();d[:,0]=c*xx-s*yy;d[:,1]=s*xx+c*yy
        # At most 32 weighted representative directions, preserving total weight.
        order=np.argsort(np.arctan2(d[:,1],d[:,0]));groups=np.array_split(order,min(32,len(order)))
        directions=[];weight=[]
        for group in groups:
            wg=w[group].sum();vec=np.sum(d[group]*w[group,None],axis=0);norm=np.linalg.norm(vec)
            if wg>0 and norm>1e-12:directions.append(vec/norm);weight.append(wg)
        if not weight:return orbital
        weight=np.array(weight);direction=np.array(directions);kernel=np.exp(8*(normals@direction.T))@(weight/weight.sum())*(8/np.sinh(8))
        mixture=min(.8,float(weight.sum()/(weight.sum()+energy_orb+1e-30)))
        return (1-mixture)*orbital+mixture*kernel
    def stellar_state(self,t,lum):
        a,dis,alive,R=self.scene.state(t);donor=self.report['donor'];base=donor['Teff']*(donor['radius']/R)**.5
        heat=float(np.interp(t,self.times,self.data['heat_code']))*G*MSUN**2/RSUN/1e7;cv=1.5*1.380649e-23/(.62*1.67262192369e-27)
        engulf=float(smooth((1.05*R-a)/(.65*R)));response=base+.06*heat/(cv*(donor['mass']-.179547731)*(MSUN/1000))*engulf
        power=min(.5*lum*engulf,4*np.pi*(R*RSUN/100)**2*SIGMA*max(0.,response**4-base**4));Teff=(base**4+power/(4*np.pi*(R*RSUN/100)**2*SIGMA))**.25
        return base,Teff,engulf,power
    def stars_native(self,t,phase,right,up,eye,L,center,lum):
        a,dis,alive,R=self.scene.state(t);xd,xc=self.scene.centers(t,phase);base,Tmean,engulf,Lstar=self.stellar_state(t,lum)
        rgb=np.zeros((self.H,self.W,3),np.float32);depth=np.full((self.H,self.W),-np.inf,np.float32);mask=np.zeros((self.H,self.W),np.float32)
        spin=float(np.interp(t,self.times,self.data['envelope_omega']));tensor=np.array([[np.interp(t,self.times,self.tidal_response[:,i,j]) for j in range(3)] for i in range(3)])
        flatten=min(.15,.35*spin**2*R**3/(self.scene.M+dis*self.scene.m))*engulf
        eigen,basis=np.linalg.eigh(np.eye(3)+tensor+np.diag([flatten/2,flatten/2,-flatten]))
        pixel=2*L/self.W;theta=float(self.spin(t));c=np.cos(theta);s=np.sin(theta);rotation=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
        e=np.array([np.cos(phase),np.sin(phase),0.])
        for k,(pos,rad,Teff) in enumerate([(xd,R,Tmean),(xc,self.scene.r2*max(1-dis,0)**(1/3),self.report['companion']['Teff'])]):
            if rad<1e-7:continue
            axes=rad*eigen/np.prod(eigen)**(1/3) if k==0 else np.full(3,rad);B=basis if k==0 else np.eye(3)
            cx=self.W/2+((pos-center)@right)/pixel;cy=self.H*.64-((pos-center)@up)/pixel;extent=max(axes)/pixel+2
            x0=max(0,int(cx-extent));x1=min(self.W,int(cx+extent+1));y0=max(0,int(cy-extent));y1=min(self.H,int(cy+extent+1))
            if x0>=x1 or y0>=y1:continue
            yy,xx=np.mgrid[y0:y1,x0:x1]
            def trace(x,y):
                plane=center+((x+.5-self.W/2)*pixel)[...,None]*right+((self.H*.64-y-.5)*pixel)[...,None]*up
                q=(plane-pos)@B/axes;direction=eye@B/axes;aa=direction@direction;bb=2*np.einsum('...i,i->...',q,direction);cc=np.sum(q*q,axis=-1)-1;disc=bb*bb-4*aa*cc
                front=(-bb+np.sqrt(np.maximum(0,disc)))/(2*aa);return q,direction,disc,front
            q,direction,disc,front=trace(xx,yy);grad=np.gradient(disc);edge=abs(disc)<.85*np.maximum(1e-30,np.hypot(grad[0],grad[1]))
            regular=(disc>=0)&~edge;iy,ix=np.nonzero(regular);px=[xx[regular].astype(float)];py=[yy[regular].astype(float)];factor=[np.ones(regular.sum())]
            ex,ey=xx[edge],yy[edge]
            for dy in [-.375,-.125,.125,.375]:
                for dx in [-.375,-.125,.125,.375]:px.append(ex+dx);py.append(ey+dy);factor.append(np.full(len(ex),1/16))
            px=np.concatenate(px);py=np.concatenate(py);weight=np.concatenate(factor);q,direction,disc,front=trace(px,py);valid=disc>=0
            px=px[valid];py=py[valid];weight=weight[valid];front=front[valid];local=q[valid]+front[:,None]*direction
            raw_world=(local*axes)@B.T;rr=np.linalg.norm(raw_world,axis=1);normal=(local/axes)@B.T;normal/=np.linalg.norm(normal,axis=1)[:,None];mu=np.clip(normal@eye,0,1)
            temp=np.full(len(mu),Teff);boost=np.ones(len(mu))
            if k==0:
                tidal=(self.scene.m/self.scene.M)*(R/max(a,R))**3
                cos_t=(raw_world@e)/np.maximum(rr,1e-12);g=np.clip((R/rr)**2-2*tidal*(1.5*cos_t*cos_t-.5),.55,1.5)
                pattern=self.thermal_pattern(t,normal) if engulf else np.ones(len(mu));temp=(base**4+max(0.,Tmean**4-base**4)*pattern)**.25*g**(.08*(1-dis))
                strength=min(.8,spin**2*R**3/(self.scene.M+dis*self.scene.m));beta=.08+.17*dis;zz,wz=np.polynomial.legendre.leggauss(32);norm=(np.sum(wz*(1-strength*(1-zz*zz))**(4*beta))/2)**.25
                temp*=np.maximum(.15,1-strength*(1-(raw_world[:,2]/rr)**2))**beta/norm
                material=(raw_world/rr[:,None])@rotation;observer=eye@rotation
                # Dimensionless source scale is mapped to a 15-Hp donor
                # granule and 6x spatial display magnification, not solar size.
                hp=self.surface_report['stars']['donor']['Hp_Mm'];source_feature_Mm=1.3
                mapping=(15*hp/source_feature_Mm)*6
                scale=R*RSUN/1e8/mapping
                source_t=self.mt[0]+t*self.report['physical_time_unit_seconds']/(self.meta['convection_clock_slowdown']*self.meta['donor_turnover_s'])*300
                samples=surface_sample(material.astype('f8'),observer.astype('f8'),mu.astype('f4'),self.morph,self.mt,source_t,scale,self.mmu,self.meso,self.active,self.active_time,self.active_time[0]+source_t-self.mt[0],self.spot_center,self.spot_east,self.spot_north)
                residual=(1-dis)+dis*.60*np.exp(-max(0.,t-300)/self.config['thermal_diffusion_time_code'])
                logI=.18*residual*samples[:,0]+.035*residual*samples[:,1];logI=np.clip(logI,-.6,.6)
                physical_factor=np.exp(logI-.5*(.18*residual)**2);temp*=physical_factor**.25
                # Illustrative activity fades as the hot surface loses its
                # assumed convective magnetic-spot signature.
                spot_weight=(1-.8*dis)*(1-.65*smooth((Tmean-6000)/3000))
                spot=1+spot_weight*(samples[:,2]-1);temp*=spot
                boost=physical_factor**(self.meta['granulation_display_contrast_exponent']-1)*spot**3.2
            I=self.atmosphere(temp,mu)*boost[:,None]
            if k==0:I=grade_photosphere(I,float(color(self.report['donor']['Teff'],self.lut)@YROW),self.grade)
            # Accumulate exact edge supersamples, preserving the companion's
            # true projected size. Other star depth is respected per sample.
            nx=np.floor(px+.5).astype(int);ny=np.floor(py+.5).astype(int);good=(nx>=0)&(nx<self.W)&(ny>=0)&(ny<self.H)
            nx=nx[good];ny=ny[good];front=front[good];weight=weight[good];I=I[good]
            good=front>depth[ny,nx];nx=nx[good];ny=ny[good];front=front[good];weight=weight[good];I=I[good]
            if k==1:
                # All samples of a secondary pixel replace the primary only
                # where the resolved secondary surface is actually in front.
                cover=np.zeros((self.H,self.W),np.float32);np.add.at(cover,(ny,nx),weight)
                rgb*=1-np.clip(cover,0,1)[...,None]
            np.add.at(rgb,(ny,nx),I*weight[:,None]);np.add.at(mask,(ny,nx),weight)
            np.maximum.at(depth,(ny,nx),front)
        return rgb,depth,np.clip(mask,0,1),lum-Lstar
    def gas_volume(self,t,x,v,T,width,w,Lgas,right,up,eye,L,center):
        selected=w>1e-7;ids=np.flatnonzero(selected)
        if not len(ids):return None
        emitting,meta=self.gas_sources(t,x,v,width,w,Lgas);p=x[ids];vel=v[ids];m=w[ids]*float(self.data['parcel_mass_msun'])
        # Reuse existing source diffusion, thermal/tidal widths and neighbor
        # data, then replace the old ribbon raster by compact volume packets.
        _,sample_mass,sample_T,_=self.material_samples(p,vel,m,emitting,ids)
        sourceT=np.sum((sample_mass*sample_T).reshape(len(ids),200),axis=1)/m
        tangent=vel/np.maximum(np.linalg.norm(vel,axis=1,keepdims=True),1e-12);radial=p/np.maximum(np.linalg.norm(p,axis=1,keepdims=True),1e-12);vertical=np.cross(radial,tangent);vertical/=np.maximum(np.linalg.norm(vertical,axis=1,keepdims=True),1e-12);across=np.cross(vertical,tangent)
        near=self.neighbor_indices;dv=vel[near]-vel[:,None,:];age=np.maximum(0,t-self.births[ids]);h=self.stream_width
        v1=np.sqrt(np.mean(np.einsum('nki,ni->nk',dv,across)**2,axis=1));v2=np.sqrt(np.mean(np.einsum('nki,ni->nk',dv,vertical)**2,axis=1))
        s1=np.sqrt(h*h+np.minimum(width[ids],.3*age*v1)**2);s2=np.sqrt(h*h+np.minimum(width[ids],.3*age*v2)**2)
        # Longitudinal material spacing from nearby launch coordinates.
        gap=np.zeros(len(ids))
        for channel,maxgap in [(0,5),(1,2),(2,1)]:
            local=np.flatnonzero(self.channel[ids]==channel)
            if len(local)<2:continue
            lo=local[:-1];hi=local[1:];good=(ids[hi]-ids[lo]<=maxgap)&(abs(self.birth_phase[ids[hi]]-self.birth_phase[ids[lo]])<.8)
            lo=lo[good];hi=hi[good];dist=np.linalg.norm(p[hi]-p[lo],axis=1);gap[lo]=np.maximum(gap[lo],dist);gap[hi]=np.maximum(gap[hi],dist)
        long=np.sqrt(h*h+gap*gap/12)
        sigma=np.stack([long,s1,s2],axis=-1)
        basis=np.stack([tangent,across,vertical],axis=-1);camera=np.column_stack([right,up,eye]);basis=np.einsum('ij,njk->nik',camera.T,basis)
        proj=(p-center)@camera;gw=self.gas_width;gh=round(gw*16/9);pixel=2*L/gw;pad=4
        core=1.8*self.scene.state(t)[3]+2.;Z=max(core*1.1,float(np.max(abs(proj[:,2])+self.support*np.max(sigma,axis=1))))
        outer=24;inner=80;negative=-np.geomspace(Z,core,outer+1);middle=np.linspace(-core,core,inner+1);positive=np.geomspace(core,Z,outer+1)
        edges=np.r_[negative[:-1],middle,positive[1:]];zs=np.diff(edges);zc=(edges[:-1]+edges[1:])/2
        shape=(len(zc),gh+2*pad,gw+2*pad);origin=np.array([-L-pad*pixel,(-.36*gh-pad)*pixel,edges[0]]);spacing=np.array([pixel,pixel,1.])
        # Resolve undersampled kernels rather than letting them miss z slabs.
        broaden=np.interp(proj[:,2],zc,zs)/np.sqrt(12)
        cov=np.einsum('nij,nkj->nik',basis*sigma[:,None,:],basis*sigma[:,None,:]);cov[:,0,0]+=pixel*pixel/12;cov[:,1,1]+=pixel*pixel/12;cov[:,2,2]+=broaden*broaden
        material_basis=basis.copy()
        values,basis=np.linalg.eigh(cov);sigma=np.sqrt(np.maximum(values,1e-8))
        sound=np.sqrt((5/3)*1.380649e-23*np.maximum(700,emitting)/(1.3*1.67262192369e-27))/((RSUN/100)/self.report['physical_time_unit_seconds'])
        dispersion=np.sqrt(np.mean(np.sum(dv*dv,axis=-1),axis=1));mach=np.minimum(4.,dispersion/np.maximum(sound,1e-8))
        forcing=1/3+.47*self.material_compression
        clumping=np.sqrt(np.log1p((forcing*mach)**2))*.90
        phase=self.birth_phase[ids]*2.7+self.channel[ids]*1.73
        mass,heat,in_view=deposit_packets(proj.astype('f8'),basis.astype('f8'),sigma.astype('f8'),m.astype('f8'),sourceT.astype('f8'),clumping.astype('f8'),phase.astype('f8'),shape,origin,spacing,zc,zs,self.support,material_basis.astype('f8'),self.clump_field)
        kappa=self.config['grey_opacity_cm2_g']*.1*(MSUN/1000)/(pixel*RSUN/100)**2
        upper=bolometric_upper(mass,heat,kappa,pixel*RSUN/100);dilution=min(1.,Lgas/max(upper,1e-30))
        metrics=dict(source_mass_msun=float(m.sum()),mass_in_volume_msun=float(in_view),gas_source_upper_W=float(upper*dilution),assigned_gas_power_W=float(Lgas),clump_log_density_sigma_max=float(clumping.max()),volume_shape=list(shape),volume_z_nonuniform=True)
        return mass,heat,origin,spacing,zc,zs,kappa,dilution,metrics
    def render_linear(self,playback):
        f=np.clip(playback/32,0,1);t=float(self.clock(f));phase=float(self.phase(t));x,v,T,h,w,lum=self.at(t);R=self.scene.state(t)[3]
        L=6.85+.30*(R-self.scene.R);inc=.82-.12*smooth((t-215)/175);az=-np.pi/2+.35*f
        eye=np.array([np.sin(inc)*np.cos(az),np.sin(inc)*np.sin(az),np.cos(inc)]);right=np.array([-np.sin(az),np.cos(az),0.]);up=np.cross(eye,right);center=np.zeros(3)
        star,depth,mask,Lgas=self.stars_native(t,phase,right,up,eye,L,center,lum);volume=self.gas_volume(t,x,v,T,h,w,Lgas,right,up,eye,L,center)
        survival=1.;metrics=dict(source_mass_msun=0.,mass_in_volume_msun=0.,gas_source_upper_W=0.,assigned_gas_power_W=Lgas)
        if volume is not None:
            mass,heat,origin,spacing,zc,zs,kappa,dilution,metrics=volume
            result,surviving=raycast(mass,heat,star,depth,origin,spacing,2*L/self.W,L,.64*self.H,zc,zs,kappa,np.log(self.lut[0]),self.lut[1],self.report['donor']['Teff'],dilution)
            survival=float(surviving.sum(dtype=np.float64)/max(star.sum(dtype=np.float64),1e-30))
        else:result=star
        _,dis,alive,_=self.scene.state(t);gas=float(w.sum()*self.data['parcel_mass_msun'])
        self.records.append(dict(playback_s=float(playback),t_code=t,remnant_mass_msun=self.scene.M+dis*self.scene.m-gas,companion_mass_msun=alive,gas_mass_msun=gas,stellar_display_flux_survival=survival,remnant_diameter_fraction=R/L,**metrics))
        return result,mask
    def frame(self,playback,fps=30,temporal_samples=2):
        # Two native renders are averaged in linear light over half a frame;
        # no frame interpolation or image upscaling is used.
        image=np.zeros((self.H,self.W,3),np.float32);mask=np.zeros((self.H,self.W),np.float32)
        for q in range(temporal_samples):
            offset=((q+.5)/temporal_samples-.5)*.5/fps
            rgb,m=self.render_linear(float(np.clip(playback+offset,0,32)));image+=rgb/temporal_samples;mask+=m/temporal_samples
        mapped=self.optics.process_masked(image,mask)
        mapped*=1-self.title_outline[...,None]*.85;mapped=mapped*(1-self.title_core[...,None])+self.title_core[...,None]*.91
        srgb=np.where(mapped<=.0031308,12.92*mapped,1.055*np.maximum(mapped,0)**(1/2.4)-.055)
        return np.uint16(np.clip(srgb,0,1)*65535+.5)
