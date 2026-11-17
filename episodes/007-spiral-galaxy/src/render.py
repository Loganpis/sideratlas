"""Final native linear-light/shutter pipeline; selected implementation only."""
import argparse,json,time
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.ndimage import gaussian_filter,gaussian_filter1d
from PIL import Image,ImageDraw,ImageFont
from . import ROOT,CONFIG
from .history import History
from .camera import Camera3
from .populations import FrozenPop,YROW,LSUN,smooth
from .transfer import dust_column,cic_rgb,compact_sources,bulge_column,columns_at_points,slab_transfer_parallel as slab_transfer,EXT,LINE_EXT,GAS_COLUMN_TO_G_CM2
from .optics import white_balance_matrix,radial_tone as optical_tone
from .video_export import encoder,metadata
CACHE=ROOT/'data/optical.npz'
FPS=30

class Renderer:
    def __init__(self,run,width=1080,nz=64):
        self.settings=dict(neighbors=64,line_gain=4.,tone='radial-balance')
        self.history=History(run);self.w=width;self.h=2*round(width*8/9);self.metrics=[];self.pop=FrozenPop()
        self.camera=Camera3('A');self.center_y=.53;self.nz=nz;self.partial=False;self.source_interval=1.
        self.view_start=560.;self.time_rate=5.;self.line_gain=self.settings['line_gain'];self.old_quadrature_fraction=.01
        self.members=4
        samples=np.random.default_rng(17007).normal(0,1,(len(self.history.birth),8,5))
        self.vel_samples=samples[:,:self.members,:2];self.pos_samples=samples[:,:self.members,2:]
        a=np.load(self.history.run/'cohorts.npz');self.cohort=a['id'];self.mass=a['mass']
        with np.load(CACHE) as a:self.rv=a['rv'].copy();self.omega=a['omega'].copy();self.kappa=a['kappa'].copy()
        self.ng=int(self.history.config['gas_grid']);self.dx=36/self.ng;self.source_ng=896;self.source_half=self.source_ng*self.dx/2
        self.gridlevels=np.geomspace(.10,.85,8);self.last_layers={};self.last_world={}
        yy,xx=np.mgrid[:self.h,:self.w];self.xx=xx+.5-self.w*.5;self.yy=self.h*self.center_y-yy-.5
        self.font_title=ImageFont.truetype(str(ROOT/'fonts/lmroman12-bold.otf'),round(self.w*.064))
        self.font_small=ImageFont.truetype(str(ROOT/'fonts/lmroman10-regular.otf'),round(self.w*.025))
        layer=Image.new('RGBA',(self.w,self.h));d=ImageDraw.Draw(layer)
        d.text((self.w*.075,self.h*.130),'@sideratlas',font=self.font_small,fill=(235,233,227),stroke_width=1,stroke_fill=(0,0,0),anchor='la')
        d.text((self.w*.925,self.h*.130),'NO.007',font=self.font_small,fill=(235,233,227),stroke_width=1,stroke_fill=(0,0,0),anchor='ra')
        d.text((self.w*.5,self.h*.157),'THE SPIRAL GALAXY',font=self.font_title,fill=(242,240,233),stroke_width=1,stroke_fill=(0,0,0),anchor='ma')
        ar=np.asarray(layer).astype(float)/255;self.text_alpha=ar[:,:,3:];self.text_rgb=ar[:,:,:3]
        self.wb=white_balance_matrix(self.pop.solar_rgb)
        self.frames=[]
        if len(self.mass)!=100000 or self.ng!=512:raise ValueError("Selected production requires100kstars/512gas cells")

    @lru_cache(maxsize=4)
    def cached_sources(self,t):
        points,unused,sigma=self.history.sample(t)
        old=self.source_grid(points[self.cohort==0],self.mass[self.cohort==0],0)
        mid=self.source_grid(points[self.cohort==1],self.mass[self.cohort==1],1)
        return old,mid

    def mass_fields(self,t,points=None,direct=False):
        if direct:
            if points is None:points=self.history.sample(t)[0]
            return (self.source_grid(points[self.cohort==0],self.mass[self.cohort==0],0),
                    self.source_grid(points[self.cohort==1],self.mass[self.cohort==1],1))
        a=self.history.start+np.floor((t-self.history.start)/self.source_interval)*self.source_interval
        b=min(a+self.source_interval,float(self.history.config['run_end_Myr']))
        if b<=a:return self.cached_sources(a)
        f=(t-a)/(b-a);o0,m0=self.cached_sources(a)
        if f<1e-10:return o0,m0
        o1,m1=self.cached_sources(b)
        return (1-f)*o0+f*o1,(1-f)*m0+f*m1

    def vertical_table(self,points):
        tab=[];profiles=[]
        for k in (0,1):
            p=points[self.cohort==k];r=np.hypot(p[:,0],p[:,1]);b=np.minimum((r/.5).astype(int),39)
            counts=np.bincount(b,minlength=40);z2=np.bincount(b,weights=p[:,2]**2,minlength=40)
            good=counts>20;h=np.sqrt(np.divide(z2,counts,out=np.zeros(40),where=counts>0))
            h=np.interp(np.arange(40),np.flatnonzero(good),h[good]);h=gaussian_filter1d(h,1.)
            hr=np.interp(np.arange(261)*.1,(np.arange(40)+.5)*.5,h)
            hr=np.clip(hr,.04,.8);z=np.linspace(-1.2,1.2,self.nz,endpoint=False)+1.2/self.nz
            tab.append(np.exp(-.5*(z[:,None]/hr[None,:])**2)/(np.sqrt(2*np.pi)*hr[None,:]));profiles.append(hr)
        return np.array(tab),profiles

    def source_grid(self,points,mass,k):
        radius=np.hypot(points[:,0],points[:,1]);rd=2.7 if k==0 else 3.8
        den=len(points)/(2*np.pi*rd**2)*np.exp(-radius/rd)
        hs=np.clip(np.sqrt(64/(2*np.pi*np.maximum(den,1e-10))),.10,.85)
        u=np.interp(np.log(hs),np.log(self.gridlevels),np.arange(len(self.gridlevels)))
        xy=(points[:,:2]+self.source_half)/self.dx
        result=np.zeros((self.source_ng,self.source_ng))
        for j,h in enumerate(self.gridlevels):
            wt=np.maximum(0,1-abs(u-j));mask=wt>0
            light=np.repeat((mass[mask]*wt[mask])[:,None],3,axis=1)
            comp=cic_rgb(xy[mask],light,self.source_ng,self.source_ng)[:,:,0]
            result+=gaussian_filter(comp,h/self.dx,mode='constant')
        return result/self.dx**2

    def project(self,points,er,eu,scale):
        return np.c_[points@er*scale+self.w*.5,self.h*self.center_y-points@eu*scale]

    def old_quadrature(self,points,c0,c1,sigma,n,er,eu,scale):
        # The same calculated orbital samples redistribute a small fraction
        # of existing light. No new masses, positions or velocities are made.
        intrinsic=self.mass[:,None]*np.where(self.cohort[:,None]==0,c0,c1)
        intrinsic*=self.old_quadrature_fraction
        tau=dust_column(points,sigma,self.dx,n)*220*GAS_COLUMN_TO_G_CM2
        received=intrinsic*np.exp(-tau[:,None]*EXT)
        img=cic_rgb(self.project(points,er,eu,scale),received,self.h,self.w)*scale**2
        img=gaussian_filter(img,(.45,.45,0),mode='constant')
        return img,float(intrinsic.sum(axis=0)@YROW),float(received.sum(axis=0)@YROW)

    def source(self,s):
        t=self.view_start+self.time_rate*s;points,young,sigma=self.history.sample(t)
        view,er,eu,n=self.camera.at(s);scale=self.w/view;area=1/scale**2
        old,mid=self.mass_fields(t,points)
        c0,ml0,age0=self.pop.cohort(0,t);c1,ml1,age1=self.pop.cohort(1,t)
        vertical,heights=self.vertical_table(points)
        age=t-self.history.birth[:len(young)];mass=self.history.mass[:len(young)]
        light,Q,Lbol=self.pop.young(age,mass);blend=smooth((age-80)/40)
        aged_xy=(young[:,:2]+self.source_half)/self.dx
        aged_img=cic_rgb(aged_xy,light*blend[:,None],self.source_ng,self.source_ng)/self.dx**2
        aged=np.moveaxis(gaussian_filter(aged_img,(.18/self.dx,.18/self.dx,0),mode='constant'),-1,0).copy()
        dxp=self.xx/scale;dyp=self.yy/scale
        plane=dxp[:,:,None]*er+dyp[:,:,None]*eu
        x0=plane[:,:,0]-plane[:,:,2]*n[0]/n[2];y0=plane[:,:,1]-plane[:,:,2]*n[1]/n[2]
        B=bulge_column(np.maximum(np.hypot(dxp,dyp),.3/scale))
        fields,tau=slab_transfer(old,mid,aged,sigma,x0,y0,n,self.dx,vertical,c0,c1,B,self.pop.bulge_rgb,1.2)
        fine,quad_intrinsic,quad_received=self.old_quadrature(points,c0,c1,sigma,n,er,eu,scale)
        f=self.old_quadrature_fraction
        old_continuum=(1-f)*(fields[0]+fields[1])+fine+fields[3]
        active=blend<1.;ids=np.flatnonzero(active);p=young[active];ages=age[active]
        rr=np.maximum(np.hypot(p[:,0],p[:,1]),.2);om=np.interp(rr,self.rv,self.omega)[:,None];kap=np.interp(rr,self.rv,self.kappa)[:,None]
        aa=ages[:,None];u=self.vel_samples[ids,:,0]*(2*.001022712165045695);v=self.vel_samples[ids,:,1]*(2*.001022712165045695)
        dr=u/kap*np.sin(kap*aa)+2*om*v/kap**2*(1-np.cos(kap*aa))+.005*self.pos_samples[ids,:,0]
        dp=2*om*u/kap**2*(np.cos(kap*aa)-1)+v*((1-4*om**2/kap**2)*aa+4*om**2/kap**3*np.sin(kap*aa))+.005*self.pos_samples[ids,:,1]
        er0=np.c_[p[:,0]/rr,p[:,1]/rr,np.zeros(len(p))];et0=np.c_[-p[:,1]/rr,p[:,0]/rr,np.zeros(len(p))]
        member=p[:,None,:]+dr[:,:,None]*er0[:,None,:]+dp[:,:,None]*et0[:,None,:]
        member[:,:,2]+=.005*self.pos_samples[ids,:,2];member=member.reshape(-1,3)
        emitted=np.repeat((light[active]*(1-blend[active,None]))/self.members,self.members,axis=0)
        av=dust_column(member,sigma,self.dx,n)*220*GAS_COLUMN_TO_G_CM2
        received=emitted*np.exp(-av[:,None]*EXT)
        compact=cic_rgb(self.project(member,er,eu,scale),received,self.h,self.w)/area
        compact=gaussian_filter(compact,(.4,.4,0),mode='constant')
        gas_at=columns_at_points(young,sigma,self.dx)
        cell_xy=np.clip(((young[:,:2]+18)/self.dx).astype(int),0,self.ng-1);cell_ids=cell_xy[:,1]*self.ng+cell_xy[:,0]
        available_mass=sigma[cell_xy[:,1],cell_xy[:,0]]*self.dx**2*.3
        available_mass=np.where(np.max(abs(young[:,:2]),axis=1)<18,available_mass,0.)
        line_rgb,line_lum,radius,Qused=self.pop.gas_limited_lines(Q,gas_at,available_mass,cell_ids)
        line_col=dust_column(young,sigma,self.dx,n)*220*GAS_COLUMN_TO_G_CM2
        line_received=(line_lum['ha'][:,None]/LSUN*self.pop.line_rgb['ha']*np.exp(-LINE_EXT[0]*line_col[:,None])
                      +line_lum['hb'][:,None]/LSUN*self.pop.line_rgb['hb']*np.exp(-LINE_EXT[1]*line_col[:,None]))
        lines=compact_sources(self.project(young,er,eu,scale),line_received,radius/1000*scale/np.sqrt(5),self.h,self.w)/area
        lines=.9*lines+.1*gaussian_filter(lines,(.12*scale,.12*scale,0),mode='constant')
        young_image=fields[2]+compact;continuum=old_continuum+young_image
        raw=continuum+self.line_gain*lines
        linear=.975*gaussian_filter(raw,(.38,.38,0),mode='constant')+.022*gaussian_filter(raw,(2.5*self.w/540,2.5*self.w/540,0),mode='constant')+.003*gaussian_filter(raw,(7*self.w/540,7*self.w/540,0),mode='constant')
        intrinsic_Y=float((self.mass[:,None]*np.where(self.cohort[:,None]==0,c0,c1)).sum(axis=0)@YROW)
        rec=dict(playback_s=float(s),physical_Myr=float(t),composition=self.camera.composition,view_width_kpc=view,
            inclination_deg=float(np.rad2deg(np.arccos(n[2]))),roll_deg=float(np.rad2deg(self.camera.roll(s))),
            camera_azimuth_deg=float(np.rad2deg(self.camera.az(s))),center_y_fraction=self.center_y,
            intrinsic_preexisting_Y=intrinsic_Y,quad_intrinsic_Y=quad_intrinsic,quad_received_Y=quad_received,
            quad_light_fraction=quad_intrinsic/intrinsic_Y,continuum_Y_in_frame=float(np.sum(continuum@YROW)*area),
            reconstructed_old_mid_mass_fraction=float((old.sum()+mid.sum())*self.dx**2/self.mass.sum()),
            tauV_max=float(tau.max()),tauV_nonempty_p50_p90=np.quantile(tau[tau>.02],[.5,.9]).tolist(),
            intrinsic_Halpha_erg_s=float(line_lum['ha'].sum()),intrinsic_Hbeta_erg_s=float(line_lum['hb'].sum()),
            used_Q_s=float(Qused.sum()),ionizing_Q_s=float(Q.sum()),
            text_foreground_max_error=0.)
        self.metrics.append(rec)
        self.last_layers=dict(linear=linear,continuum=continuum,dust_tau=tau,old_continuum=old_continuum,young=young_image,lines=lines)
        rec.update(dynamical_model='trailing22',source_neighbors=64,line_display_gain=4.)
        return linear,self.last_layers,rec

    def frame(self, s, shutter=True, keep_layers=False):
        r = self
        offset = min(1/(8*FPS), s, 32-s) if shutter else 0.
        times = [s-offset, s+offset] if offset else [s]
        linear = None; records = []; photometry = []; saved = {}
        for t in times:
            _, layers, rec = r.source(t)
            records.append(rec.copy())
            area = (rec['view_width_kpc']/r.w)**2
            photometry.append({k: float((layers[k]@YROW).sum()*area)
                for k in ['old_continuum', 'young', 'continuum', 'lines', 'linear']})
            a = layers['linear']
            linear = a.copy()/len(times) if linear is None else linear+a/len(times)
            if keep_layers:
                for k in ['continuum', 'old_continuum', 'young', 'lines', 'dust_tau']:
                    saved[k] = (layers[k].astype('f4')/len(times) if k not in saved
                                else saved[k]+layers[k].astype('f4')/len(times))
        view, er, eu, n = r.camera.at(s)
        # Preserve even the accepted display radius's floating operation order.
        plane = r.xx[:, :, None]*(view/r.w)*er+r.yy[:, :, None]*(view/r.w)*eu
        X = plane[:, :, 0]-plane[:, :, 2]*n[0]/n[2]
        Y = plane[:, :, 1]-plane[:, :, 2]*n[1]/n[2]
        display = optical_tone(linear, view/r.w, r.wb, r.settings['tone'],
                               projected_radius=np.hypot(X, Y))
        image = display*(1-r.text_alpha)+r.text_rgb*r.text_alpha
        record = dict(frame_index=round(s*FPS), playback_s=s, physical_Myr=560+5*s,
            view_width_kpc=view, inclination_deg=float(np.rad2deg(np.arccos(n[2]))),
            roll_deg=float(np.rad2deg(r.camera.roll(s))),
            camera_azimuth_deg=float(np.rad2deg(r.camera.az(s))), center_y_fraction=.53,
            shutter_sample_s=times, shutter_sample_records=records, component_photometry=photometry,
            foreground_text_max_error=float(np.max(abs(image[r.text_alpha[:, :, 0]>.999]-r.text_rgb[r.text_alpha[:, :, 0]>.999]))),
            finite_pixels=bool(np.isfinite(image).all()), source_members=4,
            scene_clipped_fraction=float(np.mean(np.max(display, axis=2)>=1-1e-12)),
            central_clipped_fraction=float(np.mean(np.max(display, axis=2)[np.hypot(X,Y)<1.]>=1-1e-12)))
        self.frames.append(record)
        if keep_layers:
            saved.update(linear=linear.astype('f4'), pre_text=display.astype('f4'))
            return image, saved, record
        return image


def main():
    p=argparse.ArgumentParser(description='Reconstruct the selected movie from original or regenerated numerical histories')
    p.add_argument('--run',type=Path,default=ROOT/'work/calculation')
    p.add_argument('--output',type=Path,default=ROOT/'work/spiral-galaxy.mp4')
    p.add_argument('--start-frame',type=int,default=0);p.add_argument('--end-frame',type=int,default=960)
    a=p.parse_args()
    if not 0<=a.start_frame<a.end_frame<=960:raise ValueError('Absolute production indices must be within0..960')
    if a.output.exists():raise FileExistsError('Preserve existing media; choose another output')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    r=Renderer(a.run);enc=encoder(a.output,1080,1920,30)
    for i in range(a.start_frame,a.end_frame):
        image=r.frame(i/30);enc.stdin.write(np.rint(image*65535).astype('<u2').tobytes())
        if i%30==0:print(i,a.end_frame,flush=True)
    enc.stdin.close()
    if enc.wait():raise RuntimeError('Encoder failed')
    a.output.with_suffix('.json').write_text(json.dumps(dict(frames=r.frames,encoding=metadata()),indent=2)+'\n')

if __name__=='__main__':main()
