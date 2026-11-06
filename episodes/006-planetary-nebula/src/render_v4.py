"""Imaging/transition pass on the unchanged accepted v3 physical snapshots.

The tested continual-seeding closure is NOT the default: its visual difference
was negligible. Actual No004/003 stellar sources and the v3 ray projector remain.
Changes here are explicitly camera, detector HDR and monotonic display response.
"""
import numpy as np
import hashlib
from .render_v3 import Renderer as Original,project_material
from .evolution import PC,RSUN,smooth
from .reuse_no003 import CompactOptics

Y=np.array([.2126729,.7151522,.0721750])

class Renderer(Original):
    def __init__(self,w=288,h=512,depth=112,dynamics='radial-v3t-12x24-144.npz'):
        super().__init__(w,h,depth,dynamics)
        self.sharp_optics=CompactOptics((h,w),100000.,scatter_fraction=.01)
        scale=w/1080
        self.sharp_optics.cfg.update(airy_first_zero_px_at_550nm_1080=max(2.2,.95/scale),
                                    jitter_sigma_px_at_1080=max(.4,.3/scale))
        self.sharp_optics.pixel_integrated_psf(100000.)

    def geometry_v3(self,t,star,state):
        if t<=8:return super().geometry_v3(t,star,state)
        giant=float(smooth(np.log10(max(star['R']*PC/RSUN,1.))/2.2))
        near=star['R']/(.50+.27*giant)
        power=state['escaped_ir_power'].sum(axis=(0,1));cum=np.cumsum(power)
        ir95=float(np.interp(.95*cum[-1],cum,state['r'])) if cum[-1]>0 else star['R']
        env=max(near,ir95/.65)
        thick=float(smooth((state['report']['dust_IR_tau_median']-.3)/.5))
        # The IR channel reveals an extended wind. Frame that physical volume
        # as the channel comes up, rather than remaining inside it for six s.
        w=max(float(smooth((t-8)/2)),thick)
        half=np.exp((1-w)*np.log(near)+w*np.log(env))
        if t>=13.5:
            if 'line_fields' in state:
                p=state['line_fields'].sum(axis=-1).sum(axis=(0,1))*np.diff(state['edges']**3)
                c=np.cumsum(np.maximum(p,0));rr=float(np.interp(.90*c[-1],c,state['r'])) if c[-1]>0 else 0.
                target=max(.003,rr/.78)
            else:target=.003
            w=float(smooth((t-13.5)/6))
            half=np.exp((1-w)*np.log(half)+w*np.log(target))
        # Continuous fade at8s preserves the original opening's actual camera.
        w=float(smooth((t-8)/1.5));old=super().geometry_v3(t,star,state)
        half=float(np.exp((1-w)*np.log(old)+w*np.log(half)))
        # Encoded review revealed a camera-induced apparent size reversal as
        # the first bright shell appeared. Balance that view only, with C1
        # log-scale offsets; physical material positions/velocities stay fixed.
        if 14<t<22:
            knots=[14.,15.5,16.5,18.,20.,22.];factors=[1.,1.6,2.3,.8,.8,1.]
            k=int(np.clip(np.searchsorted(knots,t)-1,0,len(knots)-2))
            z=float(smooth((t-knots[k])/(knots[k+1]-knots[k])))
            half*=np.exp((1-z)*np.log(factors[k])+z*np.log(factors[k+1]))
        return half

    def frame(self,t,fps=10,mode='full'):
        original=super().frame(t,fps,mode)
        record=self.records[-1];record['physical_model']='accepted v3; experimental response not selected'
        record['physical_snapshot_sha256']=hashlib.sha256(self.last_layers['state']['fields'].tobytes()).hexdigest()
        if t<=8:return original
        state=self.last_layers['state'];star=state['star'].copy();basis,_=self.basis(t)
        half=self.last_layers['halfwidth'];pixel=2*half/self.w
        pulse=record['AGB_radius_pulsation_factor'];star['R']*=pulse;star['T']/=np.sqrt(pulse)
        rpx=record['stellar_radius_pixels'];hdr=self.last_layers['hdr'].copy()
        center=self.compact.surface(np.array([1.]),star)[0]@Y
        fade=float(smooth((t-13)/5));tr=record['post_AGB_reference_yr']
        nebula=3e12*(max(100.,tr)/16000)**3
        exposure=float(np.exp((1-fade)*np.log(.95/max(center,1e-20))+fade*np.log(nebula)))
        exposure=max(exposure,record['intrinsic_point_reference_exposure'])
        old_ir=record['IR_display_gain']
        ir_window=float(smooth((t-8.8)/1.2))*(1-float(smooth((t-14.8)/2.8)))
        power=state['escaped_ir_power'].sum(axis=(0,1));cumulative=np.cumsum(power)
        ir95=float(np.interp(.95*cumulative[-1],cumulative,state['r'])) if cumulative[-1]>0 else star['R']
        mean_IR=float(cumulative[-1])/(4*np.pi**2*(max(ir95,star['R'])*PC)**2*self.track.Yref)
        # A compact, bright dusty photosphere must not receive the same long
        # IR exposure as the later dilute envelope. This changes detector gain,
        # not the escaped/reprocessed luminosity or its opacity.
        normalized_IR=4000./(1+4000.*mean_IR/.18)
        ir_exposure=(1-ir_window)*.12*exposure+ir_window*normalized_IR
        p99=0.
        if t<18:
            fields=state['fields'].copy();fields[...,:3]=0
            unresolved=state['r']<pixel;fields[...,unresolved,4]=0
            ir=project_material(fields,np.log(state['r'][0]),np.log(state['r'][1]/state['r'][0]),self.w,self.h,half,
                                star['R'] if rpx>=.65 else 0.,np.zeros_like(hdr),basis,self.radial_nodes(state),1.)
            irpoint=float(state['escaped_ir_power'][...,unresolved].sum())/(4*np.pi*(pixel*PC)**2*self.track.Yref)
            col=np.array([1.,.42,.16]);col/=col@Y
            ir[round(self.h*.575),self.w//2]+=irpoint*col
            p99=float(np.quantile(ir@Y,.99))
            normalized_IR=4000./(1+4000.*(mean_IR/.18+p99/.25))
            ir_exposure=(1-ir_window)*.12*exposure+ir_window*normalized_IR
            visible=hdr-old_ir*ir
            assert visible.min()>-1e-7*max(hdr.max(),1e-99)
            image=np.maximum(visible,0)*exposure+ir*ir_exposure
        else:image=hdr*exposure
        line_exposure=exposure
        if 15<t<18 and 'line_fields' in state:
            f=state['fields'].copy();f[...,:3]=100*state['line_fields'];f[...,4:]=0
            lines=project_material(f,np.log(state['r'][0]),np.log(state['r'][1]/state['r'][0]),self.w,self.h,half,
                                   star['R'] if rpx>=.65 else 0.,np.zeros_like(hdr),basis,self.radial_nodes(state),0.)
            w=float(smooth((t-15)/2))
            line_exposure=float(np.exp((1-w)*np.log(exposure)+w*np.log(nebula)))
            image+=lines*(line_exposure-exposure)
        point_exposure=exposure
        if rpx<.65:
            mu,w=np.polynomial.legendre.leggauss(24);mu=(mu+1)/2;w/=2
            flux=2*np.pi*rpx*rpx*np.sum(self.compact.surface(mu,star)*(mu*w)[:,None],axis=0)
            tau=record['central_components']['stellar_radial_tau_V'];atten=flux*np.exp(-min(700,tau))
            gate=float(smooth((t-14)/4))*float(smooth((star['T']-6500)/3500))
            reference=8/max(flux@Y,1e-30)
            point_exposure=float(np.exp((1-gate)*np.log(exposure)+gate*np.log(reference)))
            image[round(self.h*.575),self.w//2]+=atten*(point_exposure-exposure)
            assert image.min()>-1e-7*max(image.max(),1e-99)
            image=np.maximum(image,0)
        linear=self.sharp_optics.linear_process(image)
        exponent=1+.35*float(smooth((t-14)/5))
        maximum=np.max(linear,axis=-1)
        mapped=np.divide(maximum**exponent,1+maximum**exponent)
        display=linear*np.divide(mapped,maximum,out=np.zeros_like(mapped),where=maximum>0)[...,None]
        blend=float(smooth((t-8)/1))
        old=original.astype(float)/65535
        old=np.where(old<=.04045,old/12.92,((old+.055)/1.055)**2.4)
        display=display*(1-self.alpha[...,None])+self.alpha[...,None]*.9
        display=old*(1-blend)+display*blend
        srgb=np.where(display<=.0031308,12.92*display,1.055*np.maximum(display,0)**(1/2.4)-.055)
        record.update(v4_visible_exposure=exposure,v4_IR_detector_exposure=normalized_IR,v4_IR_exposure_ceiling=4000.,v4_IR_reference_mean_Y=.18,v4_IR_window=ir_window,
                      v4_IR_reference_p99_Y=.25,v4_IR_projected_p99_Y=p99,v4_IR95_radius_pc=ir95,
                      v4_line_exposure=line_exposure,v4_effective_line_to_continuum_gain=100*line_exposure/exposure,
                      v4_point_exposure=point_exposure,v4_point_intrinsic_reference_Y=8.,v4_tone_exponent=exponent,
                      v4_camera='IR95 wind view8–10s, line90 handoff13.5–19.5s; C1 transition scale balance14–22s; physical geometry unchanged',
                      v4_imaging='virtual positive PSF2.2px first zero at1080/.95px preview floor;1percent wing; independent short stellar-continuum HDR exposure; no changed luminosity')
        self.last_layers['v4_detector_input']=image
        return np.uint16(np.clip(srgb,0,1)*65535+.5)

    def evidence(self):
        r=super().evidence();r.update(physical_model='unchanged selected v3 radial/atomic/photon model',
                response_experiment='not selected: insufficient visual benefit',v4_detector_sampling=self.sharp_optics.sampling_check,
                v4_display='fixed line colors/base100x sensitivity; independent15–18s line exposure; IR reference meanY.18/p99Y.25/ceiling4000 window8.8–17.6s; visible exposure13–18s; attenuated point8Y reference; late contrast exponent1.35')
        return r
