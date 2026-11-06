"""Accepted radius/camera/nebula pipeline with a new statistical giant surface."""
import numpy as np
from .production_publication import PublicationRenderer
from .giant_convection import GiantConvection
from .evolution import PC,RSUN,smooth
from .reuse_no003 import circle_rectangle_area

class GiantCellRenderer(PublicationRenderer):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.convection=GiantConvection()
    def stellar_image(self,t,star,basis,pixel_radius,source_only=False):
        R=star['R']*PC/RSUN;weight=float(smooth(np.log(max(R,1)/1.6)/np.log(3./1.6)))
        compact=float(smooth((star['T']-9000)/14000))
        if weight==0 or pixel_radius<.65 or compact>.999:
            return super().stellar_image(t,star,basis,pixel_radius,source_only)
        # No coordinate scaling by R, T, M, pulsation or camera. Cell coordinates
        # are persistent world/star-frame directions from the independent model.
        x=(self.xx+.5-self.w/2)/pixel_radius;y=(self.h*.575-self.yy-.5)/pixel_radius;q=x*x+y*y;mask=q<1
        mu=np.sqrt(np.maximum(0,1-q[mask]));roll=np.deg2rad(-20)
        nx=np.cos(roll)*x[mask]+np.sin(roll)*y[mask];ny=-np.sin(roll)*x[mask]+np.cos(roll)*y[mask]
        normals=np.stack([nx,ny,mu],axis=-1)@basis
        f=self.convection.sample(normals,t);_,meta,_=self.convection.field(t)
        base=self.compact.surface(mu,star);temperature=star['T']*f**.25
        chrom=self.track.rgb(temperature)/np.maximum(self.track.rgb(star['T']),1e-30)
        values=base*chrom
        if compact>0:values=(1-compact)*values+compact*base
        image=np.zeros((self.h,self.w,3));image[mask]=values
        # Original circle/pixel coverage and nearest-interior edge source; no
        # changed radius, limb profile, optical bloom, detector gain or exposure.
        edge=np.abs(np.sqrt(q)-1)*pixel_radius<.8;ey,ex=np.nonzero(edge)
        for iy,ix in zip(ey,ex):
            area=float(circle_rectangle_area(ix-self.w/2,ix+1-self.w/2,self.h*.575-iy-1,self.h*.575-iy,pixel_radius))
            if area<=0:continue
            xx=(ix+.5-self.w/2)/pixel_radius;yy=(self.h*.575-iy-.5)/pixel_radius;rad=max(np.hypot(xx,yy),1e-10)
            xx*=min(1,(1-1e-8)/rad);yy*=min(1,(1-1e-8)/rad)
            iy0=int(np.clip(self.h*.575-yy*pixel_radius,0,self.h-1));ix0=int(np.clip(self.w/2+xx*pixel_radius,0,self.w-1))
            source=image[iy0,ix0] if image[iy0,ix0].sum()>0 else self.compact.surface(np.array([np.sqrt(max(0,1-xx*xx-yy*yy))]),star)[0]
            image[iy,ix]=source*area
        if weight<1:
            original,old_info=super().stellar_image(t,star,basis,pixel_radius,source_only)
            # This is a gradual source-physics handoff during departure from the
            # main sequence, in linear radiance; never an encoded image dissolve.
            image=(1-weight)*original+weight*image
        info=dict(resolved=True,source='Statistical spherical giant cells; original spectral/limb/circle-area rendering',
            statistical_surface_weight=weight,compact_weight=compact,convection=meta,
            bolometric_identity='T_local^4=T_eff^4*f; exact full-sphere integral of sampled f equals1',
            scope='Representative stochastic surface, not radiation-hydrodynamic convection or tracked real gas')
        return image,info
    def production_frame(self,index):
        image,ledger=super().production_frame(index)
        ledger['surface_model']='src.giant_convection.GiantConvection'
        ledger['source_clock_rule']='No004 main sequence uses original20:1clock; giant cells have staggered individual representative turnover clocks; neither is retimed by radius/camera/pulsation'
        return image,ledger
