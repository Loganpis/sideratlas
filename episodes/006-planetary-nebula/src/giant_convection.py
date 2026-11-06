"""Persistent, statistical spherical convection cells, not a fluid solver.

Cell centres and identities live in stellar coordinates. Radius never enters
their coordinate transform. Stellar parameters constrain *newborn* cells only.
Three populations represent granules, weaker long-lived larger structures,
and shorter-lived fine surface cells within the larger pattern.
No image atlas, cube face, spherical harmonic, frame noise or global zoom is used.
"""
from pathlib import Path
import json
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
from .evolution import PC,RSUN,smooth
from .evolution_v3 import Track
from .production_temporal_fix import PlaybackCoordinate,ResolvedEvolutionClock
ROOT=Path(__file__).resolve().parents[1]

def adopted_track():
    t=Track();p=json.loads((ROOT/'production/TEMPORAL-SAMPLING.json').read_text());p['windows']=[w for w in p['windows'] if w['end_s']<8]
    t.clock=ResolvedEvolutionClock(t.clock,PlaybackCoordinate(p));return t

def constraints(star):
    R=star['R']*PC/RSUN;T=star['T'];M=star['M'];g=27400*M/R**2;lg=np.log10(g)
    valid=(3600<=T<=6750 and 1<=lg<=4.5 and T>=600*lg)
    # Tremblay+2013 eqs15/16, [Fe/H]=0. Outside the giant grid use
    # an explicitly flagged H_p / sound-speed continuation from its boundary.
    if T>=3600 and lg>=1:
        term=max(100,T-300*lg)**1.75;length_km=13.5/g*term;decay_s=2.08/g*term
        method='Tremblay2013 fit' if valid else 'fit outside convective-grid domain'
    else:
        anchor=3300.**1.75
        length_km=13.5/10*anchor*(T/3600)*(10/g)
        decay_s=2.08/10*anchor*np.sqrt(T/3600)*(10/g)
        method='flagged H_p/c_s continuation from3600K/logg1; not validated giant RHD'
    theta=length_km/(R*695700.)
    # Resolved display sampling floor; explicit spatial enhancement at small R.
    shown=float(np.clip(theta,.040,.28))
    tau_display=float(.70+.20*np.log10(1+decay_s/1000.))
    return dict(R_Rsun=R,T_K=T,M_Msun=M,logg_cgs=lg,fit_in_grid=bool(valid),scale_method=method,
        reference_length_km=length_km,reference_decay_s=decay_s,reference_angular_diameter=theta,
        birth_angular_diameter=shown,spatial_display_factor=shown/theta,
        birth_decay_display_s=tau_display,reference_drift_km_s=.23*length_km/decay_s)

def normalize(v):return v/np.linalg.norm(v,axis=-1,keepdims=True)

@njit(cache=True)
def paint_cells(centres,east,north,size,activity,ellipticity,xyz,nlat,nlon):
    first=np.full((nlat,nlon),1e8);second=first.copy();owner=np.full((nlat,nlon),-1,np.int32)
    dth=np.pi/(nlat-1);dph=2*np.pi/nlon
    for k in range(len(centres)):
        if activity[k]<1e-6:continue
        c=centres[k];theta=np.arccos(min(1.,max(-1.,c[2])));phi=np.arctan2(c[1],c[0])%(2*np.pi)
        support=min(1.25,2.2*size[k]);coscut=np.cos(support)
        lo=max(0,int(np.floor((theta-support)/dth)));hi=min(nlat-1,int(np.ceil((theta+support)/dth)))
        span=np.pi if theta<support or theta>np.pi-support else np.arcsin(min(1.,np.sin(support)/np.sin(theta)))
        mid=int(np.floor(phi/dph));count=int(np.ceil(span/dph))+1
        left=mid-count;right=mid+count
        if 2*count+1>=nlon:left=0;right=nlon-1
        penalty=-2*np.log(activity[k]);sigma=.50*size[k]
        for iy in range(lo,hi+1):
            for jj in range(left,right+1):
                ix=jj%nlon;p=xyz[iy,ix];dot=p[0]*c[0]+p[1]*c[1]+p[2]*c[2]
                if dot<coscut:continue
                a=p[0]*east[k,0]+p[1]*east[k,1]+p[2]*east[k,2]
                b=p[0]*north[k,0]+p[1]*north[k,1]+p[2]*north[k,2]
                shape=1+ellipticity[k]*(a*a-b*b)/(a*a+b*b+1e-20)
                score=2*max(0.,1-dot)/sigma**2*shape+penalty
                if score<first[iy,ix]:
                    second[iy,ix]=first[iy,ix];first[iy,ix]=score;owner[iy,ix]=k
                elif score<second[iy,ix]:second[iy,ix]=score
    heat=np.zeros_like(first)
    for iy in range(nlat):
        for ix in range(nlon):
            if first[iy,ix]>1e7:heat[iy,ix]=.45;continue
            gap=min(20.,second[iy,ix]-first[iy,ix])
            # Warm interiors, narrow dark boundaries; no prescribed polygon grid.
            lane=1-np.exp(-(gap/.90)**2)
            dome=np.exp(-.085*max(0.,first[iy,ix]))
            heat[iy,ix]=.30+.70*lane*dome
    return heat,owner,first

@njit(cache=True)
def sample_sphere(field,normals):
    nlat,nlon=field.shape;out=np.empty(len(normals))
    for k in range(len(normals)):
        n=normals[k];theta=np.arccos(min(1.,max(-1.,n[2])));phi=np.arctan2(n[1],n[0])%(2*np.pi)
        y=theta/np.pi*(nlat-1);x=phi/(2*np.pi)*nlon;iy=min(nlat-2,int(y));ix=int(x)%nlon
        fy=y-iy;fx=x-np.floor(x);jx=(ix+1)%nlon
        out[k]=(1-fy)*((1-fx)*field[iy,ix]+fx*field[iy,jx])+fy*((1-fx)*field[iy+1,ix]+fx*field[iy+1,jx])
    return out

class GiantConvection:
    def __init__(self,nlat=769,nlon=1536,seed=6006,reference_star=None):
        self.track=adopted_track();self.nlat=nlat;self.nlon=nlon;self.seed=seed;self.reference_star=reference_star
        th=np.linspace(0,np.pi,nlat);ph=np.arange(nlon)*2*np.pi/nlon
        self.xyz=np.stack(np.broadcast_arrays(np.sin(th[:,None])*np.cos(ph),np.sin(th[:,None])*np.sin(ph),np.cos(th[:,None])+np.zeros_like(ph)),axis=-1)
        h=np.pi/(nlat-1);self.row_weight=np.sin(th)*2*(1-np.cos(h))/h;self.row_weight[[0,-1]]=1-np.sin(h)/h
        assert abs(self.row_weight.sum()-2)<1e-10
        self.populations=[self._population(family) for family in [0,1,2]]
        self.cache={};self.last_meta={}
    def mean(self,a):return float(np.sum(a*self.row_weight[:,None])/(2*self.nlon))
    def _population(self,family):
        rng=np.random.default_rng(self.seed+family*100003);rows=[]
        # Continuous random birth times are drawn once. There is no random call
        # during field evaluation and no periodic sequence or per-frame reset.
        for lo in np.arange(-8.,16.,.08):
            t=lo+.04;s=self.track.state(float(self.track.clock(max(0,t)))) if self.reference_star is None else self.reference_star;c=constraints(s)
            size=c['birth_angular_diameter'];tau=c['birth_decay_display_s']
            if family==1:size=min(.65,max(.14,3.2*size));tau*=3.3
            elif family==2:size=max(.024,.50*size);tau*=.60
            lifetime=2.6*tau;number=16/size**2
            count=rng.poisson(number/lifetime*.08)
            if count==0:continue
            births=lo+rng.uniform(0,.08,count);n=normalize(rng.normal(size=(count,3)))
            e=rng.normal(size=(count,3));e=normalize(e-n*np.sum(e*n,axis=1)[:,None])
            sizes=size*np.exp(np.clip(rng.normal(0,.16,count),-.36,.36));lives=lifetime*np.exp(np.clip(rng.normal(0,.22,count),-.45,.45))
            phase=rng.uniform(0,2*np.pi,count);ellip=rng.uniform(.10,.28,count)
            for j in range(count):rows.append(np.r_[births[j],lives[j],tau,sizes[j],n[j],e[j],phase[j],ellip[j]])
        return np.asarray(rows)
    def cells(self,t,family):
        p=self.populations[family];age=t-p[:,0];keep=(age>0)&(age<p[:,1]);ids=np.flatnonzero(keep);a=p[keep];age=age[keep];x=age/a[:,2];u=age/a[:,1]
        n=a[:,4:7];e=a[:,7:10];north=np.cross(n,e);theta=a[:,3]
        # Individual tangent drift follows cell size/turnover, not R(t), camera
        # or the radial pulsation. Direction bends smoothly within a cell life.
        c=normalize(n+theta[:,None]*(.23*x[:,None]*e+.065*(np.sin(1.3*x+a[:,10])-np.sin(a[:,10]))[:,None]*north))
        east=normalize(e-c*np.sum(c*e,axis=1)[:,None]);nn=np.cross(c,east)
        rotation=.30*np.sin(.9*x+a[:,10]);east2=np.cos(rotation)[:,None]*east+np.sin(rotation)[:,None]*nn;north2=np.cross(c,east2)
        activity=np.sin(np.pi*u)**2
        ellipticity=a[:,11]*(.8+.2*np.sin(1.7*x+a[:,10]))
        return c,east2,north2,theta,activity,ellipticity,ids
    def field(self,t):
        # Exact requested time; cache only repeated requests from disk/global
        # normalization. No shutter-dependent or camera-dependent evolution.
        key=round(float(t),10)
        if key in self.cache:return self.cache[key]
        maps=[];meta={};owners=[]
        for family in [0,1,2]:
            c,e,n,s,w,ellip,ids=self.cells(t,family)
            heat,owner,score=paint_cells(c,e,n,s,w,ellip,self.xyz,self.nlat,self.nlon)
            heat[0]=self.mean_row(heat[0]);heat[-1]=self.mean_row(heat[-1])
            mean=self.mean(heat);std=np.sqrt(self.mean((heat-mean)**2));z=(heat-mean)/max(std,1e-9)
            maps.append(z);owners.append(np.where(owner>=0,ids[np.maximum(0,owner)],-1))
            meta[f'population{family}']=dict(active_cells=len(c),median_diameter_rad=float(np.median(s)),median_lifetime_display_s=float(np.median(self.populations[family][ids,1])),uncovered_area_fraction=self.mean((score>1e7).astype(float)))
        # Approximately18–22% bolometric granulation contrast, with a weaker
        # large-cell component. Normalization conserves the assigned luminosity.
        star=self.track.state(float(self.track.clock(max(0,t)))) if self.reference_star is None else self.reference_star;c=constraints(star)
        amplitude=.18+.020*float(smooth((2.5-c['logg_cgs'])/1.5))
        bol=np.exp(amplitude*maps[0]+.045*maps[1]+.105*maps[2]);bol/=self.mean(bol)
        meta.update(constraints=c,bolometric_surface_mean=self.mean(bol),bolometric_surface_RMS=float(np.sqrt(self.mean((bol-1)**2))),
            physical_claim='Literature-informed statistical cells; no calculated giant fluid flow or shock field',source='Persistent spherical sites; newborn sizes only depend on stellar parameters',representative_atmospheric_clock=True)
        result=(bol,meta,owners)
        self.cache[key]=result
        while len(self.cache)>3:self.cache.pop(next(iter(self.cache)))
        self.last_meta=meta;return result
    @staticmethod
    def mean_row(a):return float(np.mean(a))
    def sample(self,normals,t):return sample_sphere(self.field(t)[0],np.asarray(normals,dtype=float))
