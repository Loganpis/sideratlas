"""Render physical-mass parcel columns with a declared grey thermal closure.

This preview uses actual Planck spectral color and formal layer attenuation,
but prescribed temperatures, opacity, parcel widths and limited source powers.
It is not a synthetic telescope image or resolved radiation hydrodynamics.
"""

from pathlib import Path

import argparse,json,hashlib,sys,time

import numpy as np

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

from PIL import Image,ImageDraw,ImageFont,ImageFilter

from src.video_export import encoder

from src.stellar_structure import G,MSUN,RSUN

H=6.62607015e-34;C=299792458.;K=1.380649e-23;SIGMA=5.670374419e-8

H=6.62607015e-34;C=299792458.;K=1.380649e-23;SIGMA=5.670374419e-8

H=6.62607015e-34;C=299792458.;K=1.380649e-23;SIGMA=5.670374419e-8

H=6.62607015e-34;C=299792458.;K=1.380649e-23;SIGMA=5.670374419e-8

XYZ_RGB=np.array([[3.2404542,-1.5371385,-.4985314],[-.969266,1.8760108,.041556],[.0556434,-.2040259,1.0572252]])

def smooth(x):x=np.clip(x,0,1);return x*x*(3-2*x)

def make_lut():
    cmf=np.loadtxt(ROOT/'data/spectral/CIE_xyz_1931_2deg.csv',delimiter=',');wave=cmf[:,0]*1e-9
    temp=np.geomspace(300,200000,700);B=2*H*C*C/wave[None,:]**5/np.expm1(H*C/(wave[None,:]*K*temp[:,None]))
    xyz=np.trapz(B[:,:,None]*cmf[None,:,1:4],wave,axis=1);rgb=np.maximum(0,xyz@XYZ_RGB.T)
    B0=2*H*C*C/wave**5/np.expm1(H*C/(wave*K*5018.53142376455));base=np.trapz(B0*cmf[:,2],wave)
    return temp,rgb/base

def color(T,lut):
    temps,rgb=lut;T=np.asarray(T);return np.stack([np.interp(np.log(T),np.log(temps),rgb[:,j]) for j in range(3)],axis=-1)
