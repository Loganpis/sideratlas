"""Selected trailing22 Cox–Gomez kernel and exact gas-grid force basis."""
import numpy as np
from numba import njit
from . import CONFIG
from .dynamics import BASE, KM_S, VREF2, field as axis_field, envelope
EMPTY=np.empty((0,9))
WAVES2=np.array(CONFIG['model']['spiral_parameters'])

@njit(cache=True)
def spatial(R,z,w):
    """F(R)*sech(Kz/B)^B and its exact logarithmic R/z derivatives."""
    k=np.sqrt(w[0]*w[0]+w[1]*w[1])/R; kp=-k/R
    q=k*w[7]; qp=-q/R
    B=q*(1+.4*q); Bp=qp*(1+.8*q)
    den=1+.3*q; num=1+q+.3*q*q
    D=num/den; Dp=((1+.6*q)*den-.3*num)/den**2*qp
    k0=np.sqrt(w[0]*w[0]+w[1]*w[1])/4; q0=k0*w[7]
    D0=(1+q0+.3*q0*q0)/(1+.3*q0)
    u0=(R/w[9])**4; inner=-np.expm1(-u0)
    f=np.exp(-(R-4)/w[8]-(R/w[10])**6)*k0*D0/(k*D)*inner
    li=4*R**3/w[9]**4*np.exp(-u0)/max(inner,1e-90)
    dlogf=-1/w[8]-6*R**5/w[10]**6+1/R-Dp/D+li
    u=k*z/B; au=abs(u)
    logcosh=au+np.log1p(np.exp(-2*au))-np.log(2.)
    V=np.exp(-B*logcosh); th=np.tanh(u)
    dlogv=-Bp*logcosh-B*th*u*(-1/R-Bp/B)
    dlogz=-k*th
    return f*V,dlogf+dlogv,dlogz

@njit(cache=True)
def field2(x,y,z,t,base=BASE,waves=WAVES2):
    ax,ay,az,phi,unused1,unused2=axis_field(x,y,z,t,base,EMPTY)
    R=max(np.hypot(x,y),1e-8); power=0.; torque=0.
    if R>.02:
        theta=np.arctan2(y,x)
        for j in range(len(waves)):
            w=waves[j]; A,Adot=envelope(t,w)
            if A==0 and Adot==0:continue
            F,Fr,Fz=spatial(R,z,w)
            phase=w[0]*(theta-w[2]*t)-w[1]*np.log(R/4)+w[3]
            c=np.cos(phase); s=np.sin(phase)
            ar=A*F*(Fr*c+w[1]/R*s); at=-A*F*w[0]/R*s
            ax+=ar*x/R-at*y/R; ay+=ar*y/R+at*x/R
            az+=A*F*Fz*c; phi-=A*F*c
            tau=-A*F*w[0]*s
            torque+=tau; power+=w[2]*tau-Adot*F*c
    return ax,ay,az,phi,power,torque

class Gravity2:
    """Precompute the same scalar-potential modes on the gas midplane."""
    def __init__(self,X,Y,waves=WAVES2):
        from .dynamics import evaluate
        self.X=X;self.Y=Y;self.waves=waves
        p=np.c_[X.ravel(),Y.ravel(),X.ravel()*0]
        a=evaluate(p,0.,BASE,EMPTY).reshape(*X.shape,6)
        self.base=np.moveaxis(a[:,:,[0,1,3]],-1,0)
        R=np.maximum(np.hypot(X,Y),1e-8); theta=np.arctan2(Y,X)
        self.basis=[]
        for w in waves:
            raw=np.array([spatial(r,0.,w) for r in R.ravel()]).reshape(*R.shape,3)
            F=raw[:,:,0]; dlog=raw[:,:,1]
            ph=w[0]*theta-w[1]*np.log(R/4)+w[3]; C=np.cos(ph);S=np.sin(ph)
            pair=[]
            for c,s in ((C,S),(S,-C)):
                ar=F*(dlog*c+w[1]/R*s); at=-F*w[0]/R*s
                pair.append(np.array([ar*X/R-at*Y/R,ar*Y/R+at*X/R,-F*c,-F*w[0]*s]))
            self.basis.append(pair)

    def at(self,t,forcing=True):
        out=np.zeros((5,*self.X.shape));out[:3]=self.base
        if forcing:
            for w,(b0,b1) in zip(self.waves,self.basis):
                A,Adot=envelope(t,w);phase=w[0]*w[2]*t
                b=b0*np.cos(phase)+b1*np.sin(phase)
                out[:3]+=A*b[:3];out[3]+=Adot*b[2]+w[2]*A*b[3];out[4]+=A*b[3]
        return out
