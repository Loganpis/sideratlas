"""Stable evaluation of the unchanged early-wind cumulative mass function.

Before the native AGB source table, only differences of cohort values are used.
Changing their additive zero and evaluating the same PCHIP in shifted form avoids
subtracting ~0.46Msun cumulative values to recover a ~5e-9Msun solar wind.
No mass rescaling, track fit, physical parameter or cached trajectory is changed.
"""
import numpy as np

def stable_early_state(model,original_state,t):
    if t>=model.birth_edges[0]:return original_state(t)
    previous=model.cohort;poly=model.track.massfn
    def cohort(birth,rate=False):
        if rate:return previous(birth,rate=True)
        x=np.asarray(birth,dtype=float)
        assert x.ndim>=1 and x.shape[0]==model.mu.size and np.all(x<model.birth_edges[0])
        j=np.clip(np.searchsorted(poly.x,x,side='right')-1,0,len(poly.x)-2)
        h=x-poly.x[j]
        # M_initial - M(x), evaluated without forming M(x) near one solar mass.
        c=poly.c
        loss=(c[3,0]-c[3,j])-h*(c[2,j]+h*(c[1,j]+h*c[0,j]))
        factor=model.A.ravel().reshape((model.mu.size,)+(1,)*(x.ndim-1))*model.omega/(4*np.pi)
        return loss*factor
    model.cohort=cohort
    try:state=original_state(t)
    finally:model.cohort=previous
    state['report']['early_wind_mass_evaluation']='same PCHIP/source angular allocation, cumulative zero shifted to initial MESA mass; stable polynomial evaluation; no renormalization'
    return state
