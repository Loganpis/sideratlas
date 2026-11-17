"""Read the original snapshots with Hermite stellar and linear gas interpolation."""
from pathlib import Path
from functools import lru_cache
import json
import numpy as np

def hermite(x0,v0,x1,v1,u,dt):
    return ((2*u**3-3*u*u+1)*x0+(u**3-2*u*u+u)*dt*v0
            +(-2*u**3+3*u*u)*x1+(u**3-u*u)*dt*v1)

class History:
    def __init__(self,run):
        self.run=Path(run)
        r=json.loads((self.run/'run-report.json').read_text())
        if not r['passed']: raise RuntimeError('Run did not pass its declared preview numerical checks')
        review=self.run/'independent-review.json'
        if review.exists() and not json.loads(review.read_text()).get('accepted_for_rendering',False):
            raise RuntimeError('Independent scientific review rejected this run for rendering')
        self.config=json.loads((self.run/'config.json').read_text())
        self.start=float(self.config.get('record_start_Myr',200.));self.step=float(self.config.get('record_step_Myr',2.));self.gas_step=float(self.config.get('gas_record_step_Myr',.5))
        self.snap=json.loads((self.run/'snapshots.json').read_text())
        self.gas_snap=json.loads((self.run/'gas-snapshots.json').read_text())
        a=np.load(self.run/'associations.npz'); self.birth=a['birth']; self.mass=a['mass']

    @lru_cache(maxsize=3)
    def load(self,i):
        with np.load(self.run/'snapshots'/self.snap[i]['file']) as a:
            return {k:a[k].astype(float) for k in a.files}

    @lru_cache(maxsize=3)
    def load_gas(self,i):
        return np.load(self.run/'snapshots'/self.gas_snap[i]['file']).astype(float)

    def sample(self,t):
        u0=np.clip((t-self.start)/self.step,0,len(self.snap)-1)
        i=min(int(np.floor(u0)),len(self.snap)-2); u=float(u0-i)
        a=self.load(i); b=self.load(i+1)
        old=hermite(a['old_x'],a['old_v'],b['old_x'],b['old_v'],u,self.step)
        count=np.searchsorted(self.birth,t+1e-8,side='right'); common=min(len(a['young_x']),count)
        young=np.empty((count,3))
        young[:common]=hermite(a['young_x'][:common],a['young_v'][:common],b['young_x'][:common],b['young_v'][:common],u,self.step)
        if count>common:
            dt=(self.start+self.step*(i+1))-t
            young[common:]=b['young_x'][common:count]-dt*b['young_v'][common:count]
        ug=np.clip((t-self.start)/self.gas_step,0,len(self.gas_snap)-1)
        j=min(int(np.floor(ug)),len(self.gas_snap)-2); f=ug-j
        sigma=(1-f)*self.load_gas(j)+f*self.load_gas(j+1)
        return old,young,sigma
