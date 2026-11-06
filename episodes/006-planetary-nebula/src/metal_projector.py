"""Same accepted ray quadrature on Metal; physical states remain CPU results."""
from pathlib import Path
import subprocess,struct,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

class MetalProjector:
    def __init__(self,spatial_samples=1):
        self.log=(ROOT/'validation/production-metal-worker.log').open('ab')
        self.process=subprocess.Popen([str(ROOT/'native/build/projector'),str(ROOT/'native/projector.metal')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log)
        line=self.process.stdout.readline()
        if not line:raise RuntimeError('GPU worker did not initialize; inspect production-metal-worker.log')
        self.info=json.loads(line);self.spatial_samples=spatial_samples;self.calls=[]
    def read(self,n):
        chunks=[];remaining=n
        while remaining:
            b=self.process.stdout.read(remaining)
            if not b:raise RuntimeError('GPU worker stopped')
            chunks.append(b);remaining-=len(b)
        return b''.join(chunks)
    def __call__(self,field,loglo,dlog,w,h,halfwidth,starR,star_image,basis,rnodes,ir_gain):
        begin=time.monotonic();f=np.array(field,dtype='<f4',order='C',copy=True);star=np.array(star_image,dtype='<f4',order='C',copy=True);nodes=np.ascontiguousarray(rnodes,dtype='<f4')
        # GPU arithmetic flushes subnormal floats. Change emission units during
        # transport so even the faint IR channel stays normal; opacity is untouched.
        unit_scale=1e20
        f[...,:3]*=np.float32(unit_scale);f[...,4]*=np.float32(unit_scale);star*=np.float32(unit_scale)
        assert np.all(np.isfinite(f)) and np.all(np.isfinite(star))
        params=np.zeros(24,dtype='<f4');params[:6]=[loglo,dlog,halfwidth,starR,ir_gain,3.085677581491367e18]
        params[8:17]=np.asarray(basis).ravel()
        samples=self.spatial_samples if starR==0 else 1
        header=struct.pack('<12I',0x53414e36,w,h,*f.shape[:3],len(nodes),samples,f.nbytes,star.nbytes,nodes.nbytes,params.nbytes)
        self.process.stdin.write(header)
        for a in [f,nodes,star,params]:self.process.stdin.write(memoryview(a).cast('B'))
        self.process.stdin.flush();magic,size=struct.unpack('<2I',self.read(8));assert magic==0x53414e36 and size==w*h*4*4
        raw=np.frombuffer(self.read(size),dtype='<f4').reshape(h,w,4)
        image=raw[...,:3].astype(float)/unit_scale
        color=np.array([1.,.42,.16]);color/=color@np.array([.2126729,.7151522,.0721750])
        image+=float(ir_gain)*(raw[...,3,None].astype(float)/unit_scale)*color
        assert np.all(np.isfinite(image)) and image.min()>=0
        self.calls.append(dict(elapsed_s=time.monotonic()-begin,width=w,height=h,spatial_samples=samples,transport_emission_unit_scale=unit_scale))
        return image
    def close(self):
        self.process.stdin.close();code=self.process.wait();self.log.close()
        if code:raise RuntimeError(f'GPU worker exit{code}')
