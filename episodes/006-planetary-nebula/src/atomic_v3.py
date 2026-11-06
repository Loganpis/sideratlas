"""Actual PyNeb line coefficients + reduced multigroup direct-ray ionization.

Not Cloudy: Te=10000K is prescribed; metals do not compete for photon absorption,
He uses approximate Case-B demand, diffuse photons are local. O/N recombination
uses published historical RR+low-temperature DR fits, not invented constants.
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from .evolution import KB,PC,MH,MSUN
ROOT=Path(__file__).resolve().parents[1]
EV=1.602176634e-12

def cross_section(z,electrons,energy):
    tab=np.loadtxt(ROOT/'data/atomic/verner-1996-photo.dat')
    row=tab[(tab[:,0]==z)&(tab[:,1]==electrons)][0]
    eth,emax,E0,sig0,ya,P,yw,y0,y1=row[2:]
    E=np.asarray(energy);x=E/E0-y0;y=np.sqrt(x*x+y1*y1)
    sig=sig0*((x-1)**2+yw*yw)*y**(.5*P-5.5)*(1+np.sqrt(y/ya))**(-P)*1e-18
    return np.where((E>=eth)&(E<=emax),sig,0.)

def recombination(z,final_electrons,T=10000.):
    # The table labels the recombined/final ion's number of electrons.
    pl=np.loadtxt(ROOT/'data/atomic/pl.dat');rr=pl[(pl[:,0]==z)&(pl[:,1]==final_electrons)][0]
    t=T/1e4;value=rr[2]*t**(-rr[3])
    dr=np.loadtxt(ROOT/'data/atomic/ldr.dat');sel=dr[(dr[:,0]==z)&(dr[:,1]==final_electrons)]
    if len(sel):
        a,b,c,d,f=sel[0 if T<20000 else -1,2:]
        value+=max(0,1e-12*(a/t+b+c*t+d*t*t)*t**(-1.5)*np.exp(-f/t))
    return float(value)

class Radiation:
    def __init__(self,track):
        self.track=track
        a=np.load(ROOT/'data/atomic/pyneb-v3.npz')
        self.T=a['T'];self.ne=a['ne'];self.coeff={k:RegularGridInterpolator((np.log(self.T),np.log(self.ne)),np.log(a[k]),bounds_error=True) for k in ['ha','hb','o3','n2']}
        # Boundaries include each relevant ionization threshold. No oxygen line
        # is assigned to a region merely because its density is high.
        self.edges=np.array([13.6,24.59,29.6,35.12,47.45,54.42,54.94,77.41,110.,200.,400.])
        g,w=np.polynomial.legendre.leggauss(12)
        self.E=.5*(self.edges[:-1,None]+self.edges[1:,None])+.5*np.diff(self.edges)[:,None]*g
        self.w=.5*np.diff(self.edges)[:,None]*w
        self.sigma=np.array([cross_section(z,n,self.E) for z,n in [(1,1),(2,2),(2,1),(8,7),(8,6),(7,6),(7,5)]])
        self.rec=np.array([recombination(z,n) for z,n in [(8,7),(8,6),(7,6),(7,5)]])
    def spectrum(self,star):
        u=self.E*EV/(KB*star['T']);phot=self.E**2*np.exp(-u)/(-np.expm1(-u))
        norm=15*star['L']/np.pi**4/(KB*star['T']/EV)**4/EV
        counts=norm*np.sum(phot*self.w,axis=-1)
        mean=np.sum(phot*self.w*self.E,axis=-1)/np.maximum(np.sum(phot*self.w,axis=-1),1e-200)
        sig=np.sum(self.sigma*phot[None,:,:]*self.w[None,:,:],axis=-1)/np.maximum(np.sum(phot*self.w,axis=-1)[None,:],1e-200)
        return counts,mean,sig
    def coefficients(self,ne,T=10000.):
        n=np.clip(ne,self.ne[0],self.ne[-1]);temp=np.broadcast_to(T,n.shape)
        points=np.c_[np.log(temp.ravel()),np.log(n.ravel())]
        return {k:np.exp(f(points)).reshape(n.shape) for k,f in self.coeff.items()}
    def solve(self,star,r,vol,nH,rho,kap,dr,omega):
        shape=nH.shape;nr=shape[-1];n=int(np.prod(shape[:-1]));h=nH.reshape(n,nr)
        rad=np.broadcast_to(r,shape).reshape(n,nr);volume=np.broadcast_to(vol,shape).reshape(n,nr)
        opacity=np.broadcast_to(kap,shape).reshape(n,nr);width=np.broadcast_to(dr,shape).reshape(n,nr)
        q,energy,sigma=self.spectrum(star);budget=np.broadcast_to(q*omega/(4*np.pi),(n,len(q))).copy()
        x=np.zeros_like(h);ne=np.zeros_like(h);of=np.zeros_like(h);nf=np.zeros_like(h)
        qdust=0.;qH=0.;qHe=0.;energyH=0.;heat=0.;alpha=2.59e-13
        rhoi=rho.reshape(n,nr)
        for j in range(nr):
            before=budget.copy();budget*=np.exp(-3*opacity[:,j,None]*rhoi[:,j,None]*width[:,j,None]*PC)
            qdust+=float(np.sum(before-budget))
            incident=budget.copy();area=omega*(rad[:,j]*PC)**2
            # H and He share the incident groups according to their cross sections;
            # cap each species' absorption by its cell recombination demand.
            hweight=sigma[0]/np.maximum(sigma[0]+.1*sigma[1],1e-99)
            avail=(budget*hweight[None,:]).sum(axis=1)
            need=alpha*1.1*h[:,j]**2*volume[:,j]
            used=np.minimum(avail,need)
            allocated=budget*hweight[None,:]*np.divide(used,avail,out=np.zeros_like(used),where=avail>0)[:,None]
            budget-=allocated
            x[:,j]=np.sqrt(np.divide(used,need,out=np.zeros_like(used),where=need>0))
            ne[:,j]=1.1*x[:,j]*h[:,j]
            qH+=float(used.sum());energyH+=float(np.sum(allocated*energy[None,:]))*EV
            heat+=float(np.sum(allocated*np.maximum(0,energy[None,:]-13.6)))*EV
            # Approximate He I/II absorption; not an independently solved full He atom.
            heweight=np.where(sigma[1]>0,1.,0.)
            avhe=(budget*heweight[None,:]).sum(axis=1)
            needhe=4.3e-13*ne[:,j]*(.1*x[:,j]*h[:,j])*volume[:,j]
            uhe=np.minimum(avhe,needhe)
            budget-=budget*heweight[None,:]*np.divide(uhe,avhe,out=np.zeros_like(uhe),where=avhe>0)[:,None]
            qHe+=float(uhe.sum())
            mid=.5*(incident+budget)/np.maximum(area[:,None],1e-99)
            rates=mid@sigma[3:].T
            ratios=rates/(np.maximum(ne[:,j,None],1e-30)*self.rec[None,:])
            o1,o2,n1,n2=ratios.T
            of[:,j]=x[:,j]*o1/(1+o1+o1*o2)
            nf[:,j]=x[:,j]/(1+n1+n1*n2)
        ne=ne.reshape(shape);x=x.reshape(shape);of=of.reshape(shape);nf=nf.reshape(shape)
        c=self.coefficients(ne)
        line=np.stack([c['ha']*ne*nH*x,c['hb']*ne*nH*x,c['o3']*ne*nH*self.track.O_H*of,c['n2']*ne*nH*self.track.N_H*nf],axis=-1)
        metal=float(np.sum(line[...,2:].sum(axis=-1)*vol));mcap=min(1,heat/max(metal,1e-99));line[...,2:]*=mcap
        lum=np.sum(line*np.broadcast_to(vol,shape)[...,None],axis=(0,1,2));ecap=min(1,energyH/max(lum.sum(),1e-99));line*=ecap
        Q=float(q.sum())
        report=dict(Q_H_s=Q,Q_dust_s=qdust,Q_recomb_s=qH,Q_He_absorbed_s=qHe,Q_escape_s=float(budget.sum()),photon_relative_closure=abs(Q-qdust-qH-qHe-budget.sum())/max(Q,1),
                    absorbed_ionizing_luminosity_erg_s=energyH,photoelectron_heat_erg_s=heat,metal_heating_cap=mcap,line_energy_cap=ecap,
                    electron_temperature_K=10000.,temperature_scope='prescribed photoionized thermostat, not solved thermal balance',
                    line_luminosities_erg_s=dict(zip(['Halpha6563','Hbeta4861','OIII5007','NII6584'],[float(i*ecap) for i in lum])),
                    ne_cm3_p05_p50_p95=np.quantile(ne[x>.1],[.05,.5,.95]).tolist() if np.any(x>.1) else [0,0,0],
                    O2_fraction_range=[float(of.min()),float(of.max())],N1_fraction_range=[float(nf.min()),float(nf.max())],H1_fraction_range=[float(x.min()),float(x.max())])
        return line,x,ne,report,of,nf

    def solve_timed(self,star,r,vol,nH,rho,kap,dr,omega,filling):
        """Line radiation on a separately conserved, time-dependent H front.

        Filling is ionized VOLUME fraction of a sharp-front material cell. The
        ionized part has ne=1.1nH; recombination is proportional to filling, not
        to filling squared. Ionic rates use threshold-resolved blackbody groups,
        but H transport is grey. Recombination can release stored ionization
        energy when the front retreats; it is not capped by current UV alone.
        """
        shape=nH.shape;nr=shape[-1];n=int(np.prod(shape[:-1]))
        h=nH.reshape(n,nr);f=filling.reshape(n,nr);rad=r.reshape(n,nr);volume=vol.reshape(n,nr)
        tau=np.sum((kap*rho*dr*PC).reshape(n,nr),axis=-1)
        q,energy,sigma=self.spectrum(star);Q=float(q.sum());trans=np.exp(-np.minimum(700,3*tau))
        incident=np.broadcast_to(q*omega/(4*np.pi),(n,len(q))).copy()
        qdust=float(np.sum(incident*(1-trans[:,None])))
        high=q*(self.edges[:-1]>=24.59);he_fraction=min(.12,float(high.sum())/max(Q,1e-99))
        budget=incident*trans[:,None]
        he_remove=trans[:,None]*omega/(4*np.pi)*high[None,:]*he_fraction*Q/max(high.sum(),1e-99)
        qhe=float(he_remove.sum());budget-=he_remove
        budget=np.maximum(budget,0);of=np.zeros_like(h);nf=np.zeros_like(h)
        rec=2.59e-13*1.1*h*h*f*volume
        # He photons are a separate travelling packet, not a screen placed ahead
        # of every emitting cell. Distribute their absorption through the HII
        # volume in proportion to its recombination demand. This approximate
        # closure retains hard photons in the inner nebula and avoids erasing
        # OIII/NIII everywhere when the total He allowance uses those photons.
        he_packet=he_remove.copy();he_norm=np.maximum(rec.sum(axis=1),1e-99)
        absorbed=0.;labs=0.;heat=0.
        for j in range(nr):
            before=budget.copy();avail=budget.sum(axis=1)
            # The first partially ionized material cell is the ionization front.
            # Its photons ionize new H as well as replace recombinations.
            want=np.where(f[:,j]<1-1e-9,avail,rec[:,j])
            use=np.minimum(want,avail)
            allocate=budget*np.divide(use,avail,out=np.zeros_like(use),where=avail>0)[:,None]
            budget-=allocate;absorbed+=float(use.sum())
            labs+=float(np.sum(allocate*energy[None,:]))*EV
            heat+=float(np.sum(allocate*np.maximum(0,energy[None,:]-13.6)))*EV
            area=omega*(rad[:,j]*PC)**2
            he_here=he_remove*np.divide(rec[:,j],he_norm)[:,None]
            travelling_before=before+he_packet
            he_packet=np.maximum(0,he_packet-he_here)
            travelling_after=budget+he_packet
            rates=(.5*(travelling_before+travelling_after)/np.maximum(area[:,None],1e-99))@sigma[3:].T
            ratios=rates/(np.maximum(1.1*h[:,j,None],1e-30)*self.rec[None,:])
            o1,o2,n1,n2=ratios.T
            of[:,j]=f[:,j]*o1/(1+o1+o1*o2)
            nf[:,j]=f[:,j]/(1+n1+n1*n2)
        ni=1.1*nH;F=filling
        c=self.coefficients(ni)
        line=np.stack([c['ha']*ni*nH*F,c['hb']*ni*nH*F,c['o3']*ni*nH*self.track.O_H*of.reshape(shape),c['n2']*ni*nH*self.track.N_H*nf.reshape(shape)],axis=-1)
        metal=float(np.sum(line[...,2:].sum(axis=-1)*vol));mcap=min(1,heat/max(metal,1e-99));line[...,2:]*=mcap
        qr=float(rec.sum());storage=absorbed-qr
        available_power=labs+max(0,-storage)*13.6*EV
        lum=np.sum(line*vol[...,None],axis=(0,1,2));ecap=min(1,available_power/max(lum.sum(),1e-99));line*=ecap
        esc=float(budget.sum())
        report=dict(Q_H_s=Q,Q_dust_s=qdust,Q_He_absorbed_s=qhe,Q_recomb_s=qr,Q_net_ionizations_s=storage,Q_escape_s=esc,
                    photon_relative_closure=abs(Q-qdust-qhe-qr-storage-esc)/max(Q,1),
                    instantaneous_absorbed_H_ionizing_luminosity_erg_s=labs,absorbed_ionizing_luminosity_erg_s=available_power,
                    recombination_storage_release_erg_s=max(0,-storage)*13.6*EV,photoelectron_heat_erg_s=heat,metal_heating_cap=mcap,line_energy_cap=ecap,
                    electron_temperature_K=10000.,temperature_scope='prescribed ionized thermostat,100K neutral; no solved thermal balance',
                    ionization_scope='time-dependent sharp H front; grey dusty screen/12percent He allowance; threshold-resolved metal rates, historical RR+DR, no charge exchange',
                    line_luminosities_erg_s=dict(zip(['Halpha6563','Hbeta4861','OIII5007','NII6584'],[float(i*ecap) for i in lum])),
                    O2_fraction_range=[float(of.min()),float(of.max())],N1_fraction_range=[float(nf.min()),float(nf.max())],H1_fraction_range=[float(F.min()),float(F.max())])
        return line,F,ni,report,of.reshape(shape),nf.reshape(shape)
