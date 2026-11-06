#include <metal_stdlib>
using namespace metal;

inline float attenuation_integral(float tau) {
    if(tau<1e-6f)return 1.f;
    // Stable (1-exp(-tau))/tau; Metal has no expm1 intrinsic.
    if(tau<.05f)return 1.f-tau*.5f+tau*tau*(1.f/6.f-tau/24.f+tau*tau/120.f);
    return (1.f-exp(-min(700.f,tau)))/tau;
}

inline void sample_field(device const float *field, float3 pos,
                         constant float *p, constant uint *n, thread float *out) {
    for (uint c=0;c<6;c++) out[c]=0.f;
    float r=length(pos); if(r<=0.f)return;
    float fr=(log(r)-p[0])/p[1];
    if(fr<0.f || fr>float(n[4]-1))return;
    float fu=(pos.z/r+1.f)*.5f*float(n[2])-.5f;
    float phi=atan2(pos.y,pos.x);if(phi<0.f)phi+=2.f*M_PI_F;
    float fp=phi/(2.f*M_PI_F)*float(n[3])-.5f;
    int iu=int(floor(fu)),ip=int(floor(fp)),ir=min(int(fr),int(n[4])-2);
    float du=fu-float(iu),dp=fp-float(ip),dr=fr-float(ir);
    for(int a=0;a<2;a++){
        uint u=uint(clamp(iu+a,0,int(n[2])-1));float wa=a?du:1.f-du;
        for(int b=0;b<2;b++){
            uint v=uint((ip+b+int(n[3]))%int(n[3]));float wb=b?dp:1.f-dp;
            for(int c=0;c<2;c++){
                float wc=c?dr:1.f-dr;
                uint base=((u*n[3]+v)*n[4]+uint(ir+c))*6;
                for(uint k=0;k<6;k++)out[k]+=wa*wb*wc*field[base+k];
            }
        }
    }
}

kernel void project_material(device const float *field [[buffer(0)]],
                             device const float *nodes [[buffer(1)]],
                             device const float *star [[buffer(2)]],
                             constant float *p [[buffer(3)]],
                             constant uint *n [[buffer(4)]],
                             device float *out [[buffer(5)]],
                             uint2 id [[thread_position_in_grid]]) {
    if(id.x>=n[0] || id.y>=n[1])return;
    float3 total=0.f;float totalIR=0.f;
    float pixel=2.f*p[2]/float(n[0]);
    float3 bx=float3(p[8],p[9],p[10]);
    float3 by=float3(p[11],p[12],p[13]);
    float3 bz=float3(p[14],p[15],p[16]);
    uint count=n[6];
    for(uint ray=0;ray<count;ray++) {
        float ox=0.f,oy=0.f;
        if(count==4){ox=(ray%2)? .25f:-.25f;oy=(ray/2)? .25f:-.25f;}
        float y=(.575f*float(n[1])-float(id.y)-.5f+oy)*pixel;
        float x=(float(id.x)+.5f-float(n[0])*.5f+ox)*pixel;
        float b2=x*x+y*y;float3 I=0.f;float Iir=0.f,zstar=-INFINITY;
        if(b2<p[3]*p[3]){
            zstar=sqrt(p[3]*p[3]-b2);
            uint k=(id.y*n[0]+id.x)*3;I=float3(star[k],star[k+1],star[k+2]);
        }
        for(uint side=0;side<2;side++){
            for(uint j=0;j<n[5]-1;j++){
                uint k=side==0?n[5]-2-j:j;
                float lo=nodes[k],hi=nodes[k+1];if(hi*hi<=b2)continue;
                float zl=sqrt(max(0.f,lo*lo-b2)),zh=sqrt(max(0.f,hi*hi-b2));
                float z0=side==0?-zh:zl,z1=side==0?-zl:zh;
                z0=max(z0,zstar);if(z1<=z0)continue;
                for(uint q=0;q<2;q++){
                    float zz=z0+(float(q)+.5f)*.5f*(z1-z0),ds=.5f*(z1-z0);
                    float f[6];sample_field(field,x*bx+y*by+zz*bz,p,n,f);
                    float tau=max(0.f,f[3]*ds),att=exp(-min(700.f,tau));
                    float fac=ds*attenuation_integral(tau);
                    I=I*att+float3(f[0],f[1],f[2])*p[5]*fac;
                    float ti=max(0.f,f[5]*ds),ai=exp(-min(700.f,ti));
                    float fi=ds*attenuation_integral(ti);
                    Iir=Iir*ai+f[4]*p[5]*fi;
                }
            }
        }
        total+=I;totalIR+=Iir;
    }
    total/=float(count);totalIR/=float(count);uint k=(id.y*n[0]+id.x)*4;
    out[k]=total.x;out[k+1]=total.y;out[k+2]=total.z;out[k+3]=totalIR;
}
