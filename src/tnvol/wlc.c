/* SPDX-License-Identifier: AGPL-3.0-only
 * Copyright (c) 2026 Yackov Mozgori
 *
 * tnvol / wlc: 3D CDF 9/7 wavelet + context-modelled binary range coder.
 * This file is part of tnvol. tnvol is free software under AGPL-3.0-only;
 * see the LICENSE file. For a commercial license see COMMERCIAL.md.
 */
// wlc: 3D CDF9/7 wavelet + context-modelled binary range coder (tnvol invention, 2026-11-17)
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#if defined(__AVX2__)
#include <immintrin.h>
#define WLC_AVX2 1
#else
#define WLC_AVX2 0
#endif
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <windows.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif
#else
#include <sys/mman.h>
#endif
#ifdef _MSC_VER
#define restrict __restrict
#endif
#ifdef _WIN32
#define WLC_API __declspec(dllexport)
#else
#define WLC_API
#endif
#ifndef POPW
#define POPW 0
#endif
#ifndef HUGEP
#define HUGEP 0
#endif
/* no MADV_HUGEPAGE by default: with THP defrag=madvise a fragmented host does synchronous compaction in the page fault (one gate point burned 467 s of sys CPU) */
static void* balloc(size_t b){
  if(!b) b=1;
#ifdef _WIN32
  return malloc(b);
#else
  size_t A=1<<21; size_t r=(b+A-1)/A*A; void* p=NULL;
  if(posix_memalign(&p,A,r)) return NULL;
#ifdef MADV_HUGEPAGE
  if(p && HUGEP) madvise(p,r,MADV_HUGEPAGE);
#endif
#if defined(MADV_POPULATE_WRITE)
  if(p && POPW && r>=(1u<<22)) madvise(p,r,MADV_POPULATE_WRITE);
#else
  if(p && POPW && r>=(1u<<22)) madvise(p,r,23 /* MADV_POPULATE_WRITE */);
#endif
  return p;
#endif
}
WLC_API void* wlc_alloc(long b){ return balloc((size_t)b); }
WLC_API void wlc_free(void* p){ free(p); }
#ifdef _WIN32
static double now(void){
  static LARGE_INTEGER f; static int init=0; LARGE_INTEGER c;
  if(!init){ QueryPerformanceFrequency(&f); init=1; }
  QueryPerformanceCounter(&c); return (double)c.QuadPart/(double)f.QuadPart;
}
double TM[8];
static double cpun(void){
  FILETIME cr,ex,ke,us; ULARGE_INTEGER k,u;
  if(!GetProcessTimes(GetCurrentProcess(),&cr,&ex,&ke,&us)) return 0;
  k.LowPart=ke.dwLowDateTime; k.HighPart=ke.dwHighDateTime;
  u.LowPart=us.dwLowDateTime; u.HighPart=us.dwHighDateTime;
  return (k.QuadPart+u.QuadPart)*1e-7;
}
#else
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
double TM[8];
static double cpun(void){ struct timespec t; clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t); return t.tv_sec+1e-9*t.tv_nsec; }
#endif
static inline float hsum8s(const float t[8]){ return ((t[0]+t[1])+(t[2]+t[3]))+((t[4]+t[5])+(t[6]+t[7])); }
WLC_API int wlc_build_avx2(void){ return WLC_AVX2; }
WLC_API int wlc_cpu_avx2(void){
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#ifdef _MSC_VER
  int info[4]; __cpuid(info, 1);
  if((info[2] & (1<<27))==0 || (info[2] & (1<<28))==0) return 0;
  unsigned long long xcr = _xgetbv(0);
  if((xcr & 6ull)!=6ull) return 0;
  __cpuidex(info, 7, 0);
  return (info[1] & (1<<5))!=0;
#else
  unsigned int eax,ebx,ecx,edx;
  __asm__ __volatile__("cpuid":"=a"(eax),"=b"(ebx),"=c"(ecx),"=d"(edx):"a"(1),"c"(0));
  if((ecx&(1u<<27))==0 || (ecx&(1u<<28))==0) return 0;
  unsigned int xeax,xedx; unsigned int xcr_leaf=0;
  __asm__ __volatile__(".byte 0x0f,0x01,0xd0":"=a"(xeax),"=d"(xedx):"c"(xcr_leaf));
  if((xeax&6u)!=6u) return 0;
  __asm__ __volatile__("cpuid":"=a"(eax),"=b"(ebx),"=c"(ecx),"=d"(edx):"a"(7),"c"(0));
  return (ebx&(1u<<5))!=0;
#endif
#else
  return 0;
#endif
}
double PF[16]; static double ZR=0; WLC_API double wlc_zr(void){ return ZR; } static double RHO=0; static float MAXE=0;
#define MK(k) {double t2_=cpun(); PF[k]+=t2_-pt_; pt_=t2_;}
WLC_API void wlc_prof(double* o){ memcpy(o,PF,sizeof(PF)); memset(PF,0,sizeof(PF)); }

/* ---------------- range coder (LZMA-style, 12-bit probs from 2-rate mix) ---------------- */
typedef struct { uint64_t low; uint32_t range; uint8_t cache; uint64_t csize; uint8_t* out; long pos, cap; } Enc;
typedef struct { uint32_t code, range; const uint8_t* in; long pos, n; } Dec;
typedef struct { uint16_t a, b; uint16_t n, pad; } Prob;
static int RT[256]; static int rt_init=0;
static void rtinit(void){ if(rt_init) return; for(int i=0;i<256;i++) RT[i]=(int)(65536.0/(i+1.5)); rt_init=1; }
static void e_init(Enc* e, uint8_t* out, long cap){ e->low=0; e->range=0xFFFFFFFFu; e->cache=0; e->csize=1; e->out=out; e->pos=0; e->cap=cap; }
static inline void e_byte(Enc* e, uint8_t b){ if(e->pos<e->cap) e->out[e->pos]=b; e->pos++; }
static inline void e_shift(Enc* e){
  if((uint32_t)e->low < 0xFF000000u || (e->low>>32)!=0){ uint8_t carry=(uint8_t)(e->low>>32); uint8_t t=e->cache;
    do{ e_byte(e,(uint8_t)(t+carry)); t=0xFF; }while(--e->csize); e->cache=(uint8_t)(e->low>>24); }
  e->csize++; e->low=(e->low & 0x00FFFFFFu)<<8; }
static inline int pget(const Prob* p){ int v=(int)(((int)p->a+(int)p->b)>>1); if(v<12)v=12; if(v>65524)v=65524; return v; }
#ifndef NF
#define NF 20
#endif
#define NS 255
static inline void pupd(Prob* p,int bit){ int t= bit? 65535:0; int n=p->n; int rf= n<NF? RT[n]:RT[NF]; int rs=RT[n]; /* n<=255=NS */
  int a=p->a, b=p->b; a += (int)(((long)(t-a)*rf)>>16); b += (int)(((long)(t-b)*rs)>>16); p->a=(uint16_t)a; p->b=(uint16_t)b; if(n<255) p->n=n+1; }
#ifndef NWU
#define NWU 0
#endif
#ifndef P1R
#define P1R 0
#endif
static inline void pupd1(Prob* p,int bit){ int t= bit? 65535:0; int n=p->n; int rf= n<NF? RT[n]:RT[NF]; int a=p->a; a += (int)(((long)(t-a)*rf)>>16); p->a=p->b=(uint16_t)a; if(n<255) p->n=n+1; }
static inline void e_bitp(Enc* e,int bit,int p1){ uint32_t bound=(e->range>>16)*(uint32_t)p1;
  uint32_t mk=(uint32_t)0-(uint32_t)(bit==0); e->low+=bound&mk; e->range= (bound&~mk) | ((e->range-bound)&mk);
  while(e->range<(1u<<24)){ e->range<<=8; e_shift(e);} }
static inline void e_bit(Enc* e,Prob* p,int bit){ e_bitp(e,bit,pget(p)); pupd(p,bit); }
static long e_finish(Enc* e){ for(int i=0;i<5;i++) e_shift(e); return e->pos; }
static inline uint8_t d_byte(Dec* d){ return d->pos<d->n ? d->in[d->pos++] : 0; }
static void d_init(Dec* d,const uint8_t* in,long n){ d->in=in; d->n=n; d->pos=0; d->code=0; d->range=0xFFFFFFFFu; for(int i=0;i<5;i++) d->code=(d->code<<8)|d_byte(d); }
static inline int d_bitp(Dec* d,int p1){ uint32_t bound=(d->range>>16)*(uint32_t)p1; int bit;
  bit= d->code<bound; uint32_t mk=(uint32_t)0-(uint32_t)(bit==0); d->code-=bound&mk; d->range= (bound&~mk) | ((d->range-bound)&mk);
  while(d->range<(1u<<24)){ d->range<<=8; d->code=(d->code<<8)|d_byte(d);} return bit; }
static inline void e_bypass(Enc* e,unsigned val,int k){ while(k>0){ int kk=k>8?8:k; k-=kk; unsigned v=(val>>k)&((1u<<kk)-1); e->range>>=kk; e->low+=(uint64_t)v*e->range; while(e->range<(1u<<24)){ e->range<<=8; e_shift(e);} } }
static inline unsigned d_bypass(Dec* d,int k){ unsigned val=0; while(k>0){ int kk=k>8?8:k; k-=kk; d->range>>=kk; unsigned v=d->code/d->range; if(v>=(1u<<kk)) v=(1u<<kk)-1; d->code-=v*d->range; val=(val<<kk)|v; while(d->range<(1u<<24)){ d->range<<=8; d->code=(d->code<<8)|d_byte(d);} } return val; }
static inline int d_bit(Dec* d,Prob* p){ int b=d_bitp(d,pget(p)); pupd(p,b); return b; }
static void pinit(Prob* p,long n){ rtinit(); for(long i=0;i<n;i++){ p[i].a=32768; p[i].b=32768; p[i].n=0; } }
/* ---------------- CDF 9/7 lifting ---------------- */
#define LA (-1.586134342f)
#define LB (-0.05298011854f)
#define LG (0.8829110762f)
#define LD (0.4435068522f)
#define LZ (1.149604398f)
static void lift_fwd(float* x, float* t, int n){
  if(n<2) return; int i;
  for(i=1;i<n;i+=2){ float r=(i+1<n)?x[i+1]:x[i-1]; x[i]+=LA*(x[i-1]+r); }
  for(i=0;i<n;i+=2){ float l=(i>0)?x[i-1]:x[i+1]; float r=(i+1<n)?x[i+1]:x[i-1]; x[i]+=LB*(l+r); }
  for(i=1;i<n;i+=2){ float r=(i+1<n)?x[i+1]:x[i-1]; x[i]+=LG*(x[i-1]+r); }
  for(i=0;i<n;i+=2){ float l=(i>0)?x[i-1]:x[i+1]; float r=(i+1<n)?x[i+1]:x[i-1]; x[i]+=LD*(l+r); }
  int h=(n+1)/2; for(i=0;i<n;i+=2) t[i/2]=x[i]*(1.0f/LZ); for(i=1;i<n;i+=2) t[h+i/2]=x[i]*LZ;
  memcpy(x,t,sizeof(float)*n);
}
static void lift_inv(float* x, float* t, int n){
  if(n<2) return; int i; int h=(n+1)/2;
  for(i=0;i<n;i+=2) t[i]=x[i/2]*LZ; for(i=1;i<n;i+=2) t[i]=x[h+i/2]*(1.0f/LZ);
  memcpy(x,t,sizeof(float)*n);
  for(i=0;i<n;i+=2){ float l=(i>0)?x[i-1]:x[i+1]; float r=(i+1<n)?x[i+1]:x[i-1]; x[i]-=LD*(l+r); }
  for(i=1;i<n;i+=2){ float r=(i+1<n)?x[i+1]:x[i-1]; x[i]-=LG*(x[i-1]+r); }
  for(i=0;i<n;i+=2){ float l=(i>0)?x[i-1]:x[i+1]; float r=(i+1<n)?x[i+1]:x[i-1]; x[i]-=LB*(l+r); }
  for(i=1;i<n;i+=2){ float r=(i+1<n)?x[i+1]:x[i-1]; x[i]-=LA*(x[i-1]+r); }
}
/* split-array lifting: evens E[0..ne) and odds O[0..no) are rows of width w (vectorised over the row);
   arithmetic and operand order identical to the interleaved scalar version (bit-exact). */
#define TW 64
static void lrows(float* restrict E,float* restrict O,int ne,int no,int w,int inv){
  #define ODD(A,SG) for(int k=0;k<no;k++){ float* o=O+(long)k*w; const float* l=E+(long)k*w; const float* r=(k+1<ne)?E+(long)(k+1)*w:l; for(int j=0;j<w;j++) o[j] SG (A)*(l[j]+r[j]); }
  #define EVN(A,SG) for(int k=0;k<ne;k++){ float* e=E+(long)k*w; const float* l=O+(long)(k>0?k-1:0)*w; const float* r=(k<no)?O+(long)k*w:O+(long)(k-1)*w; for(int j=0;j<w;j++) e[j] SG (A)*(l[j]+r[j]); }
  if(!inv){ ODD(LA,+=) EVN(LB,+=) ODD(LG,+=) EVN(LD,+=) } else { EVN(LD,-=) ODD(LG,-=) EVN(LB,-=) ODD(LA,-=) }
  #undef ODD
  #undef EVN
}
static void l1d(float* restrict E,float* restrict O,int ne,int no,int inv){
  #define ODD1(A,SG) { int k=0; for(;k+1<ne&&k<no;k++) O[k] SG (A)*(E[k]+E[k+1]); for(;k<no;k++) O[k] SG (A)*(E[k]+E[k]); }
  #define EVN1(A,SG) { E[0] SG (A)*(O[0]+(no>0?O[0]:O[0])); int k=1; for(;k<no;k++) E[k] SG (A)*(O[k-1]+O[k]); for(;k<ne;k++) E[k] SG (A)*(O[k-1]+O[k-1]); }
  if(!inv){ ODD1(LA,+=) EVN1(LB,+=) ODD1(LG,+=) EVN1(LD,+=) } else { EVN1(LD,-=) ODD1(LG,-=) EVN1(LB,-=) ODD1(LA,-=) }
  #undef ODD1
  #undef EVN1
}
/* transform region [0,cz)x[0,cy)x[0,cx) of array (nz,ny,nx) along axis */
static void axis_pass(float* a,int ny,int nx,int cz,int cy,int cx,int axis,int inv){
  int n = axis==0?cz:(axis==1?cy:cx); if(n<2) return; int ne=(n+1)/2, no=n/2;
  long sz=(long)ny*nx; float* t=(float*)malloc(sizeof(float)*(long)n*TW);
  if(axis==2){
    for(int z=0;z<cz;z++) for(int y=0;y<cy;y++){ float* p=a+z*sz+(long)y*nx; float* E=t; float* O=t+ne;
      if(!inv){ for(int k=0;k<ne;k++) E[k]=p[2*k]; for(int k=0;k<no;k++) O[k]=p[2*k+1]; l1d(E,O,ne,no,0);
        for(int k=0;k<ne;k++) p[k]=E[k]*(1.0f/LZ); for(int k=0;k<no;k++) p[ne+k]=O[k]*LZ; }
      else { for(int k=0;k<ne;k++) E[k]=p[k]*LZ; for(int k=0;k<no;k++) O[k]=p[ne+k]*(1.0f/LZ); l1d(E,O,ne,no,1);
        for(int k=0;k<ne;k++) p[2*k]=E[k]; for(int k=0;k<no;k++) p[2*k+1]=O[k]; } }
  } else {
    long stride = axis==0? sz : nx; int outer = axis==0? cy : cz;
    for(int o=0;o<outer;o++) for(int x0=0;x0<cx;x0+=TW){ int w = cx-x0<TW? cx-x0: TW;
      float* base = a + (axis==0? (long)o*nx + x0 : (long)o*sz + x0); float* E=t; float* O=t+(long)ne*w;
      if(!inv){ for(int i=0;i<n;i++){ float* dst=(i&1)? O+(long)(i>>1)*w : E+(long)(i>>1)*w; memcpy(dst,base+i*stride,sizeof(float)*w); }
        lrows(E,O,ne,no,w,0);
        for(int k=0;k<ne;k++){ float* d=base+k*stride; const float* s_=E+(long)k*w; for(int j=0;j<w;j++) d[j]=s_[j]*(1.0f/LZ); }
        for(int k=0;k<no;k++){ float* d=base+(ne+k)*stride; const float* s_=O+(long)k*w; for(int j=0;j<w;j++) d[j]=s_[j]*LZ; } }
      else { for(int k=0;k<ne;k++){ const float* s_=base+k*stride; float* d=E+(long)k*w; for(int j=0;j<w;j++) d[j]=s_[j]*LZ; }
        for(int k=0;k<no;k++){ const float* s_=base+(ne+k)*stride; float* d=O+(long)k*w; for(int j=0;j<w;j++) d[j]=s_[j]*(1.0f/LZ); }
        lrows(E,O,ne,no,w,1);
        for(int i=0;i<n;i++){ const float* src=(i&1)? O+(long)(i>>1)*w : E+(long)(i>>1)*w; memcpy(base+i*stride,src,sizeof(float)*w); } }
    }
  }
  free(t);
}
static void dwt3_core(float* a,int nz,int ny,int nx,int Lz,int Ly,int Lx,int inv){
  int L = Lz>Ly?Lz:Ly; if(Lx>L) L=Lx;
  int dz[32],dy[32],dx[32]; int cz=nz,cy=ny,cx=nx;
  for(int l=0;l<L;l++){ dz[l]=cz;dy[l]=cy;dx[l]=cx; if(l<Lz&&cz>1) cz=(cz+1)/2; if(l<Ly&&cy>1) cy=(cy+1)/2; if(l<Lx&&cx>1) cx=(cx+1)/2; }
  if(!inv){ for(int l=0;l<L;l++){ if(l<Lx) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],2,0); if(l<Ly) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],1,0); if(l<Lz) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],0,0);} }
  else { for(int l=L-1;l>=0;l--){ if(l<Lz) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],0,1); if(l<Ly) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],1,1); if(l<Lx) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],2,1);} }
}
/* r8p: z-first wavelet packet (as in SPERR / JPEG2000 Part 2): PKZ levels of 1D DWT along z on the whole volume, z bands
   rescaled by sqrt(synthesis gain), then the 2D dyadic xy transform (Lz=0 path) on every plane. Selected per volume. */
static int PKZ=0; static double G1[2][16]; static void init_gains(void);
static void pk_scale(float* a,int nz,int ny,int nx,int pz,int inv){ init_gains(); long sz=(long)ny*nx; int c[32]; c[0]=nz; for(int l=0;l<pz;l++) c[l+1]=(c[l]+1)/2;
  for(int z=0;z<nz;z++){ double g; if(z<c[pz]) g=G1[0][pz]; else { int l=0; while(!(z>=c[l+1] && z<c[l])) l++; g=G1[1][l]; }
    float f=(float)(inv? 1.0/sqrt(g) : sqrt(g)); float* p=a+(long)z*sz; for(long i=0;i<sz;i++) p[i]*=f; } }
void dwt3_raw(float* a,int nz,int ny,int nx,int Lz,int Ly,int Lx,int inv){
  if(PKZ>0 && Lz==0 && nz>1){ if(!inv){ dwt3_core(a,nz,ny,nx,PKZ,0,0,0); pk_scale(a,nz,ny,nx,PKZ,0); dwt3_core(a,nz,ny,nx,0,Ly,Lx,0); }
    else { dwt3_core(a,nz,ny,nx,0,Ly,Lx,1); pk_scale(a,nz,ny,nx,PKZ,1); dwt3_core(a,nz,ny,nx,PKZ,0,0,1); } }
  else dwt3_core(a,nz,ny,nx,Lz,Ly,Lx,inv); }
void dwt3(float* a,int nz,int ny,int nx,int Lz,int Ly,int Lx,int inv){ double t_=cpun(); dwt3_raw(a,nz,ny,nx,Lz,Ly,Lx,inv); PF[0]+=cpun()-t_;  }
/* ---------------- subband boxes ---------------- */
typedef struct { int z0,z1,y0,y1,x0,x1,lev,o,parent; } Box;
int boxes(int nz,int ny,int nx,int Lz,int Ly,int Lx,int* out){ /* out: 9 ints per box, coding order */
  int L = Lz>Ly?Lz:Ly; if(Lx>L) L=Lx; Box b[256]; int nb=0; int cz=nz,cy=ny,cx=nx;
  Box lv[32][8]; int has[32][8]; memset(has,0,sizeof(has));
  for(int l=0;l<L;l++){
    int tz=(l<Lz&&cz>1), ty=(l<Ly&&cy>1), tx=(l<Lx&&cx>1);
    int lz=tz?(cz+1)/2:cz, ly=ty?(cy+1)/2:cy, lx=tx?(cx+1)/2:cx;
    for(int o=1;o<8;o++){ int hz=o&4,hy=o&2,hx=o&1; if((hz&&!tz)||(hy&&!ty)||(hx&&!tx)) continue;
      Box q={hz?lz:0,hz?cz:lz,hy?ly:0,hy?cy:ly,hx?lx:0,hx?cx:lx,l,o,-1};
      if(q.z1>q.z0&&q.y1>q.y0&&q.x1>q.x0){ lv[l][o]=q; has[l][o]=1; } }
    cz=lz;cy=ly;cx=lx; }
  Box ll={0,cz,0,cy,0,cx,L,0,-1}; b[nb++]=ll;
  int idx[32][8]; for(int l=L-1;l>=0;l--) for(int o=1;o<8;o++) if(has[l][o]){ Box q=lv[l][o]; q.parent=(l+1<L&&has[l+1][o])?idx[l+1][o]:-1; idx[l][o]=nb; b[nb++]=q; }
  for(int i=0;i<nb;i++){ int* p=out+9*i; p[0]=b[i].z0;p[1]=b[i].z1;p[2]=b[i].y0;p[3]=b[i].y1;p[4]=b[i].x0;p[5]=b[i].x1;p[6]=b[i].lev;p[7]=b[i].o;p[8]=b[i].parent; }
  return nb;
}
/* ---------------- coefficient coder: 3-model logistic mixing ---------------- */
#define NBK 13
#define NPB 9
#define NCLS 16
#define NMB 14
#define QMAX 12
#ifndef MIXG
#define MIXG 8
#endif
#ifndef QZ
#define QZ 0
#endif
#ifndef MIXB
#define MIXB MIXG
#endif
#ifndef MIXT
#define MIXT 1024 /* R8: zero flag coded by the primary model alone when it is already confident */
#endif
#ifndef RDOA
#define RDOA 3
#endif
#ifndef MIXJ
#define MIXJ 1
#endif
#define NFB 24
static int bucket(int s){ static const int th[NBK-1]={1,2,3,5,7,10,14,20,28,40,60,90}; int b=0; while(b<NBK-1 && s>=th[b]) b++; return b; }
static int fbucket(int s){ static const int th[NFB-1]={1,2,3,4,5,6,8,10,12,15,18,22,27,33,40,50,62,77,95,120,150,200,280}; int b=0; while(b<NFB-1 && s>=th[b]) b++; return b; }
typedef struct {
  Prob zf[NCLS][NBK][NPB]; Prob z2[NCLS][27][9]; Prob z3[NCLS][NFB][3];
  Prob mg[NCLS][NBK][NMB]; Prob mg2[NCLS][NPB][NMB]; Prob mg3[NCLS][27][6];
  Prob sg[NCLS][27]; Prob eg[NCLS][24]; Prob mk[NCLS][16][QMAX];
  int zw[NCLS][NBK][4]; int mw[NCLS][NMB][4]; int ra[NCLS][NBK][2];
} Model;
static int clsof(int lev,int o,int L){ if(o==0) return 15; int pc=((o>>2)&1)+((o>>1)&1)+(o&1); int l=lev>4?4:lev; return l*3+pc-1; }
/* anisotropic data (GZW==1): coding classes also split z-highpass vs not (offsets/dequant keep the 16 base classes) */
#ifndef XCLS
#define XCLS 0 /* tried: -0.1% on hurr_TC, needs NCLS 32 */
#endif
static int clsx(int lev,int o,int L){ if(o==0) return 31; int pc=((o>>2)&1)+((o>>1)&1)+(o&1); int zh=(o>>2)&1; int l=lev>4?4:lev;
  int id= pc==3? 4 : (pc-1)*2+(zh?0:1); return l*5+id; }
static inline int sgn(int v){ return v>0?1:(v<0?2:0); }
typedef struct { int nbk,pbk,sc,nb3,fb,cs3,pq; } Ctx;
typedef struct { int hasP, ncz; long *pz,*py,*px; long *cz[64][3]; long* mem; uint16_t* cs; int dy,dx; } BT;
static uint8_t BKL[8192], FBL[8192]; static int bkl_init=0;
static void bkl(void){ if(bkl_init) return; for(int s=0;s<8192;s++){ BKL[s]=(uint8_t)bucket(s); FBL[s]=(uint8_t)fbucket(s); } bkl_init=1; }
#ifndef CSUM
#define CSUM 1
#endif
static void bt_make(BT* T,const int* B,const int* P,const int* cz,int ncz,long sz,int nx,const uint8_t* M_){
  int dz=B[1]-B[0], dy=B[3]-B[2], dx=B[5]-B[4]; long per=(long)dz+dy+dx; T->mem=(long*)malloc(sizeof(long)*per*(1+ncz)); long* p=T->mem;
  T->hasP= P!=NULL; T->ncz=ncz; T->pz=p; T->py=p+dz; T->px=p+dz+dy; p+=per;
  if(P){ for(int z=0;z<dz;z++) T->pz[z]=(long)(P[0]+(int)((long)z*(P[1]-P[0])/dz))*sz; for(int y=0;y<dy;y++) T->py[y]=(long)(P[2]+(int)((long)y*(P[3]-P[2])/dy))*nx; for(int x=0;x<dx;x++) T->px[x]=P[4]+(int)((long)x*(P[5]-P[4])/dx); }
  for(int j=0;j<ncz;j++){ const int* C=cz+9*j; T->cz[j][0]=p; T->cz[j][1]=p+dz; T->cz[j][2]=p+dz+dy; p+=per;
    for(int z=0;z<dz;z++){ int zz=C[0]+z; if(zz>=C[1]) zz=C[1]-1; T->cz[j][0][z]=(long)zz*sz; }
    for(int y=0;y<dy;y++){ int yy=C[2]+y; if(yy>=C[3]) yy=C[3]-1; T->cz[j][1][y]=(long)yy*nx; }
    for(int x=0;x<dx;x++){ int xx=C[4]+x; if(xx>=C[5]) xx=C[5]-1; T->cz[j][2][x]=xx; } }
  T->cs=NULL; T->dy=dy; T->dx=dx;
#if CSUM
  if(ncz>0 && M_){ /* R8: cousin magnitude sum precomputed per band (identical value, one load per voxel instead of 3*ncz) */
    T->cs=(uint16_t*)calloc((long)dz*dy*dx,sizeof(uint16_t));
    for(int j=0;j<ncz;j++){ const int* C=cz+9*j; int xl=C[5]-C[4]; if(xl>dx) xl=dx;
      for(int z=0;z<dz;z++) for(int y=0;y<dy;y++){ uint16_t* o=T->cs+((long)z*dy+y)*dx; const uint8_t* r=M_+T->cz[j][0][z]+T->cz[j][1][y];
        for(int x=0;x<xl;x++) o[x]+=r[C[4]+x]; for(int x=xl;x<dx;x++) o[x]+=r[T->cz[j][2][x]]; } } }
#endif
}
/* anisotropy-adaptive neighbourhood (R7): weight of the previous-plane (z) neighbours in the magnitude context.
   Chosen by the encoder from the log-magnitude correlation of neighbouring finest-scale HHH coefficients along z vs x,y;
   e.g. hurricane fields (z = vertical levels) have uncorrelated z-detail magnitudes -> weight 1 and no z-diagonals. Stored in the header. */
static int GZW=2;
#ifndef ZTH
#define ZTH 0.6
#endif
static double lcorr(const double* a,const double* b,long n){ double sa=0,sb=0,saa=0,sbb=0,sab=0; for(long i=0;i<n;i++){ sa+=a[i]; sb+=b[i]; } sa/=n; sb/=n;
  for(long i=0;i<n;i++){ double x=a[i]-sa, y=b[i]-sb; saa+=x*x; sbb+=y*y; sab+=x*y; } return (saa>0&&sbb>0)? sab/sqrt(saa*sbb) : 0; }
static int zw_choose(const void* vp,int isf,int nz,int ny,int nx){
  if(nz<12||ny<12||nx<12) return 2; long sz=(long)ny*nx; const float* vf=(const float*)vp; const uint16_t* vu=(const uint16_t*)vp;
  #define VAL(i) (isf? (double)vf[i] : (double)vu[i])
  int nzs=(nz-6)/2, nys=(ny-6)/2, nxs=(nx-6)/2; double tot=(double)nzs*nys*nxs; int st=1; while(tot/((double)st*st*st)>40000) st++;
  long cap=(long)(nzs/st+1)*(nys/st+1)*(nxs/st+1); double* A=(double*)malloc(sizeof(double)*cap*4); long n=0;
  for(int z=1;z+4<nz-1;z+=2*st) for(int y=1;y+4<ny-1;y+=2*st) for(int x=1;x+4<nx-1;x+=2*st){
    double dv[4]; int pz[4]={z,z+2,z,z}, py[4]={y,y,y+2,y}, px[4]={x,x,x,x+2};
    for(int t=0;t<4;t++){ double acc=0; for(int a=-1;a<=1;a++) for(int b=-1;b<=1;b++) for(int c=-1;c<=1;c++){ double wgt=(a?-0.5:1.0)*(b?-0.5:1.0)*(c?-0.5:1.0);
          acc+=wgt*VAL((long)(pz[t]+a)*sz+(long)(py[t]+b)*nx+(px[t]+c)); } dv[t]=log(fabs(acc)+1e-30); }
    if(n<cap){ for(int t=0;t<4;t++) A[n*4+t]=dv[t]; n++; } }
  #undef VAL
  double* a=(double*)malloc(sizeof(double)*n*4); for(long i=0;i<n;i++){ a[i]=A[i*4]; a[n+i]=A[i*4+1]; a[2*n+i]=A[i*4+2]; a[3*n+i]=A[i*4+3]; }
  double cz=lcorr(a,a+n,n), cy=lcorr(a,a+2*n,n), cx=lcorr(a,a+3*n,n); free(A); free(a);
  double cxy=0.5*(cx+cy); ZR=cxy>0.05? cz/cxy : 1; return (cxy>0.05 && cz<ZTH*cxy)? 1 : 2; }
typedef struct { int hz,hy,ey,ym; long pzy; const uint16_t* csr; } RowC;
static inline void ctx_row(const uint8_t* m,const int32_t* q,int ny,int nx,const int* B,const BT* T,long i,int x,const RowC* R,Ctx* C){
  long sz=(long)ny*nx; int s=0,w=0,nn=0,mw=0,mn=0,mu=0;
  int hx=x>B[4], hy=R->hy, hz=R->hz, ex=x+1<B[5], ey=R->ey;
  if(hx){ mw=m[i-1]; w=q[i-1]; } if(hy){ mn=m[i-nx]; nn=q[i-nx]; } if(hz) mu=m[i-sz];
  s=2*(mw+mn)+GZW*mu;
  if(hy){ if(hx) s+=m[i-nx-1]; if(ex) s+=m[i-nx+1]; }
  if(hz && GZW>=2){ if(hx) s+=m[i-sz-1]; if(ex) s+=m[i-sz+1]; if(hy) s+=m[i-sz-nx]; if(ey) s+=m[i-sz+nx]; }
  if(x-1>B[4]) s+=m[i-2]; if(R->ym) s+=m[i-2*nx];
  C->nbk=BKL[(s+1)/2];
  int lx=x-B[4];
  int pb=0,pm=0,ps=0; if(R->pzy>=0){ long pi=R->pzy+T->px[lx]; pm=m[pi]; ps=sgn(q[pi]); pb=pm==0?0:(pm==1?1:2); }
  int cs= R->csr? R->csr[lx] : 0;
  C->cs3=cs==0?0:(cs<=2?1:2); C->pbk=pb*3+C->cs3; C->sc=(sgn(w)*3+sgn(nn))*3+ps; C->pq=pb;
  C->nb3=(mw>2?2:mw)*9+(mn>2?2:mn)*3+(mu>2?2:mu);
  { int t=s+3*pm+cs; C->fb=FBL[t<8191?t:8191]; }
}
static inline void row_init(const int* B,const BT* T,int z,int y,RowC* R){ int lz=z-B[0], ly=y-B[2];
  R->hz=z>B[0]; R->hy=y>B[2]; R->ey=y+1<B[3]; R->ym=y-1>B[2]; R->pzy= T->hasP? T->pz[lz]+T->py[ly] : -1; R->csr= T->cs? T->cs+((long)lz*T->dy+ly)*T->dx : NULL; }
static inline void ctx_at(const uint8_t* m,const int32_t* q,int ny,int nx,const int* B,const BT* T,int z,int y,int x,Ctx* C){
  long sz=(long)ny*nx; long i=(long)z*sz+(long)y*nx+x; int s=0,w=0,nn=0,mw=0,mn=0,mu=0;
  int hx=x>B[4], hy=y>B[2], hz=z>B[0], ex=x+1<B[5], ey=y+1<B[3];
  if(hx){ mw=m[i-1]; w=q[i-1]; } if(hy){ mn=m[i-nx]; nn=q[i-nx]; } if(hz) mu=m[i-sz];
  s=2*(mw+mn)+GZW*mu;
  if(hx&&hy) s+=m[i-nx-1]; if(ex&&hy) s+=m[i-nx+1];
  if(hz && GZW>=2){ if(hx) s+=m[i-sz-1]; if(ex) s+=m[i-sz+1]; if(hy) s+=m[i-sz-nx]; if(ey) s+=m[i-sz+nx]; }
  if(x-1>B[4]) s+=m[i-2]; if(y-1>B[2]) s+=m[i-2*nx];
  C->nbk=BKL[(s+1)/2];
  int lz=z-B[0], ly=y-B[2], lx=x-B[4];
  int pb=0,pm=0,ps=0; if(T->hasP){ long pi=T->pz[lz]+T->py[ly]+T->px[lx]; pm=m[pi]; ps=sgn(q[pi]); pb=pm==0?0:(pm==1?1:2); }
  int cs=0; if(T->cs) cs=T->cs[((long)lz*T->dy+ly)*T->dx+lx]; else for(int j=0;j<T->ncz;j++) cs+=m[T->cz[j][0][lz]+T->cz[j][1][ly]+T->cz[j][2][lx]];
  C->cs3=cs==0?0:(cs<=2?1:2); C->pbk=pb*3+C->cs3; C->sc=(sgn(w)*3+sgn(nn))*3+ps; C->pq=pb;
  C->nb3=(mw>2?2:mw)*9+(mn>2?2:mn)*3+(mu>2?2:mu);
  { int t=s+3*pm+cs; C->fb=FBL[t<8191?t:8191]; }
}
static float CT[4096]; static int ct_init=0; static short STR[4096]; static short SQ[4096];
static void ctinit(void){ if(ct_init) return; for(int i=0;i<4096;i++){ CT[i]=(float)(-log2((i*16+8)/65536.0));
    double p=(i+0.5)/4096.0; int v=(int)lrint(256*log(p/(1-p))); if(v>2047)v=2047; if(v<-2047)v=-2047; STR[i]=(short)v;
    double x=(i-2048)/256.0; int u=(int)lrint(4096/(1+exp(-x))); if(u<1)u=1; if(u>4095)u=4095; SQ[i]=(short)u; } ct_init=1; }
static inline float cbit(const Prob* p,int bit){ int v=pget(p); return bit? CT[v>>4] : CT[(65536-v)>>4]; }
#ifndef LR
#define LR 6
#endif
#ifndef MS
#define MS 14
#endif
static inline int mixc(Enc* e,Dec* d,Prob* p0,Prob* p1,Prob* p2,int* W,int bit){
  /* STR index = clamp(avg)>>4 == avg>>4 (the clamp [12,65524] never changes the top 12 bits) */
  int s0=STR[((int)p0->a+(int)p0->b)>>5], s1=STR[((int)p1->a+(int)p1->b)>>5], s2=STR[((int)p2->a+(int)p2->b)>>5];
  long dot=(long)W[0]*s0+(long)W[1]*s1+(long)W[2]*s2+(long)W[3]*256; int x=(int)(dot>>16); if(x>2047)x=2047; if(x<-2047)x=-2047;
  int p12=SQ[x+2048]; int p16=p12<<4; if(p16<32)p16=32; if(p16>65504)p16=65504;
  if(e) e_bitp(e,bit,p16); else bit=d_bitp(d,p16);
#if NWU
  (void)p12;
#else
  int err=((bit<<12)-p12)*LR; W[0]+=(s0*err)>>MS; W[1]+=(s1*err)>>MS; W[2]+=(s2*err)>>MS; W[3]+=(256*err)>>MS;
#endif
#if P1R
  pupd(p0,bit); pupd1(p1,bit); pupd1(p2,bit); return bit;
#else
  pupd(p0,bit); pupd(p1,bit); pupd(p2,bit); return bit;
#endif
}
#ifndef Z2X
#define Z2X 1
#endif
#ifndef MB2
#define MB2 0
#endif
#ifndef ZM
#define ZM 3
#endif
/* R8 2-input mixer (fewer models on the hot path): p0 + one secondary */
static inline int mix2(Enc* e,Dec* d,Prob* p0,Prob* p1,int* W,int bit){
  int s0=STR[((int)p0->a+(int)p0->b)>>5], s1=STR[((int)p1->a+(int)p1->b)>>5];
  long dot=(long)W[0]*s0+(long)W[1]*s1+(long)W[3]*256; int x=(int)(dot>>16); if(x>2047)x=2047; if(x<-2047)x=-2047;
  int p12=SQ[x+2048]; int p16=p12<<4; if(p16<32)p16=32; if(p16>65504)p16=65504;
  if(e) e_bitp(e,bit,p16); else bit=d_bitp(d,p16);
  int err=((bit<<12)-p12)*LR; W[0]+=(s0*err)>>MS; W[1]+=(s1*err)>>MS; W[3]+=(256*err)>>MS;
  pupd(p0,bit); pupd(p1,bit); return bit; }
static void model_init(Model* M){ pinit((Prob*)M,(long)(offsetof(Model,zw)/sizeof(Prob)));
  for(int c=0;c<NCLS;c++){ for(int k=0;k<NBK;k++){ M->zw[c][k][0]=26000; M->zw[c][k][1]=20000; M->zw[c][k][2]=20000; M->zw[c][k][3]=0; }
    for(int k=0;k<NBK;k++){ M->ra[c][k][0]=1; M->ra[c][k][1]=1; }
    for(int k=0;k<NMB;k++){ M->mw[c][k][0]=26000; M->mw[c][k][1]=20000; M->mw[c][k][2]=20000; M->mw[c][k][3]=0; } } }
/* unified coder: e!=NULL encode (q may be modified by RDO), else decode into q.
   Each subband is tiled in 4x4x4 cells; an 'all-zero cell' flag (ctx: class, coded neighbour cells, parent, cousins)
   skips empty cells; voxels inside non-empty cells are context-mixed coded. */
#define CS 4
#ifndef CCF
#define CCF 1
#endif
typedef struct { Prob cf[NCLS][4][3][2]; Prob sf[NCLS][4][3][2][2]; } CModel;
#ifndef SUB
#define SUB 0
#endif
#ifndef RCAP
#define RCAP 0
#endif
#ifndef RKB
#define RKB 1
#endif
static inline int ricek(const int* ra){ int k=0; while((ra[1]<<(k+RKB))<ra[0] && k<15) k++; return k; }
static inline float magbits(Model* M,int c,int nbk,int k,int r){ /* approx bits of magnitude r=|v|-1 (RDO) */
  float b=cbit(&M->mg[c][nbk][0],r>0); if(r<=0) return b; int r1=r-1, q=r1>>k; b+=k;
  for(int j=0;j<QMAX;j++){ int bit=q>j; b+=cbit(k? &M->mk[c][k][j] : &M->mg[c][nbk][j+1],bit); if(!bit) return b; }
  return b+2*log2f((float)(q-QMAX+1))+1; }
static inline void code_one(int32_t* q,uint8_t* m,Model* M,Enc* e,Dec* d,const float* cn,const float* off,float lam,int c,int cb,const int* B,const BT* T,int z,int y,int x,int ny,int nx,const RowC* RC){
  long sz=(long)ny*nx; long i=(long)z*sz+(long)y*nx+x; Ctx C; if(RC) ctx_row(m,q,ny,nx,B,T,i,x,RC,&C); else ctx_at(m,q,ny,nx,B,T,z,y,x,&C);
  int v= e? q[i] : 0; int* RA=M->ra[c][C.nbk]; int k=0;
  if(e && cn && v && v<=RDOA && v>=-RDOA){ k=ricek(RA); int a=v<0?-v:v; float aa=fabsf(cn[i]); float o=off[cb];
    float b_a=cbit(&M->zf[c][C.nbk][C.pbk],1)+1.0f; float b_m=b_a; int ra=a-1, rm=a-2;
    b_a+=magbits(M,c,C.nbk,k,ra);
    if(rm<0) b_m=cbit(&M->zf[c][C.nbk][C.pbk],0); else b_m+=magbits(M,c,C.nbk,k,rm);
    float da=aa-(a+o), dm= rm<0? aa : aa-(a-1+o);
    if(dm*dm+lam*b_m < da*da+lam*b_a){ v= v<0? v+1 : v-1; q[i]=v; } }
  int nzf; if(C.nbk>=MIXG || (QZ && C.nbk==0 && C.pbk==0) || (MIXT && (unsigned)(pget(&M->zf[c][C.nbk][C.pbk])-MIXT) > (unsigned)(65536-2*MIXT))){ Prob* pz=&M->zf[c][C.nbk][C.pbk]; if(e){ nzf=v!=0; e_bit(e,pz,nzf);} else nzf=d_bit(d,pz); }
  else if(ZM==3) nzf=mix2(e,d,&M->zf[c][C.nbk][C.pbk],&M->z2[c][C.nb3][Z2X? C.pbk : C.pq],M->zw[c][C.nbk],v!=0);
  else if(ZM==5) nzf=mix2(e,d,&M->zf[c][C.nbk][C.pbk],&M->z3[c][C.fb][C.cs3],M->zw[c][C.nbk],v!=0);
  else nzf=mixc(e,d,&M->zf[c][C.nbk][C.pbk],&M->z2[c][C.nb3][C.pq],&M->z3[c][C.fb][C.cs3],M->zw[c][C.nbk],v!=0);
  if(!nzf) return;
  int neg; Prob* ps=&M->sg[c][C.sc]; if(e){ neg=v<0; e_bit(e,ps,neg); } else neg=d_bit(d,ps);
  int r= e? (v<0?-v:v)-1 : 0;
  int b0; if(C.nbk>=MIXB){ Prob* pb=&M->mg[c][C.nbk][0]; if(e){ b0=r>0; e_bit(e,pb,b0);} else b0=d_bit(d,pb); }
  else if(MB2) b0=mix2(e,d,&M->mg[c][C.nbk][0],&M->mg2[c][C.pbk][0],M->mw[c][0], e? (r>0):0);
  else b0=mixc(e,d,&M->mg[c][C.nbk][0],&M->mg2[c][C.pbk][0],&M->mg3[c][C.nb3][0],M->mw[c][0], e? (r>0):0);
  if(b0){ k=ricek(RA); int r1= e? r-1 : 0; int qq= r1>>k; int j=0;
    for(;j<QMAX;j++){ Prob* pp= k? &M->mk[c][k][j] : &M->mg[c][C.nbk][j+1]; int bit; if(e){ bit=qq>j; e_bit(e,pp,bit);} else bit=d_bit(d,pp); if(!bit) break; }
    if(j==QMAX){ if(e){ unsigned u=(unsigned)(qq-QMAX)+1; int nbits=0; while((u>>nbits)>1) nbits++;
        for(int t=0;t<nbits;t++) e_bit(e,&M->eg[c][t],1); e_bit(e,&M->eg[c][nbits],0); e_bypass(e,u,nbits); }
      else { int nbits=0; while(d_bit(d,&M->eg[c][nbits])) nbits++; unsigned u=(1u<<nbits)|d_bypass(d,nbits); qq=QMAX+(int)(u-1); } }
    else if(d) qq=j;
    if(k){ if(e) e_bypass(e,(unsigned)(r1&((1<<k)-1)),k); else r1=(qq<<k)|(int)d_bypass(d,k); } else if(d) r1=qq;
    r=r1+1; { int cap_= RCAP? ((RA[0]/RA[1])+1)*RCAP : 1<<30; RA[0]+= r1<cap_? r1 : cap_; } RA[1]++; if(RA[1]>=64){ RA[0]>>=1; RA[1]>>=1; } }
  else r=0;
  int a=r+1; m[i]=a>255?255:a; if(d) q[i]=neg?-a:a;
}
static long EABORT=0; /* r10: >0 -> stop coefficient coding once the stream passes EABORT bytes (palette already smaller) */
#define ABT (e && EABORT && e->pos>EABORT)
static long coef_code_raw(int32_t* q,int nz,int ny,int nx,int Lz,int Ly,int Lx,Enc* e,Dec* d,const float* cn,const float* off,float lam){
  ctinit(); bkl(); int bx[256*9]; int nb=boxes(nz,ny,nx,Lz,Ly,Lx,bx); uint8_t* CFS[256]; int CFG[256][3]; for(int k=0;k<256;k++) CFS[k]=NULL; int L=Lz>Ly?Lz:Ly; if(Lx>L)L=Lx;
  long n=(long)nz*ny*nx; uint8_t* m=(uint8_t*)memset(balloc(n),0,n); Model* M=(Model*)malloc(sizeof(Model)); model_init(M);
  CModel* CM=(CModel*)malloc(sizeof(CModel)); pinit((Prob*)CM,sizeof(CModel)/sizeof(Prob));
  if(d) memset(q,0,sizeof(int32_t)*n); long sz=(long)ny*nx;
  for(int k=0;k<nb && !ABT;k++){ int* B=bx+9*k; int* P=B[8]>=0? bx+9*B[8]:NULL; int cb_=clsof(B[6],B[7],L); int c= (GZW==1 && XCLS)? clsx(B[6],B[7],L) : cb_;
    int cz[64*9]; int czi[64]; int ncz=0; for(int j=0;j<k;j++){ int* Cb=bx+9*j; if(Cb[7]!=0 && Cb[6]==B[6] && B[7]!=0){ memcpy(cz+9*ncz,Cb,9*sizeof(int)); czi[ncz]=j; ncz++; } }
    BT T; bt_make(&T,B,P,cz,ncz,sz,nx,m);
    if(B[7]==0){ for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++) for(int x=B[4];x<B[5];x++) code_one(q,m,M,e,d,cn,off,lam,c,cb_,B,&T,z,y,x,ny,nx,NULL); free(T.mem); continue; }
    int gz=(B[1]-B[0]+CS-1)/CS, gy=(B[3]-B[2]+CS-1)/CS, gx=(B[5]-B[4]+CS-1)/CS; uint8_t* cf=(uint8_t*)calloc((long)gz*gy*gx,1); CFS[k]=cf; CFG[k][0]=gz; CFG[k][1]=gy; CFG[k][2]=gx;
    for(int Z=0;Z<gz && !ABT;Z++) for(int Y=0;Y<gy;Y++) for(int X=0;X<gx;X++){ long ci=((long)Z*gy+Y)*gx+X;
      int z0=B[0]+Z*CS, y0=B[2]+Y*CS, x0=B[4]+X*CS, z1=z0+CS<B[1]?z0+CS:B[1], y1=y0+CS<B[3]?y0+CS:B[3], x1=x0+CS<B[5]?x0+CS:B[5];
      int nbf=(X>0&&cf[ci-1])+(Y>0&&cf[ci-gx])+(Z>0&&cf[ci-(long)gy*gx]);
      int pmx=0; if(P){ for(int z=z0;z<z1;z+=2) for(int y=y0;y<y1;y+=2){ long pb_=T.pz[z-B[0]]+T.py[y-B[2]]; for(int x=x0;x<x1;x+=2){ int pm=m[pb_+T.px[x-B[4]]]; if(pm>pmx) pmx=pm; } } }
      int csx=0;
#if CCF
      for(int j=0;j<ncz;j++){ int jj=czi[j]; const int* G=CFG[jj]; if(CFS[jj] && Z<G[0] && Y<G[1] && X<G[2] && CFS[jj][((long)Z*G[1]+Y)*G[2]+X]){ csx=1; break; } }
#else
      for(int j=0;j<ncz&&!csx;j++){ const int* Cb=cz+9*j; for(int z=z0;z<z1&&!csx;z++){ int zz=Cb[0]+(z-B[0]); if(zz>=Cb[1]) break; for(int y=y0;y<y1&&!csx;y++){ int yy=Cb[2]+(y-B[2]); if(yy>=Cb[3]) break; for(int x=x0;x<x1;x++){ int xx=Cb[4]+(x-B[4]); if(xx>=Cb[5]) break; if(m[(long)zz*sz+(long)yy*nx+xx]){csx=1;break;} } } } }
#endif
      Prob* pf=&CM->cf[c][nbf][pmx==0?0:(pmx==1?1:2)][csx];
      int has=0;
      if(e){ if(cn && lam>0){ /* RDO may zero the cell later; decide on current q */ }
        for(int z=z0;z<z1&&!has;z++) for(int y=y0;y<y1&&!has;y++) for(int x=x0;x<x1;x++) if(q[(long)z*sz+(long)y*nx+x]){has=1;break;}
        e_bit(e,pf,has); } else has=d_bit(d,pf);
      cf[ci]=has; if(!has) continue;
#if SUB
      { /* 2x2x2 sub-cell significance inside a non-empty cell: skips coding the zeros of empty sub-cells */
        int nsub=((z1-z0+1)/2)*((y1-y0+1)/2)*((x1-x0+1)/2), si=0, any=0;
        for(int sz_=z0;sz_<z1;sz_+=2) for(int sy=y0;sy<y1;sy+=2) for(int sx=x0;sx<x1;sx+=2){ si++;
          int ez=sz_+2<z1?sz_+2:z1, ey=sy+2<y1?sy+2:y1, ex=sx+2<x1?sx+2:x1; int f;
          if(si==nsub && !any) f=1;
          else { int nb=0; long i0=(long)sz_*sz+(long)sy*nx+sx;
            for(int y=sy;y<ey;y++) for(int x=sx;x<ex;x++){ if(sz_>B[0]) nb+= m[(long)(sz_-1)*sz+(long)y*nx+x]!=0; }
            for(int z=sz_;z<ez;z++) for(int x=sx;x<ex;x++){ if(sy>B[2]) nb+= m[(long)z*sz+(long)(sy-1)*nx+x]!=0; }
            for(int z=sz_;z<ez;z++) for(int y=sy;y<ey;y++){ if(sx>B[4]) nb+= m[(long)z*sz+(long)y*nx+sx-1]!=0; }
            int nbb= nb==0?0:(nb<=2?1:(nb<=5?2:3)); int pb_=0;
            if(P){ int pm=m[T.pz[sz_-B[0]]+T.py[sy-B[2]]+T.px[sx-B[4]]]; pb_= pm==0?0:(pm==1?1:2); }
            int cs_=0; for(int j=0;j<ncz;j++){ if(m[T.cz[j][0][sz_-B[0]]+T.cz[j][1][sy-B[2]]+T.cz[j][2][sx-B[4]]]){ cs_=1; break; } }
            (void)i0; Prob* ps_=&CM->sf[c][nbb][pb_][cs_][any?1:0];
            if(e){ f=0; for(int z=sz_;z<ez&&!f;z++) for(int y=sy;y<ey&&!f;y++) for(int x=sx;x<ex;x++) if(q[(long)z*sz+(long)y*nx+x]){f=1;break;} e_bit(e,ps_,f); } else f=d_bit(d,ps_); }
          if(!f) continue; any=1;
          for(int z=sz_;z<ez;z++) for(int y=sy;y<ey;y++) for(int x=sx;x<ex;x++) code_one(q,m,M,e,d,cn,off,lam,c,cb_,B,&T,z,y,x,ny,nx,NULL); } }
#else
      for(int z=z0;z<z1;z++) for(int y=y0;y<y1;y++){ RowC RC; row_init(B,&T,z,y,&RC); for(int x=x0;x<x1;x++) code_one(q,m,M,e,d,cn,off,lam,c,cb_,B,&T,z,y,x,ny,nx,&RC); }
#endif
    }
    free(T.cs); free(T.mem); }
  for(int k=0;k<nb;k++) free(CFS[k]); free(m); free(M); free(CM); return 0;
}
static long coef_code(int32_t* q,int nz,int ny,int nx,int Lz,int Ly,int Lx,Enc* e,Dec* d,const float* cn,const float* off,float lam){ double t_=cpun(); long r_=coef_code_raw(q,nz,ny,nx,Lz,Ly,Lx,e,d,cn,off,lam); PF[1]+=cpun()-t_; return r_; }
long enc_coefs(int32_t* q,int nz,int ny,int nx,int Lz,int Ly,int Lx,uint8_t* out,long cap,const float* cn,float id,const float* off,float lam){
  Enc e; e_init(&e,out,cap); coef_code(q,nz,ny,nx,Lz,Ly,Lx,&e,NULL,cn,off,lam); return e_finish(&e); }
int dec_coefs(const uint8_t* in,long len,int nz,int ny,int nx,int Lz,int Ly,int Lx,int32_t* q){
  Dec d; d_init(&d,in,len); coef_code(q,nz,ny,nx,Lz,Ly,Lx,NULL,&d,NULL,NULL,0); return 0; }
/* ---------------- sparse outlier coder with 8^3 block flags ---------------- */
typedef struct { Prob bf[4][3]; Prob f[6][4]; Prob s; Prob mg[8]; Prob eg[24]; } OModel;
static void make_ec_raw(const float* b,int nz,int ny,int nx,float tau,uint8_t* ec){
  long sz=(long)ny*nx; long n=(long)nz*sz; memset(ec,0,n); float t= tau>0? tau : 1.0f; float th[5]={0.5f*t,t,2*t,4*t,8*t};
  for(int z=1;z<nz-1;z++) for(int y=1;y<ny-1;y++){ long i0=(long)z*sz+(long)y*nx; const float* p=b+i0;
    for(int x=1;x<nx-1;x++){ float c2=2*p[x]; float a=fabsf(c2-p[x-1]-p[x+1])+fabsf(c2-p[x-nx]-p[x+nx])+fabsf(c2-p[x-sz]-p[x+sz]);
      int k=(a>=th[0])+(a>=th[1])+(a>=th[2])+(a>=th[3])+(a>=th[4]); ec[i0+x]=(uint8_t)k; } } }
#define OB 8
static long out_code_raw(int32_t* k,const uint8_t* ec,int nz,int ny,int nx,Enc* e,Dec* d){
  OModel M; pinit((Prob*)&M,sizeof(M)/sizeof(Prob)); long sz=(long)ny*nx;
  int bz=(nz+OB-1)/OB, by=(ny+OB-1)/OB, bxn=(nx+OB-1)/OB; uint8_t* bf=(uint8_t*)calloc((long)bz*by*bxn,1);
  for(int Z=0;Z<bz;Z++) for(int Y=0;Y<by;Y++) for(int X=0;X<bxn;X++){ long bi=((long)Z*by+Y)*bxn+X;
    int has=0;
    if(e){ for(int z=Z*OB;z<nz&&z<Z*OB+OB&&!has;z++) for(int y=Y*OB;y<ny&&y<Y*OB+OB&&!has;y++) for(int x=X*OB;x<nx&&x<X*OB+OB;x++) if(k[(long)z*sz+(long)y*nx+x]){has=1;break;} }
    int nb=(X>0&&bf[bi-1])+(Y>0&&bf[bi-bxn])+(Z>0&&bf[bi-(long)by*bxn]);
    int mx=0; for(int z=Z*OB;z<nz&&z<Z*OB+OB;z++) for(int y=Y*OB;y<ny&&y<Y*OB+OB;y++) for(int x=X*OB;x<nx&&x<X*OB+OB;x++){ int t=ec[(long)z*sz+(long)y*nx+x]; if(t>mx) mx=t; }
    int hb= mx<=2?0:(mx==3?1:2);
    if(e) e_bit(e,&M.bf[nb][hb],has); else has=d_bit(d,&M.bf[nb][hb]);
    bf[bi]=has; if(!has) continue;
    for(int z=Z*OB;z<nz&&z<Z*OB+OB;z++) for(int y=Y*OB;y<ny&&y<Y*OB+OB;y++) for(int x=X*OB;x<nx&&x<X*OB+OB;x++){ long i=(long)z*sz+(long)y*nx+x;
      int nf=(x>0&&k[i-1]!=0)+(y>0&&k[i-nx]!=0)+(z>0&&k[i-sz]!=0); if(nf>3)nf=3; int v;
      if(e){ v=k[i]; e_bit(e,&M.f[ec[i]][nf],v!=0); if(!v) continue; e_bit(e,&M.s,v<0); int a=(v<0?-v:v)-1;
        int j=0; for(;j<8;j++){ int b=a>j; e_bit(e,&M.mg[j],b); if(!b) break; }
        if(j==8){ unsigned u=(unsigned)(a-8)+1; int nbits=0; while((u>>nbits)>1) nbits++; for(int t=0;t<nbits;t++) e_bit(e,&M.eg[t],1); e_bit(e,&M.eg[nbits],0); for(int t=nbits-1;t>=0;t--) e_bitp(e,(u>>t)&1,32768); }
      } else { k[i]=0; if(!d_bit(d,&M.f[ec[i]][nf])) continue; int neg=d_bit(d,&M.s); int j=0; for(;j<8;j++) if(!d_bit(d,&M.mg[j])) break;
        int a=j; if(j==8){ int nbits=0; while(d_bit(d,&M.eg[nbits])) nbits++; unsigned u=1; for(int t=0;t<nbits;t++) u=(u<<1)|d_bitp(d,32768); a=8+(int)(u-1); }
        a+=1; k[i]=neg?-a:a; } } }
  free(bf); return 0;
}
static void make_ec(const float* b,int nz,int ny,int nx,float tau,uint8_t* ec){ double t_=cpun(); make_ec_raw(b,nz,ny,nx,tau,ec); PF[5]+=cpun()-t_;  }
static long out_code(int32_t* k,const uint8_t* ec,int nz,int ny,int nx,Enc* e,Dec* d){ double t_=cpun(); long r_=out_code_raw(k,ec,nz,ny,nx,e,d); PF[4]+=cpun()-t_; return r_; }
/* ---------------- full pipeline ---------------- */
static int g_init=0;
static void init_gains(void){ if(g_init) return; int n=2048; float* a=(float*)balloc(sizeof(float)*n);
  for(int kind=0;kind<2;kind++) for(int l=0;l<12;l++){
    if(kind==0&&l==0){ G1[0][0]=1; continue; }
    int levels= kind? l+1 : l; int c=n; for(int t=0;t<l;t++) c=(c+1)/2; int h=(c+1)/2; int pos= kind? (h+c)/2 : c/2;
    memset(a,0,sizeof(float)*n); a[pos]=1; dwt3(a,1,1,n,0,0,levels,1); double s=0; for(int i=0;i<n;i++) s+=(double)a[i]*a[i]; G1[kind][l]=s; }
  free(a); g_init=1; }
static int Lmax3(int a,int b,int c){ int L=a>b?a:b; return L>c?L:c; }
static void box_w(int* bx,int nb,const int* Ls,float* w){ init_gains();
  for(int i=0;i<nb;i++){ int* B=bx+9*i; double g=1; int lev=B[6],o=B[7];
    for(int ax=0;ax<3;ax++){ int bit=4>>ax; if(o==0) g*=G1[0][Ls[ax]]; else if(o&bit) g*=G1[1][lev]; else g*=G1[0][(lev+1<Ls[ax])?lev+1:Ls[ax]]; }
    w[i]=(float)sqrt(g); } }
#define NFC 8
#ifndef FITV
#define FITV 1
#endif
#ifndef NPAIR
#define NPAIR 16  /* r10 default: 16 neighbor pairs; do not change without re-benching */
#endif
#define NFT (NPAIR+3)
typedef struct { int32_t magic,nz,ny,nx,Lz,Ly,Lx,vmin,vmax,tau; float d; float off[16]; int64_t l1,l2; int32_t nfc,pad; int32_t isf,pad2; float fmn,fmx,ftau,fc0; float fth[NFC]; float fw[NFC][NFT]; } Hdr;
/* ---- adaptive Wiener post-filter (ALF-style, 3x3x3 point-symmetric, activity classes), LS-fitted by encoder ---- */
static const int PO[28][3]={{0,0,2},{0,2,0},{2,0,0},{0,0,1},{0,1,0},{1,0,0},{0,1,1},{0,1,-1},{1,0,1},{1,0,-1},{1,1,0},{1,-1,0},{1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
 {0,1,2},{0,1,-2},{0,2,1},{0,2,-1},{1,0,2},{1,0,-2},{2,0,1},{2,0,-1},{1,2,0},{1,-2,0},{2,1,0},{2,-1,0}};
static inline void feats(const float* b,long i,const long* po,float* f){ f[0]=b[i]; for(int k=0;k<NPAIR;k++) f[k+1]=b[i+po[k]]+b[i-po[k]]; f[NPAIR+1]=1.0f; f[NPAIR+2]=0; }
static inline float act(const float* b,long i,long sz,int nx){ float c2=2*b[i]; return fabsf(c2-b[i-1]-b[i+1])+fabsf(c2-b[i-nx]-b[i+nx])+fabsf(c2-b[i-sz]-b[i+sz]); }
static inline int fcls(float a,const float* th,int K){ int k=0; while(k<K-1 && a>=th[k]) k++; return k; }
static int cmpf(const void* a,const void* b){ float x=*(const float*)a,y=*(const float*)b; return x<y?-1:(x>y); }
static void solve(double* A,double* b,int n,float* w){ /* gaussian elimination with ridge */
  double tr=0; for(int i=0;i<n;i++) tr+=A[i*n+i]; for(int i=0;i<n;i++) A[i*n+i]+=1e-7*tr/n+1e-9;
  for(int c=0;c<n;c++){ int p=c; for(int r=c+1;r<n;r++) if(fabs(A[r*n+c])>fabs(A[p*n+c])) p=r;
    if(p!=c){ for(int k=0;k<n;k++){ double t=A[c*n+k];A[c*n+k]=A[p*n+k];A[p*n+k]=t;} double t=b[c];b[c]=b[p];b[p]=t; }
    double dv=A[c*n+c]; if(fabs(dv)<1e-30) dv=1e-30;
    for(int r=c+1;r<n;r++){ double f=A[r*n+c]/dv; if(f==0) continue; for(int k=c;k<n;k++) A[r*n+k]-=f*A[c*n+k]; b[r]-=f*b[c]; } }
  for(int c=n-1;c>=0;c--){ double s=b[c]; for(int k=c+1;k<n;k++) s-=A[c*n+k]*w[k]; double dv=A[c*n+c]; if(fabs(dv)<1e-30) dv=1e-30; w[c]=(float)(s/dv); } }
static void fit_filter_raw(const float* b,const uint16_t* v,const float* vf,float c0,int nz,int ny,int nx,Hdr* h,int K){
  long sz=(long)ny*nx; long po[28]; for(int k=0;k<NPAIR;k++) po[k]=PO[k][0]*sz+PO[k][1]*(long)nx+PO[k][2];
  h->nfc=0; if(nz<5||ny<5||nx<5) return;
  long ni=(long)(nz-4)*(ny-4)*(nx-4); int stp=1; while(ni/(stp*stp)>500000) stp++;
  /* activity thresholds from sample quantiles */
  long ns=0, cap=200000; float* sa=(float*)malloc(sizeof(float)*cap); int st2=1; while(ni/st2>cap) st2*=2;
  long cnt=0; for(int z=2;z<nz-2&&ns<cap;z++) for(int y=2;y<ny-2;y++){ int w_=nx-4; /* same picks as (cnt++ % st2)==0, without a division per voxel */
    long r0=(st2-cnt%st2)%st2; for(long x=2+r0;x<nx-2 && ns<cap;x+=st2) sa[ns++]=act(b,(long)z*sz+(long)y*nx+x,sz,nx); cnt+=w_; }
  qsort(sa,ns,sizeof(float),cmpf); for(int k=0;k<K-1;k++) h->fth[k]=sa[(long)(ns*(k+1))/K]; free(sa);
  double* A=(double*)calloc((size_t)K*NFT*NFT,sizeof(double)); double* bb=(double*)calloc((size_t)K*NFT,sizeof(double)); double* e0=(double*)calloc(K,sizeof(double)); float f[NFT];
#if FITV && WLC_AVX2
  /* AVX2 accumulation: float*float products are exact in double, and every element keeps its sample order -> bit-identical to the scalar loop */
  { enum { NP=(NFT+3)&~3 };
#ifdef _MSC_VER
    double* A2=(double*)_aligned_malloc(sizeof(double)*K*NFT*NP+32,32); double* b2=(double*)_aligned_malloc(sizeof(double)*K*NP+32,32);
    __declspec(align(32)) double fd[NP];
#else
    double* A2=(double*)aligned_alloc(32,sizeof(double)*K*NFT*NP+32); double* b2=(double*)aligned_alloc(32,sizeof(double)*K*NP+32);
    double fd[NP] __attribute__((aligned(32)));
#endif
    memset(A2,0,sizeof(double)*K*NFT*NP); memset(b2,0,sizeof(double)*K*NP);
    for(int z=2;z<nz-2;z++) for(int y=2;y<ny-2;y+=stp) for(int x=2+(y/stp+z)%stp;x<nx-2;x+=stp){ long i=(long)z*sz+(long)y*nx+x;
      int k=fcls(act(b,i,sz,nx),h->fth,K); feats(b,i,po,f); double t= vf? (double)vf[i]-c0 : (double)v[i]-c0;
      for(int a=0;a<NP;a++) fd[a]= a<NFT? (double)f[a] : 0.0;
      double* Ak=A2+(size_t)k*NFT*NP; double* bk=b2+(size_t)k*NP; __m256d T=_mm256_set1_pd(t);
      for(int c=0;c<NP;c+=4) _mm256_store_pd(bk+c,_mm256_add_pd(_mm256_load_pd(bk+c),_mm256_mul_pd(_mm256_load_pd(fd+c),T)));
      for(int a=0;a<NFT;a++){ __m256d FA=_mm256_set1_pd(fd[a]); double* row=Ak+a*NP;
        for(int c=a&~3;c<NP;c+=4) _mm256_store_pd(row+c,_mm256_add_pd(_mm256_load_pd(row+c),_mm256_mul_pd(FA,_mm256_load_pd(fd+c)))); }
      e0[k]+=(t-f[0])*(t-f[0]); }
    for(int k=0;k<K;k++){ for(int a=0;a<NFT;a++){ bb[k*NFT+a]=b2[k*NP+a]; for(int c2=a;c2<NFT;c2++) A[(size_t)k*NFT*NFT+a*NFT+c2]=A2[(size_t)k*NFT*NP+a*NP+c2]; } }
#ifdef _MSC_VER
    _aligned_free(A2); _aligned_free(b2);
#else
    free(A2); free(b2);
#endif
  }

#else
  for(int z=2;z<nz-2;z++) for(int y=2;y<ny-2;y+=stp) for(int x=2+(y/stp+z)%stp;x<nx-2;x+=stp){ long i=(long)z*sz+(long)y*nx+x;
    int k=fcls(act(b,i,sz,nx),h->fth,K); feats(b,i,po,f); double t= vf? (double)vf[i]-c0 : (double)v[i]-c0; double* Ak=A+(size_t)k*NFT*NFT;
    for(int a=0;a<NFT;a++){ double fa=f[a]; bb[k*NFT+a]+=fa*t; for(int c2=a;c2<NFT;c2++) Ak[a*NFT+c2]+=fa*f[c2]; } e0[k]+=(t-f[0])*(t-f[0]); }
#endif
  double gain=0, base=0;
  for(int k=0;k<K;k++){ double* Ak=A+(size_t)k*NFT*NFT; for(int a=0;a<NFT;a++) for(int c2=0;c2<a;c2++) Ak[a*NFT+c2]=Ak[c2*NFT+a];
    double Ac[NFT*NFT], bc[NFT]; memcpy(Ac,Ak,sizeof(Ac)); memcpy(bc,bb+k*NFT,sizeof(bc)); float w[NFT]; solve(Ac,bc,NFT,w);
    /* residual energy reduction = 2 w.b - w'Aw  vs identity filter */
    double wb=0,wAw=0,ib=bb[k*NFT],iAi=Ak[0]; for(int a=0;a<NFT;a++){ wb+=w[a]*bb[k*NFT+a]; for(int c2=0;c2<NFT;c2++) wAw+=w[a]*Ak[a*NFT+c2]*w[c2]; }
    double red=(2*wb-wAw)-(2*ib-iAi); gain+=red; base+=e0[k]; memcpy(h->fw[k],w,sizeof(w)); }
  free(A); free(bb); free(e0);
  if(gain>0.005*base) h->nfc=K;
}
static void fit_filter(const float* b,const uint16_t* v,const float* vf,float c0,int nz,int ny,int nx,Hdr* h,int K){ double t_=cpun(); fit_filter_raw(b,v,vf,c0,nz,ny,nx,h,K); PF[2]+=cpun()-t_;  }
static void apply_filter(const float* b,float* o,int nz,int ny,int nx,const Hdr* h);
#if WLC_AVX2
#ifndef NOFMA
#define FMAV(a,b,c) _mm256_fmadd_ps(a,b,c)
#else
#define FMAV(a,b,c) _mm256_add_ps(c,_mm256_mul_ps(a,b))
#endif
#endif
static void apply_filter_raw(const float* b,float* o,int nz,int ny,int nx,const Hdr* h){
  long sz=(long)ny*nx; long n=(long)nz*sz; long po[28]; for(int k=0;k<NPAIR;k++) po[k]=PO[k][0]*sz+PO[k][1]*(long)nx+PO[k][2];
  memcpy(o,b,sizeof(float)*n); int K=h->nfc; int32_t* cl=(int32_t*)malloc(sizeof(int32_t)*(nx+8));
  float wt[NFT][8]; for(int a=0;a<NFT;a++) for(int k=0;k<8;k++) wt[a][k]= k<K? h->fw[k][a] : 0;
#if WLC_AVX2
  float* acc=(float*)malloc(sizeof(float)*(nx+8));
  __m256 W[NFT]; for(int a=0;a<NFT;a++) W[a]=_mm256_loadu_ps(wt[a]);
  for(int z=2;z<nz-2;z++) for(int y=2;y<ny-2;y++){ long i0=(long)z*sz+(long)y*nx; const float* p=b+i0;
    { int x=2; __m256 AB=_mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff)), TWO=_mm256_set1_ps(2.f);
      for(;x+8<=nx-2;x+=8){ __m256 c2=_mm256_mul_ps(TWO,_mm256_loadu_ps(p+x));
        __m256 a1=_mm256_and_ps(_mm256_sub_ps(_mm256_sub_ps(c2,_mm256_loadu_ps(p+x-1)),_mm256_loadu_ps(p+x+1)),AB);
        __m256 a2=_mm256_and_ps(_mm256_sub_ps(_mm256_sub_ps(c2,_mm256_loadu_ps(p+x-nx)),_mm256_loadu_ps(p+x+nx)),AB);
        __m256 a3=_mm256_and_ps(_mm256_sub_ps(_mm256_sub_ps(c2,_mm256_loadu_ps(p+x-sz)),_mm256_loadu_ps(p+x+sz)),AB);
        __m256 ac=_mm256_add_ps(_mm256_add_ps(a1,a2),a3); __m256i k=_mm256_setzero_si256();
        for(int t=0;t<K-1;t++) k=_mm256_sub_epi32(k,_mm256_castps_si256(_mm256_cmp_ps(ac,_mm256_set1_ps(h->fth[t]),_CMP_GE_OQ)));
        _mm256_storeu_si256((__m256i*)(cl+x),k); }
      for(;x<nx-2;x++){ float c2=2*p[x]; float ac=fabsf(c2-p[x-1]-p[x+1])+fabsf(c2-p[x-nx]-p[x+nx])+fabsf(c2-p[x-sz]-p[x+sz]); int k=0; for(int t=0;t<K-1;t++) k+=(ac>=h->fth[t]); cl[x]=k; } }
    int x=2;
    for(;x+8<=nx-2;x+=8){ __m256i id=_mm256_loadu_si256((const __m256i*)(cl+x));
      __m256 s=_mm256_add_ps(_mm256_mul_ps(_mm256_permutevar8x32_ps(W[0],id),_mm256_loadu_ps(p+x)),_mm256_permutevar8x32_ps(W[NPAIR+1],id));
      for(int a=0;a<NPAIR;a++){ long oo=po[a]; __m256 f=_mm256_add_ps(_mm256_loadu_ps(p+x+oo),_mm256_loadu_ps(p+x-oo)); s=FMAV(_mm256_permutevar8x32_ps(W[a+1],id),f,s); }
      _mm256_storeu_ps(o+i0+x,s); }
    for(;x<nx-2;x++){ int k=cl[x]; float s=wt[0][k]*p[x]+wt[NPAIR+1][k]; for(int a=0;a<NPAIR;a++){ long oo=po[a]; s+=wt[a+1][k]*(p[x+oo]+p[x-oo]); } o[i0+x]=s; }
  }
  free(cl); free(acc);
#else
  for(int z=2;z<nz-2;z++) for(int y=2;y<ny-2;y++){ long i0=(long)z*sz+(long)y*nx; const float* p=b+i0;
    int x=2;
    for(;x+8<=nx-2;x+=8){
      for(int k=0;k<8;k++){
        float c2=2*p[x+k];
        float a1=fabsf((c2-p[x+k-1])-p[x+k+1]);
        float a2=fabsf((c2-p[x+k-nx])-p[x+k+nx]);
        float a3=fabsf((c2-p[x+k-sz])-p[x+k+sz]);
        float ac=(a1+a2)+a3; int cls=0; for(int t=0;t<K-1;t++) cls+=(ac>=h->fth[t]); cl[x+k]=cls;
      }
    }
    for(;x<nx-2;x++){ float c2=2*p[x]; float ac=fabsf(c2-p[x-1]-p[x+1])+fabsf(c2-p[x-nx]-p[x+nx])+fabsf(c2-p[x-sz]-p[x+sz]); int k=0; for(int t=0;t<K-1;t++) k+=(ac>=h->fth[t]); cl[x]=k; }
    x=2;
    for(;x+8<=nx-2;x+=8){
      for(int k=0;k<8;k++){
        int cls=cl[x+k]; float s=wt[0][cls]*p[x+k]+wt[NPAIR+1][cls];
        for(int a=0;a<NPAIR;a++){ long oo=po[a]; float f=p[x+k+oo]+p[x+k-oo];
#ifndef NOFMA
          s=fmaf(wt[a+1][cls],f,s);
#else
          s+=wt[a+1][cls]*f;
#endif
        }
        o[i0+x+k]=s;
      }
    }
    for(;x<nx-2;x++){ int k=cl[x]; float s=wt[0][k]*p[x]+wt[NPAIR+1][k]; for(int a=0;a<NPAIR;a++){ long oo=po[a]; s+=wt[a+1][k]*(p[x+oo]+p[x-oo]); } o[i0+x]=s; }
  }
  free(cl);
#endif
}
static void recon2(const int32_t* q,float* buf,int32_t* r,int nz,int ny,int nx,const int* Ls,int* bx,int nb,const float* w,float d,const float* off,int vmin,int vmax,Hdr* h,const uint16_t* fitv,const float* fitf,int K,float c);
static void recon(const int32_t* q,float* buf,int32_t* r,int nz,int ny,int nx,const int* Ls,int* bx,int nb,const float* w,float d,const float* off,int vmin,int vmax,Hdr* h,const uint16_t* fitv,int K){ recon2(q,buf,r,nz,ny,nx,Ls,bx,nb,w,d,off,vmin,vmax,h,fitv,NULL,K,0.5f*(vmin+vmax)); }
static float* RSCR=NULL; /* caller-provided scratch (n floats) for the filter output; avoids a fresh n-float allocation (page faults = sys CPU) */
static void recon2(const int32_t* q,float* buf,int32_t* r,int nz,int ny,int nx,const int* Ls,int* bx,int nb,const float* w,float d,const float* off,int vmin,int vmax,Hdr* h,const uint16_t* fitv,const float* fitf,int K,float c){
  long sz=(long)ny*nx; int L=Lmax3(Ls[0],Ls[1],Ls[2]);
  for(int k=0;k<nb;k++){ int* B=bx+9*k; float o=off[clsof(B[6],B[7],L)]; float s=d/w[k];
    for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ long i0=(long)z*sz+(long)y*nx;
      int x=B[4];
#if WLC_AVX2
      { const __m256i Z8=_mm256_setzero_si256(); const __m256 O8=_mm256_set1_ps(o), S8=_mm256_set1_ps(s), ONE8=_mm256_set1_ps(1.0f);
        for(;x+8<=B[5];x+=8){ __m256i v8=_mm256_loadu_si256((const __m256i*)(q+i0+x)); __m256 sg=_mm256_sub_ps(_mm256_and_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(v8,Z8)),ONE8),_mm256_and_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(Z8,v8)),ONE8));
          _mm256_storeu_ps(buf+i0+x,_mm256_mul_ps(_mm256_add_ps(_mm256_cvtepi32_ps(v8),_mm256_mul_ps(sg,O8)),S8)); } } /* r10: AVX dequant, same float ops */
#endif
      for(;x<B[5];x++){ int v=q[i0+x]; float sg=(float)((v>0)-(v<0)); buf[i0+x]=((float)v+sg*o)*s; } } }
  dwt3(buf,nz,ny,nx,Ls[0],Ls[1],Ls[2],1);
  long n=(long)nz*sz;
  if(fitv||fitf){ if(K>1) fit_filter(buf,fitv,fitf,c,nz,ny,nx,h,K); else h->nfc=0; }
  if(h->nfc>0){ float* o= RSCR? RSCR : (float*)balloc(sizeof(float)*n); apply_filter(buf,o,nz,ny,nx,h); memcpy(buf,o,sizeof(float)*n); if(!RSCR) free(o); }
  if(r){ long i=0;
#if WLC_AVX2
    __m256 C=_mm256_set1_ps(c); __m256i LO=_mm256_set1_epi32(vmin), HI=_mm256_set1_epi32(vmax);
    for(;i+8<=n;i+=8){ __m256i t=_mm256_cvtps_epi32(_mm256_add_ps(_mm256_loadu_ps(buf+i),C)); t=_mm256_max_epi32(t,LO); t=_mm256_min_epi32(t,HI); _mm256_storeu_si256((__m256i*)(r+i),t); }
#endif
    for(;i<n;i++){ long t=lrintf(buf[i]+c); if(t<vmin)t=vmin; if(t>vmax)t=vmax; r[i]=(int32_t)t; } }
}
static void apply_filter(const float* b,float* o,int nz,int ny,int nx,const Hdr* h){ double t_=cpun(); apply_filter_raw(b,o,nz,ny,nx,h); PF[3]+=cpun()-t_;  }
/* ---- rate control helpers: predicted MSE of the deadzone quantizer with per-class centroids on a coefficient sample ---- */
#if WLC_AVX2
static inline float hsum8(__m256 v){ float t[8]; _mm256_storeu_ps(t,v); return hsum8s(t); }
#endif
/* samples are kept sorted by class (csort) so every class is a contiguous run: no gathers, AVX loops */
static void csort(uint8_t* sc,float* sv,long ns){ long cnt[17]={0}; for(long j=0;j<ns;j++) cnt[sc[j]+1]++; for(int c=0;c<16;c++) cnt[c+1]+=cnt[c];
  float* t=(float*)malloc(sizeof(float)*ns); long p[16]; for(int c=0;c<16;c++) p[c]=cnt[c]; for(long j=0;j<ns;j++) t[p[sc[j]]++]=sv[j];
  memcpy(sv,t,sizeof(float)*ns); free(t); for(int c=0;c<16;c++) memset(sc+cnt[c],c,cnt[c+1]-cnt[c]); }
static void csegs(const uint8_t* sc,long ns,long* sg){ for(int c=0;c<=16;c++){ long lo=0,hi=ns; while(lo<hi){ long m=(lo+hi)/2; if(sc[m]<c) lo=m+1; else hi=m; } sg[c]=lo; } }
static double rc_pred(const uint8_t* sc,const float* sv,long ns,double dd,double rnd,double addc,float* off){
  long sg[17]; csegs(sc,ns,sg); float id_=(float)(1.0/dd), fd_=(float)dd; double tot=0;
#if WLC_AVX2
  __m256 ID=_mm256_set1_ps(id_), RN=_mm256_set1_ps((float)rnd), Z=_mm256_setzero_ps(), ONE=_mm256_set1_ps(1.f), FD=_mm256_set1_ps(fd_);
#endif
  for(int c=0;c<16;c++){ long a=sg[c], b=sg[c+1]; if(b<=a){ off[c]=(float)(0.5-rnd); continue; }
    double sm=0, ct=0; long j=a;
#if WLC_AVX2
    while(j+8<=b){ __m256 S=Z, C=Z; long e=j+2048<b? j+2048 : b; for(;j+8<=e;j+=8){ __m256 aa=_mm256_mul_ps(_mm256_loadu_ps(sv+j),ID); __m256 qq=_mm256_floor_ps(_mm256_add_ps(aa,RN)); __m256 nz=_mm256_cmp_ps(qq,Z,_CMP_GT_OQ);
        S=_mm256_add_ps(S,_mm256_and_ps(nz,_mm256_sub_ps(aa,qq))); C=_mm256_add_ps(C,_mm256_and_ps(nz,ONE)); } sm+=hsum8(S); ct+=hsum8(C); }
#else
    while(j+8<=b){ float S[8]={0}, C[8]={0}; long e=j+2048<b? j+2048 : b; for(;j+8<=e;j+=8){ for(int k=0;k<8;k++){ float aa=sv[j+k]*id_; float qq=floorf(aa+(float)rnd); if(qq>0){ S[k]+=aa-qq; C[k]+=1.f; } } } sm+=hsum8s(S); ct+=hsum8s(C); }
#endif
    for(;j<b;j++){ float aa=sv[j]*id_; float qq=floorf(aa+(float)rnd); if(qq>0){ sm+=aa-qq; ct+=1; } }
    double o= ct>0? sm/ct : 0.5-rnd; if(o<-0.5)o=-0.5; if(o>0.5)o=0.5; off[c]=(float)o;
    j=a;
#if WLC_AVX2
    { __m256 OF=_mm256_set1_ps(off[c]);
    while(j+8<=b){ __m256 S=Z; long e=j+2048<b? j+2048 : b; for(;j+8<=e;j+=8){ __m256 x=_mm256_loadu_ps(sv+j); __m256 qq=_mm256_floor_ps(_mm256_add_ps(_mm256_mul_ps(x,ID),RN)); __m256 nz=_mm256_cmp_ps(qq,Z,_CMP_GT_OQ);
        __m256 rr=_mm256_and_ps(nz,_mm256_mul_ps(_mm256_add_ps(qq,OF),FD)); __m256 d=_mm256_sub_ps(x,rr); S=_mm256_add_ps(S,_mm256_mul_ps(d,d)); } tot+=hsum8(S); } }
#else
    while(j+8<=b){ float S[8]={0}; long e=j+2048<b? j+2048 : b; for(;j+8<=e;j+=8){ for(int k=0;k<8;k++){ float x=sv[j+k]; float qq=floorf(x*id_+(float)rnd); float rr= qq>0? (qq+off[c])*fd_:0; float d=x-rr; S[k]+=d*d; } } tot+=hsum8s(S); }
#endif
    for(;j<b;j++){ float aa=sv[j]; float qq=floorf(aa*id_+(float)rnd); float rr= qq>0? (qq+off[c])*fd_:0; tot+=(double)(aa-rr)*(aa-rr); } }
  return tot/ns+addc; }
static double rc_solve(const uint8_t* sc,const float* sv,long ns,double rnd,double addc,double tm,double lo,double hi,float* off){
  for(int it=0;it<80 && hi/lo>1.002;it++){ double mid=sqrt(lo*hi); double pm=rc_pred(sc,sv,ns,mid,rnd,addc,off); if(pm>tm) hi=mid; else lo=mid; } return lo; }
static double rho_at(const float* sv,long ns,double d,double rnd){ long c_=0; float id2=(float)(1.0/d); for(long j=0;j<ns;j++) if(floorf(sv[j]*id2+(float)rnd)>0) c_++; return (double)c_/ns; }
/* L-inf mode: choose kt = tau/rmse from the significance density rho of the quantised coefficients (fixed point of
   kt = g(rho(kt)), g = 0.65 - 0.79 log2 rho, clipped to [1.3,8]): dense spectra -> near-uniform error, small kt; sparse -> large kt. */
static double auto_kt(const uint8_t* sc,const float* sv,long ns,double rnd,double addc,double linf,double lo,double hi,double ga){
  float off[16]; long st=ns/150000; if(st<1)st=1; long ns2=(ns+st-1)/st; uint8_t* sc2=(uint8_t*)malloc(ns2); float* sv2=(float*)malloc(sizeof(float)*ns2);
  for(long j=0,i=0;i<ns;i+=st,j++){ sc2[j]=sc[i]; sv2[j]=sv[i]; } double klo=1.3, khi=8;
  double dlo=lo, dhi=hi; /* d(kt) is decreasing in kt: each step narrows the step-size bracket of the next solve */
  for(int it=0;it<9;it++){ double km=sqrt(klo*khi); double tm=(linf/km)*(linf/km); double d=rc_solve(sc2,sv2,ns2,rnd,addc,tm,dlo,dhi,off);
    double rho=rho_at(sv2,ns2,d,rnd); double g= rho>0? ga-0.79*log2(rho) : 8; if(g<1.3)g=1.3; if(g>8)g=8; if(km<g){ klo=km; if(d*1.01<dhi) dhi=d*1.01; } else { khi=km; if(d/1.01>dlo) dlo=d/1.01; } }
  free(sc2); free(sv2); return sqrt(klo*khi); }
/* decorrelation test: max over axes of lag-1 Pearson correlation on ~200k sampled pairs. Near-white data (e.g. NN weights)
   gains nothing from the wavelet; it is then coded in 'direct' mode: no transform, uniform step 2*tau, context-adaptive coding. */
static double lag1corr_f(const float* v,int nz,int ny,int nx){ long sz=(long)ny*nx, n=(long)nz*sz; long strides[3]={sz,nx,1}; int dims[3]={nz,ny,nx}; double best=-1;
  for(int a=0;a<3;a++){ if(dims[a]<2) continue; long s=strides[a]; long m=n-s; if(m<=0) continue; long np=m<200000? m:200000; long P=1000003L%m; if(P==0)P=1;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0; long i=0; for(long j=0;j<np;j++){ double x=v[i], y=v[i+s]; sx+=x; sy+=y; sxx+=x*x; syy+=y*y; sxy+=x*y; i+=P; if(i>=m) i-=m; }
    double cx=sxx/np-(sx/np)*(sx/np), cy=syy/np-(sy/np)*(sy/np), c=sxy/np-(sx/np)*(sy/np); double r= (cx>0&&cy>0)? c/sqrt(cx*cy) : 1; if(r>best) best=r; }
  return best; }
static void levels_for(const int* sh,int* Ls,int maxl){ for(int a=0;a<3;a++){ int ml= a==0? (maxl/10) : (maxl%10); if(maxl<10) ml=maxl; if(maxl>=100) ml= a==0? 0 : maxl%10; /* 1xy: no z transform; 2pxy: z packet with p levels (PKZ) */ int L=0; double n=sh[a]; while(L<ml && n/4.0>=2.0){ n=ceil(n/2); L++; } Ls[a]=L; } }
int FK=6;
WLC_API void wlc_setfk(int k){ FK=k; }
static void fminmax(const float* v,long n,float* mn,float* mx){
  long i=0; float lo,hi;
#if WLC_AVX2
  __m256 a=_mm256_set1_ps(v[0]), b=a;
  for(;i+8<=n;i+=8){ __m256 x=_mm256_loadu_ps(v+i); a=_mm256_min_ps(a,x); b=_mm256_max_ps(b,x); }
  float ta[8],tb[8]; _mm256_storeu_ps(ta,a); _mm256_storeu_ps(tb,b); lo=ta[0]; hi=tb[0]; for(int k=1;k<8;k++){ if(ta[k]<lo)lo=ta[k]; if(tb[k]>hi)hi=tb[k]; }
#else
  float a[8], b[8]; for(int k=0;k<8;k++) a[k]=b[k]=v[0];
  for(;i+8<=n;i+=8){ for(int k=0;k<8;k++){ float x=v[i+k]; a[k]= (a[k]<x)? a[k]:x; b[k]= (b[k]>x)? b[k]:x; } }
  lo=a[0]; hi=b[0]; for(int k=1;k<8;k++){ if(a[k]<lo)lo=a[k]; if(b[k]>hi)hi=b[k]; }
#endif
  for(;i<n;i++){ if(v[i]<lo)lo=v[i]; if(v[i]>hi)hi=v[i]; } *mn=lo; *mx=hi; }
/* quantise a row: cn<-cn*id (scaled), q=sign*floor(|cn|+rnd); returns sum of (|cn|-|q|) over nonzero q and their count */
static inline void qrow(float* restrict cn,int32_t* restrict q,int w,float id,float rnd,double* sm,long* ct){ float s=0; int c=0; int x=0;
#if WLC_AVX2
  __m256 ID=_mm256_set1_ps(id), RN=_mm256_set1_ps(rnd), AB=_mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff)), SG=_mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)), Z=_mm256_setzero_ps(), ONE=_mm256_set1_ps(1.f);
  __m256 S=Z, CC=Z;
  for(;x+8<=w;x+=8){ __m256 a=_mm256_mul_ps(_mm256_loadu_ps(cn+x),ID); _mm256_storeu_ps(cn+x,a); __m256 aa=_mm256_and_ps(a,AB);
    __m256 qq=_mm256_floor_ps(_mm256_add_ps(aa,RN)); __m256 nz=_mm256_cmp_ps(qq,Z,_CMP_GT_OQ);
    S=_mm256_add_ps(S,_mm256_and_ps(nz,_mm256_sub_ps(aa,qq))); CC=_mm256_add_ps(CC,_mm256_and_ps(nz,ONE));
    __m256 sq=_mm256_or_ps(qq,_mm256_and_ps(a,SG)); _mm256_storeu_si256((__m256i*)(q+x),_mm256_cvttps_epi32(sq)); }
  s=hsum8(S); c=(int)hsum8(CC);
#else
  float S[8]={0}, CC[8]={0};
  for(;x+8<=w;x+=8){ for(int k=0;k<8;k++){ float a=cn[x+k]*id; cn[x+k]=a; float aa=fabsf(a); float qq=floorf(aa+rnd); if(qq>0){ S[k]+=aa-qq; CC[k]+=1.f; } q[x+k]=(int32_t)(a<0?-qq:qq); } }
  s=hsum8s(S); c=(int)hsum8s(CC);
#endif
  for(;x<w;x++){ float a=cn[x]*id; cn[x]=a; float aa=fabsf(a); float qq=floorf(aa+rnd); float nzm=(qq>0)?1.f:0.f; s+=nzm*(aa-qq); c+=(int)nzm; q[x]=(int32_t)(a<0?-qq:qq); }
  *sm+=s; *ct+=c; }
static inline void orow(const float* restrict cn,const int32_t* restrict q,int w,double* sm,long* ct){ float s=0; int c=0; int x=0;
#if WLC_AVX2
  __m256 AB=_mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff)), Z=_mm256_setzero_ps(), ONE=_mm256_set1_ps(1.f), S=Z, CC=Z;
  for(;x+8<=w;x+=8){ __m256 aq=_mm256_cvtepi32_ps(_mm256_abs_epi32(_mm256_loadu_si256((const __m256i*)(q+x)))); __m256 nz=_mm256_cmp_ps(aq,Z,_CMP_GT_OQ);
    S=_mm256_add_ps(S,_mm256_and_ps(nz,_mm256_sub_ps(_mm256_and_ps(_mm256_loadu_ps(cn+x),AB),aq))); CC=_mm256_add_ps(CC,_mm256_and_ps(nz,ONE)); }
  s=hsum8(S); c=(int)hsum8(CC);
#else
  float S[8]={0}, CC[8]={0};
  for(;x+8<=w;x+=8){ for(int k=0;k<8;k++){ int qq=q[x+k]; int aq=qq<0?-qq:qq; if(aq){ S[k]+=fabsf(cn[x+k])-(float)aq; CC[k]+=1.f; } } }
  s=hsum8s(S); c=(int)hsum8s(CC);
#endif
  for(;x<w;x++){ int qq=q[x]; int aq=qq<0?-qq:qq; float nzm= aq? 1.f:0.f; s+=nzm*(fabsf(cn[x])-(float)aq); c+=aq?1:0; } *sm+=s; *ct+=c; }
/* clamp recon to [mn,mx], return sum sq err (double) and max |err| (float) */
static void finrow(const float* restrict v,float* restrict b,long n,float c0,float mn,float mx,double* ss,float* me){
  double tot=0; long i=0; float m=0;
#if WLC_AVX2
  __m256 C=_mm256_set1_ps(c0), MN=_mm256_set1_ps(mn), MX=_mm256_set1_ps(mx), M=_mm256_setzero_ps(), AB=_mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
  while(i+8<=n){ __m256 S=_mm256_setzero_ps(); long e=i+512<n? i+512: n; for(;i+8<=e;i+=8){ __m256 xr=_mm256_add_ps(_mm256_loadu_ps(b+i),C); xr=_mm256_min_ps(_mm256_max_ps(xr,MN),MX); _mm256_storeu_ps(b+i,xr);
      __m256 d=_mm256_sub_ps(_mm256_loadu_ps(v+i),xr); S=_mm256_add_ps(S,_mm256_mul_ps(d,d)); M=_mm256_max_ps(M,_mm256_and_ps(d,AB)); }
    float t[8]; _mm256_storeu_ps(t,S); for(int k=0;k<8;k++) tot+=t[k]; }
  { float t[8]; _mm256_storeu_ps(t,M); for(int k=0;k<8;k++) if(t[k]>m)m=t[k]; }
#else
  float M[8]={0};
  while(i+8<=n){ float S[8]={0}; long e=i+512<n? i+512: n; for(;i+8<=e;i+=8){ for(int k=0;k<8;k++){ float xr=b[i+k]+c0; if(xr<mn)xr=mn; if(xr>mx)xr=mx; b[i+k]=xr; float d=v[i+k]-xr; S[k]+=d*d; float ad=fabsf(d); if(ad>M[k]) M[k]=ad; } } for(int k=0;k<8;k++) tot+=S[k]; }
  for(int k=0;k<8;k++) if(M[k]>m)m=M[k];
#endif
  for(;i<n;i++){ float xr=b[i]+c0; if(xr<mn)xr=mn; if(xr>mx)xr=mx; b[i]=xr; float d=v[i]-xr; tot+=(double)d*d; if(fabsf(d)>m)m=fabsf(d); }
  *ss=tot; *me=m; }
/* ---- R8 optional high-effort mode (encoder only; decoder and bitstream format unchanged): post-filter-aware requantisation.
   One Jacobi step with the true post-filter residual: t = chat + a * w * Wfwd(x - F(A chat)) / d, requantised with the same deadzone.
   Pre-compensates what the decoder-side Wiener filter will do. Cost: one extra inverse DWT + filter fit/apply + forward DWT. ---- */
static inline float chat_of(int q,float o,float d){ return q? ((float)q+(q>0?o:-o))*d : 0.f; }
static void he_step(float* cn,int32_t* q,const uint16_t* vu,const float* vf,float c0,int nz,int ny,int nx,const int* Ls,int* bx,int nb,const float* w,int L,float d,const float* off,float rnd,float alpha,int lmax){
  long sz=(long)ny*nx, n=(long)nz*sz; float* buf=(float*)balloc(sizeof(float)*n); Hdr h; memset(&h,0,sizeof(h));
  recon2(q,buf,NULL,nz,ny,nx,Ls,bx,nb,w,d,off,0,0,&h,vu,vf,FK,c0);
  for(long i=0;i<n;i++) buf[i]= (vf? vf[i] : (float)vu[i]) - c0 - buf[i];
  dwt3_raw(buf,nz,ny,nx,Ls[0],Ls[1],Ls[2],0);
  for(int k=0;k<nb;k++){ const int* B=bx+9*k; if(B[7]==0 || B[6]>lmax) continue; float of=off[clsof(B[6],B[7],L)]; float s=w[k]/d;
    for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ long i0=(long)z*sz+(long)y*nx; for(int x=B[4];x<B[5];x++){ long i=i0+x;
      float t=chat_of(q[i],of,1.0f)+alpha*buf[i]*s; float at_=fabsf(t); int m=(int)floorf(at_+rnd); q[i]= t<0?-m:m; cn[i]=t; } } }
  free(buf); }
#ifndef HEDEF
#define HEDEF 0
#endif
static int HE=-1; static int he_on(void){ if(HE<0){ const char* e=getenv("WLC_EFFORT"); HE= e? atoi(e) : HEDEF; } return HE>=2? HE-1 : 0; }

/* ---- trial-chosen z decomposition depth (anisotropy fix). Encoder-only: Lz is already in the header, decoder unchanged.
   Codes a centred crop (all z planes, ~ny/4 x nx/4) at the same absolute quality with Lz = default / 2 / 1 and keeps the
   smallest (needs > ZTM relative gain to leave the default). The PSNR target is re-referenced to the crop's range. ---- */
#ifndef ZTRIAL
#define ZTRIAL 1
#endif
#ifndef ZTM
#define ZTM 0.02
#endif
#ifndef ZT2
#define ZT2 1
#endif
#ifndef ZTINT
#define ZTINT 0
#endif
#ifndef ZTLV
#define ZTLV 1
#endif
#ifndef ZTG
#define ZTG 0.9
#endif
#ifndef ZTC
#define ZTC 4
#endif
#ifndef ZTSKIP
#define ZTSKIP 0 /* r9: dropped - 'no z' does win on some z-smooth hurricane log fields (decision check) */
#endif
static int ZTON=1; static int FDIR;
#ifndef ZTE
#define ZTE 1
#endif
#ifndef ZEHALF
#define ZEHALF 1
#endif
static int pkz_of(int nz,int maxl){ if(maxl<200) return 0; int p=(maxl/10)%10, L=0; double n=nz; while(L<p && n/4.0>=2.0){ n=ceil(n/2); L++; } return L; }
static float* ZCACHE=NULL; static int ZCM=-1; static float* ZPRE=NULL; static float* ZPOOL=NULL; static long ZPOOLN=0;
static double zest_bits(const void* v,int isf,int nz,int ny,int nx,int maxl,double psnr,double linf,double kt,double rnd,double R,float** keep){
  long n=(long)nz*ny*nx, sz=(long)ny*nx; int sh[3]={nz,ny,nx}; int Ls[3]; levels_for(sh,Ls,maxl); int pk0=PKZ; PKZ=pkz_of(nz,maxl); int L=Lmax3(Ls[0],Ls[1],Ls[2]);
  float* cn; float c0; double t0_=cpun(), t1_;
  if(ZPRE){ cn=ZPRE; ZPRE=NULL; t1_=cpun(); goto have_cn; }
  if(ZPOOL && ZPOOLN==n){ cn=ZPOOL; ZPOOL=NULL; } else cn=(float*)balloc(sizeof(float)*n);
  if(isf){ const float* vf=(const float*)v; float a_,b_; fminmax(vf,n,&a_,&b_); c0=0.5f*(a_+b_); for(long i=0;i<n;i++) cn[i]=vf[i]-c0; } /* same c0 and values as the encoder -> winner's transform is reused */
  else { double mn=1e300,mx=-1e300; for(long i=0;i<n;i++){ double t=((const uint16_t*)v)[i]; if(t<mn)mn=t; if(t>mx)mx=t; }
    c0=(float)(0.5*(mn+mx)); for(long i=0;i<n;i++) cn[i]=(float)((const uint16_t*)v)[i]-c0; }
  t0_=cpun(); dwt3_raw(cn,nz,ny,nx,Ls[0],Ls[1],Ls[2],0); t1_=cpun();
  have_cn:;
  int bx[256*9]; int nb=boxes(nz,ny,nx,Ls[0],Ls[1],Ls[2],bx); float w[256]; box_w(bx,nb,Ls,w);
  long step=n/600000; if(step<1) step=1; long ns=(n+step-1)/step; uint8_t* sc=(uint8_t*)malloc(ns); float* sv=(float*)malloc(sizeof(float)*ns);
  /* r9 speed (bit-identical): visit only the sampled indices i = j*step box row by box row, no n-byte box map */
  for(int k=0;k<nb;k++){ int* B=bx+9*k; uint8_t c=(uint8_t)clsof(B[6],B[7],L); float wk=w[k];
    for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ long r0=(long)z*sz+(long)y*nx+B[4], r1=r0+(B[5]-B[4]);
      for(long i=(r0+step-1)/step*step;i<r1;i+=step){ long j=i/step; sc[j]=c; sv[j]=fabsf(cn[i]*wk); } } }
  csort(sc,sv,ns);
  double lo=R*1e-9, hi=4*R; float off[16]; double tm;
  if(psnr>0) tm=R*R/pow(10.0,psnr/10.0)/1.1; else { double k_=kt>0? kt : auto_kt(sc,sv,ns,rnd,0,linf,lo,hi,0.65); tm=(linf/k_)*(linf/k_); }
  double d=rc_solve(sc,sv,ns,rnd,0,tm,lo,hi,off);
  free(sc); free(sv);
  /* full pass mirroring the coder's structure: 4x4x4 all-zero cell flags, 2x2x2 sub-cell flags inside non-empty cells,
     zeroth-order magnitude entropy of the voxels in non-empty sub-cells, 1 sign bit per nonzero; all per class */
  double hc[16][2], hs[16][2], hm[16][65], esc=0, nzc=0; memset(hc,0,sizeof(hc)); memset(hs,0,sizeof(hs)); memset(hm,0,sizeof(hm));
  int qb[64];
  for(int k=0;k<nb;k++){ int* B=bx+9*k; int c=clsof(B[6],B[7],L); float s_=(float)(w[k]/d);
    for(int z0=B[0];z0<B[1];z0+=4) for(int y0=B[2];y0<B[3];y0+=4) for(int x0=B[4];x0<B[5];x0+=4){ if(ZEHALF && (((z0-B[0])>>2)+((y0-B[2])>>2)+((x0-B[4])>>2))&1) continue; /* checkerboard half of the cells */ int any=0; memset(qb,0,sizeof(qb));
      if(x0+4<=B[5]){
#if WLC_AVX2
        const __m128 S_=_mm_set1_ps(s_), RN_=_mm_set1_ps((float)rnd), AB4=_mm_castsi128_ps(_mm_set1_epi32(0x7fffffff)); __m128i AN=_mm_setzero_si128(), ONE4=_mm_set1_epi32(1);
        for(int z=z0;z<z0+4 && z<B[1];z++) for(int y=y0;y<y0+4 && y<B[3];y++){ const float* r=cn+(long)z*sz+(long)y*nx+x0;
          __m128i q4=_mm_cvttps_epi32(_mm_floor_ps(_mm_add_ps(_mm_mul_ps(_mm_and_ps(_mm_loadu_ps(r),AB4),S_),RN_))); AN=_mm_or_si128(AN,q4);
          _mm_storeu_si128((__m128i*)(qb+((z-z0)*4+(y-y0))*4),_mm_add_epi32(q4,ONE4)); }
        AN=_mm_or_si128(AN,_mm_srli_si128(AN,8)); AN=_mm_or_si128(AN,_mm_srli_si128(AN,4)); any=_mm_cvtsi128_si32(AN);
#else
        for(int z=z0;z<z0+4 && z<B[1];z++) for(int y=y0;y<y0+4 && y<B[3];y++){ const float* r=cn+(long)z*sz+(long)y*nx; for(int x=x0;x<x0+4 && x<B[5];x++){
          int q=(int)floorf(fabsf(r[x])*s_+(float)rnd); qb[((z-z0)*4+(y-y0))*4+(x-x0)]= q+1; any|=q; } }
#endif
      }
      else for(int z=z0;z<z0+4 && z<B[1];z++) for(int y=y0;y<y0+4 && y<B[3];y++){ const float* r=cn+(long)z*sz+(long)y*nx; for(int x=x0;x<x0+4 && x<B[5];x++){
        int q=(int)floorf(fabsf(r[x])*s_+(float)rnd); qb[((z-z0)*4+(y-y0))*4+(x-x0)]= q+1; any|=q; } }
      hc[c][any!=0]++; if(!any) continue;
      for(int sz_=0;sz_<2;sz_++) for(int sy=0;sy<2;sy++) for(int sx=0;sx<2;sx++){ int sa=0, cntv=0, qs[8];
        for(int dz=0;dz<2;dz++) for(int dy=0;dy<2;dy++) for(int dx=0;dx<2;dx++){ int t=qb[((sz_*2+dz)*4+sy*2+dy)*4+sx*2+dx]; if(t){ qs[cntv++]=t-1; sa|=t-1; } }
        if(!cntv) continue; hs[c][sa!=0]++; if(!sa) continue;
        for(int u=0;u<cntv;u++){ int q=qs[u]; if(q>0) nzc++; if(q>=64 || q<0){ hm[c][64]++; esc+= q<0? 63.0 : 2*log2((double)q)+1; } else hm[c][q]++; } /* r9: q<0 = int overflow at sub-ulp tau (crashed r8pg) */ } } }
  double bits=esc+nzc; if(keep) *keep=cn; else if(!ZPOOL){ ZPOOL=cn; ZPOOLN=n; } else free(cn); if(getenv("WLC_ZTDBG")) fprintf(stderr,"ZE time dwt %.3f rest %.3f\n",t1_-t0_,cpun()-t1_);
  for(int c=0;c<16;c++){ double a=hc[c][0]+hc[c][1], b=hs[c][0]+hs[c][1], m=0; for(int t=0;t<65;t++) m+=hm[c][t];
    for(int t=0;t<2;t++){ if(hc[c][t]) bits-=hc[c][t]*log2(hc[c][t]/a); if(hs[c][t]) bits-=hs[c][t]*log2(hs[c][t]/b); }
    for(int t=0;t<65;t++) if(hm[c][t]) bits-=hm[c][t]*log2(hm[c][t]/m); }
  PKZ=pk0; return bits; }


WLC_API long wlc_encode(const uint16_t* v,int nz,int ny,int nx,double psnr,double linf,double kt,double rnd,int maxl,uint8_t* out,long cap,double* stats,double lam);
WLC_API long wlc_encode_f(const float* v,int nz,int ny,int nx,double psnr,double linf,double kt,double rnd,int maxl,uint8_t* out,long cap,double* stats,double lam);
static int ZTSKIP_=ZTSKIP;
#ifndef ZSHARE
#define ZSHARE 1
#endif
static int ZSHARE_=ZSHARE; void wlc_set_zshare(int s){ ZSHARE_=s; }
#ifndef ZDEC
#define ZDEC 1
#endif
#ifndef ZDM
#define ZDM 0.10
#endif
static double ZLASTMG=-1; static int ZLASTD=-1; double wlc_zlastmg(void){ return ZLASTMG; } int wlc_zlastd(void){ return ZLASTD; }
static int ZDECM_=1; void wlc_set_zdecm(int m){ ZDECM_=m; }
static int ZDEC_=ZDEC, ZINNER=0; static double ZMRG=1e9; void wlc_set_zdec(int s){ ZDEC_=s; }
static inline void zmrg(double e,double best){ double r=fabs(e/best-(1-ZTM)); if(r<ZMRG) ZMRG=r; }
/* forward dyadic DWT identical to dwt3_core(...,0) except that the level-0 x and y passes are assumed already done */
static void dwt_fwd_from0(float* a,int nz,int ny,int nx,int Lz,int Ly,int Lx){
  int L = Lz>Ly?Lz:Ly; if(Lx>L) L=Lx; int dz[32],dy[32],dx[32]; int cz=nz,cy=ny,cx=nx;
  for(int l=0;l<L;l++){ dz[l]=cz;dy[l]=cy;dx[l]=cx; if(l<Lz&&cz>1) cz=(cz+1)/2; if(l<Ly&&cy>1) cy=(cy+1)/2; if(l<Lx&&cx>1) cx=(cx+1)/2; }
  for(int l=0;l<L;l++){ if(l>0 && l<Lx) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],2,0); if(l>0 && l<Ly) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],1,0); if(l<Lz) axis_pass(a,ny,nx,dz[l],dy[l],dx[l],0,0); } }
static float* zcenter(const void* v,int isf,long n){ float* cn=(float*)balloc(sizeof(float)*n); float c0;
  if(isf){ const float* vf=(const float*)v; float a_,b_; fminmax(vf,n,&a_,&b_); c0=0.5f*(a_+b_); for(long i=0;i<n;i++) cn[i]=vf[i]-c0; }
  else { double mn=1e300,mx=-1e300; for(long i=0;i<n;i++){ double t=((const uint16_t*)v)[i]; if(t<mn)mn=t; if(t>mx)mx=t; }
    c0=(float)(0.5*(mn+mx)); for(long i=0;i<n;i++) cn[i]=(float)((const uint16_t*)v)[i]-c0; } return cn; }
static int ztrial(const void* v,int isf,int nz,int ny,int nx,double psnr,double linf,double kt,double rnd,int maxl,double lam,const int* Ls,double R){ int zlow=0;
  if(!ZTRIAL || !ZTON) return maxl; { const char* e=getenv("WLC_ZTRIAL"); if(e && !atoi(e)) return maxl; }
  if(!ZINNER){ if(nz<16 || Ls[0]<2 || ny<64 || nx<64 || (!isf && !ZTINT && !he_on())) return maxl; } /* integer (medical) volumes: estimate only in high-effort mode (speed) */
  if(!ZINNER){ /* cheap gate: only volumes whose z-roughness is at least the in-plane roughness (rms 1-step z diff vs 2-step x/y diff on a sparse grid) */
    long sz=(long)ny*nx; int zs=nz/32>0? nz/32:1, ys=ny/64>0? ny/64:1, xs=nx/64>0? nx/64:1; double dz=0,dy=0,dx=0;
    for(int z=0;z+1<nz;z+=zs) for(int y=0;y+2<ny;y+=ys) for(int x=0;x+2<nx;x+=xs){ long i=(long)z*sz+(long)y*nx+x;
      double c= isf? ((const float*)v)[i] : ((const uint16_t*)v)[i], a1= isf? ((const float*)v)[i+sz] : ((const uint16_t*)v)[i+sz],
             a2= isf? ((const float*)v)[i+2*nx] : ((const uint16_t*)v)[i+2*nx], a3= isf? ((const float*)v)[i+2] : ((const uint16_t*)v)[i+2];
      dz+=(a1-c)*(a1-c); dy+=(a2-c)*(a2-c); dx+=(a3-c)*(a3-c); }
    double r= (dx+dy>0)? sqrt(dz)/(0.5*(sqrt(dy)+sqrt(dx))) : 0; if(getenv("WLC_ZTDBG")) fprintf(stderr,"ZT ratio %.3f\n",r); if(r<ZTG && !(ZTLV && (Ls[0]<Ls[1] || Ls[0]<Ls[2]))) return maxl; zlow=(r<ZTG); } /* r8p: also run when z gets fewer dyadic levels than x/y (packet differs from dyadic) */
#if ZTE
  /* r8ze: whole-volume estimate instead of crop trials. For each candidate depth: forward DWT of the full volume, the same
     band weights, sampled |coef| per class, the same step-size solve as the encoder (PSNR target or auto-kt L-inf target),
     then the zeroth-order entropy of the quantised magnitudes per class (+1 sign bit per nonzero). Smallest estimate wins
     (> ZTM relative gain needed to leave the default). */
  { int xy= maxl<10? maxl : maxl%10; int cand[4]={maxl, 200+10*Ls[0]+xy, 100+xy, 210+xy}; double best=-1; int bm=maxl; long n_=(long)nz*ny*nx;
    if(ZDEC_ && ZSHARE_ && !ZINNER && ny>=128 && nx>=128){
      /* r9 speed (engineering, coarse-to-fine screening): run the same estimate on the 2x xy-decimated volume (every other
         row/column); if every comparison in that decision is clear of the 2% switching threshold by >= ZDM, take it (the main
         encode then computes the winner's transform itself); otherwise fall through to the full-volume estimate. */
      int nyd=ny/2, nxd=nx/2, it_=isf?4:2; long nd=(long)nz*nyd*nxd; uint8_t* vd=(uint8_t*)malloc(nd*it_); int crop= ZDECM_==2; int cy0=(ny-nyd)/2, cx0=(nx-nxd)/2;
      for(int z=0;z<nz;z++) for(int y=0;y<nyd;y++){ const uint8_t* src=(const uint8_t*)v+(((long)z*ny+(crop? cy0+y : 2*y))*nx)*it_; uint8_t* dst=vd+(((long)z*nyd+y)*nxd)*it_;
        if(crop){ memcpy(dst,src+(long)cx0*it_,(size_t)nxd*it_); continue; }
        if(isf){ const float* s4=(const float*)src; float* d4=(float*)dst; for(int x=0;x<nxd;x++) d4[x]=s4[2*x]; }
        else { const uint16_t* s2=(const uint16_t*)src; uint16_t* d2=(uint16_t*)dst; for(int x=0;x<nxd;x++) d2[x]=s2[2*x]; } }
      int shd[3]={nz,nyd,nxd}, Lsd[3]; levels_for(shd,Lsd,maxl); ZINNER=1; ZMRG=1e9; int z0n=ZTON; ZTON=1;
      int bmd=ztrial(vd,isf,nz,nyd,nxd,psnr,linf,kt,rnd,maxl,lam,Lsd,R); double mg=ZMRG; ZLASTMG=mg; ZLASTD=bmd; ZINNER=0; ZTON=z0n; free(vd); free(ZCACHE); ZCACHE=NULL;
      if(getenv("WLC_ZTDBG")) fprintf(stderr,"ZD decimated pick %d margin %.3f %s\n",bmd,mg,mg>=ZDM?"(taken)":"(full estimate)");
      if(mg>=ZDM && mg<1e8) return bmd; }
    int sh_[3]={nz,ny,nx}, Ld[3], Lp[3]; levels_for(sh_,Ld,cand[0]); levels_for(sh_,Lp,cand[2]);
    if(ZSHARE_ && Ld[1]==Lp[1] && Ld[2]==Lp[2] && Ld[1]>=1 && Ld[2]>=1){
      /* r9 speed (engineering): one shared xy transform. Level-0 x/y passes are common to the default and 'no z' transforms
         (bit-identical); the z packet is computed as z-DWT applied after the xy transform (separable ops commute; equal up to
         float rounding to the z-first order the decoder inverts). */
      double tq=cpun(); float* A=zcenter(v,isf,n_); axis_pass(A,ny,nx,nz,ny,nx,2,0); axis_pass(A,ny,nx,nz,ny,nx,1,0);
      float* B=(float*)balloc(sizeof(float)*n_); memcpy(B,A,sizeof(float)*n_); dwt_fwd_from0(B,nz,ny,nx,Ld[0],Ld[1],Ld[2]);
      dwt_fwd_from0(A,nz,ny,nx,0,Lp[1],Lp[2]);   /* A = xy-only transform ('no z' candidate) */
      if(getenv("WLC_ZTDBG")) fprintf(stderr,"ZS shared transforms %.3f\n",cpun()-tq);
      double e[4]; float* kb[4]={0,0,0,0};
      ZPRE=B; e[0]=zest_bits(v,isf,nz,ny,nx,cand[0],psnr,linf,kt,rnd,R,&kb[0]);
      { float* C=(float*)balloc(sizeof(float)*n_); memcpy(C,A,sizeof(float)*n_); int pz=pkz_of(nz,cand[1]); dwt3_core(C,nz,ny,nx,pz,0,0,0); pk_scale(C,nz,ny,nx,pz,0);
        ZPRE=C; e[1]=zest_bits(v,isf,nz,ny,nx,cand[1],psnr,linf,kt,rnd,R,&kb[1]); }
      best=e[0]; zmrg(e[1],best); if(e[1]<best*(1-ZTM)){ best=e[1]; bm=cand[1]; }
      if(bm==cand[0]){ free(kb[1]); kb[1]=NULL; } else { free(kb[0]); kb[0]=NULL; }
      ZPRE=A; e[2]=zest_bits(v,isf,nz,ny,nx,cand[2],psnr,linf,kt,rnd,R,&kb[2]);
      zmrg(e[2],best); if(e[2]<best*(1-ZTM)){ best=e[2]; bm=cand[2]; free(kb[0]); free(kb[1]); kb[0]=kb[1]=NULL; }
      if(bm==cand[2]){ float* C=(float*)balloc(sizeof(float)*n_); memcpy(C,kb[2],sizeof(float)*n_); int pz=pkz_of(nz,cand[3]); dwt3_core(C,nz,ny,nx,pz,0,0,0); pk_scale(C,nz,ny,nx,pz,0);
        ZPRE=C; e[3]=zest_bits(v,isf,nz,ny,nx,cand[3],psnr,linf,kt,rnd,R,&kb[3]);
        zmrg(e[3],best); if(e[3]<best*(1-ZTM)){ best=e[3]; bm=cand[3]; free(kb[2]); kb[2]=NULL; } else { free(kb[3]); kb[3]=NULL; } }
      else { free(kb[2]); kb[2]=NULL; }
      free(ZCACHE); ZCACHE=NULL; for(int k=0;k<4;k++) if(kb[k]){ ZCACHE=kb[k]; ZCM=cand[k]; }
      if(getenv("WLC_ZTDBG")){ fprintf(stderr,"ZE est %.0f %.0f %.0f bytes\n",e[0]/8,e[1]/8,e[2]/8); fprintf(stderr,"ZE %dx%dx%d p%.1f l%g -> maxl %d (shared)\n",nz,ny,nx,psnr,linf,bm); }
      return bm; }
    for(int k=0;k<4;k++){ if(k==3 && bm!=100+xy) break; if(k>=2 && zlow && ZTSKIP_) break; /* r9 speed: z smoother than x/y -> 'no z transform' is not tried */
      float* kb=NULL; double e=zest_bits(v,isf,nz,ny,nx,cand[k],psnr,linf,kt,rnd,R,&kb);
      if(k==0 || e < best*(1-ZTM)){ float* old=ZCACHE; ZCACHE=kb; ZCM=cand[k]; kb=old; } if(kb){ if(!ZPOOL){ ZPOOL=kb; ZPOOLN=n_; } else free(kb); }
      if(getenv("WLC_ZTDBG")) fprintf(stderr,"ZE maxl %d est %.0f bytes\n",cand[k],e/8);
      if(k==0) best=e; else if(e < best*(1-ZTM)){ best=e; bm=cand[k]; } }
    free(ZPOOL); ZPOOL=NULL; if(getenv("WLC_ZTDBG")) fprintf(stderr,"ZE %dx%dx%d p%.1f l%g -> maxl %d\n",nz,ny,nx,psnr,linf,bm); return bm; }
#endif
  int cy=ny/ZTC, cx=nx/ZTC; if(cy<64) cy=64; if(cx<64) cx=64; int y0=(ny-cy)/2, x0=(nx-cx)/2; long cn=(long)nz*cy*cx; int it=isf?4:2;
  uint8_t* c=(uint8_t*)malloc(cn*it); for(int z=0;z<nz;z++) for(int y=0;y<cy;y++) memcpy(c+((long)z*cy+y)*cx*it,(const uint8_t*)v+(((long)z*ny+y0+y)*nx+x0)*it,(size_t)cx*it);
  double mn=1e300,mx=-1e300; for(long i=0;i<cn;i++){ double t= isf? ((float*)c)[i] : ((uint16_t*)c)[i]; if(t<mn)mn=t; if(t>mx)mx=t; } double Rc=mx-mn; if(!(Rc>0)) { free(c); return maxl; }
  double pc= psnr>0? psnr+20*log10(Rc/R) : 0; long cap=cn*5+(1<<20); uint8_t* o=(uint8_t*)malloc(cap);
  int xy= maxl<10? maxl : maxl%10; long best=-1; int bm=maxl; ZTON=0; int fd=FDIR; if(isf) FDIR=0;
  /* default vs no z transform; only if "no z" wins by > ZTM also try one z level */
  int cand[3]={maxl, 100+xy, 10+xy};
  for(int k=0;k<3;k++){ if(k==2 && bm==maxl) break;
    long eab_=EABORT; EABORT=0; long L= isf? wlc_encode_f((const float*)c,nz,cy,cx,pc,linf,kt,rnd,cand[k],o,cap,NULL,lam) : wlc_encode((const uint16_t*)c,nz,cy,cx,pc,linf,kt,rnd,cand[k],o,cap,NULL,lam); EABORT=eab_;
    if(k==0) best=L; else if(L < best*(1-ZTM)){ best=L; bm=cand[k]; } }
  ZTON=1; FDIR=fd; free(o); free(c); if(getenv("WLC_ZTDBG")) fprintf(stderr,"ZT %dx%dx%d p%.1f l%g -> maxl %d\n",nz,ny,nx,psnr,linf,bm); return bm; }
WLC_API long wlc_encode(const uint16_t* v,int nz,int ny,int nx,double psnr,double linf,double kt,double rnd,int maxl,uint8_t* out,long cap,double* stats,double lam){
  double pt_=cpun(); long n=(long)nz*ny*nx, sz=(long)ny*nx; int sh[3]={nz,ny,nx}; int Ls[3]; levels_for(sh,Ls,maxl); int L=Lmax3(Ls[0],Ls[1],Ls[2]);
  int vmin=65535,vmax=0; { uint16_t a_=65535,b_=0; for(long i=0;i<n;i++){ uint16_t t_=v[i]; a_= t_<a_? t_:a_; b_= t_>b_? t_:b_; } vmin=a_; vmax=b_; }
  double R=vmax-vmin; if(R<1)R=1; int autok= (psnr<=0 && kt<=0);
  PKZ=0; { int m2=ztrial(v,0,nz,ny,nx,psnr,linf,kt,rnd,maxl,lam,Ls,R); if(m2!=maxl){ PKZ=pkz_of(nz,m2); maxl=m2; levels_for(sh,Ls,maxl); L=Lmax3(Ls[0],Ls[1],Ls[2]); } }
  float* cn; float c0=0.5f*(vmin+vmax); MK(8) double t0=now();
  if(ZCACHE && ZCM==maxl){ cn=ZCACHE; ZCACHE=NULL; } /* estimate used the same c0 = (vmin+vmax)/2 (exact) -> reuse its transform */
  else { free(ZCACHE); ZCACHE=NULL; cn=(float*)balloc(sizeof(float)*n); for(long i=0;i<n;i++) cn[i]=(float)v[i]-c0; dwt3(cn,nz,ny,nx,Ls[0],Ls[1],Ls[2],0); } TM[0]=now()-t0; t0=now();
  int bx[256*9]; int nb=boxes(nz,ny,nx,Ls[0],Ls[1],Ls[2],bx); float w[256]; box_w(bx,nb,Ls,w);
  for(int k=0;k<nb;k++){ int* B=bx+9*k; for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ float* p=cn+(long)z*sz+(long)y*nx; for(int x=B[4];x<B[5];x++) p[x]*=w[k]; } }
  MK(9) long step = n/600000; if(step<1) step=1; long ns=(n+step-1)/step;
  uint8_t* sc=(uint8_t*)malloc(ns); float* sv=(float*)malloc(sizeof(float)*ns); int32_t* q=(int32_t*)balloc(sizeof(int32_t)*n);
  { uint8_t* cm=(uint8_t*)q; for(int k=0;k<nb;k++){ int* B=bx+9*k; uint8_t c=(uint8_t)clsof(B[6],B[7],L);
      for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++) memset(cm+(long)z*sz+(long)y*nx+B[4],c,B[5]-B[4]); }
    for(long i=0,j=0;i<n;i+=step,j++){ sc[j]=cm[i]; sv[j]=fabsf(cn[i]); } csort(sc,sv,ns); }
  MK(10) /* predicted MSE with exact quantizer + per-class centroids (+1/12 integer rounding) */
  double lo=1e-3, hi=4*R+1; float off[16];
  MK(11) if(autok) kt=auto_kt(sc,sv,ns,rnd,1.0/12,linf,lo,hi,-0.15);
  double tm = psnr>0? R*R/pow(10.0,psnr/10.0)/1.1 : (linf/kt)*(linf/kt);
  MK(12) lo=rc_solve(sc,sv,ns,rnd,1.0/12,tm,lo,hi,off); MK(13)
  double rho_; { long c_=0; float id2=(float)(1.0/lo); for(long j=0;j<ns;j++) if(floorf(sv[j]*id2+(float)rnd)>0) c_++; rho_=(double)c_/ns; } RHO=rho_;
  free(sc); free(sv);
  TM[1]=now()-t0; t0=now(); int iters=1;
  double d=lo; int32_t* r=q; float* buf=cn; /* r aliases q: q is dead once recon has copied the filter output (RSCR=q) back into buf */ /* cn is dead after RDO/centroids */
  double mse=0; long L1=0; Hdr h; memset(&h,0,sizeof(h));
  { float id=(float)(1.0/d); double sm[16]={0}; long cnt[16]={0};
    for(int k=0;k<nb;k++){ int* B=bx+9*k; int c=clsof(B[6],B[7],L);
      for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ long i0=(long)z*sz+(long)y*nx+B[4]; qrow(cn+i0,q+i0,B[5]-B[4],id,(float)rnd,&sm[c],&cnt[c]); } }
    for(int c=0;c<16;c++){ double o= cnt[c]? sm[c]/cnt[c] : 0.5-rnd; if(o<-0.5)o=-0.5; if(o>0.5)o=0.5; off[c]=(float)o; }
    for(int it_=he_on();it_>0;it_--) he_step(cn,q,v,NULL,c0,nz,ny,nx,Ls,bx,nb,w,L,(float)d,off,(float)rnd,1.0f,99);
    GZW= getenv("WLC_ZW")? atoi(getenv("WLC_ZW")) : zw_choose(v,0,nz,ny,nx); MK(14) L1=enc_coefs(q,nz,ny,nx,Ls[0],Ls[1],Ls[2],out+sizeof(Hdr),cap-sizeof(Hdr), lam>0?cn:NULL, 1.0f, off, (float)lam);
    if(lam>0){ double sm2[16]={0}; long ct2[16]={0};
      for(int k=0;k<nb;k++){ int* B=bx+9*k; int c=clsof(B[6],B[7],L);
        for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ long i0=(long)z*sz+(long)y*nx+B[4]; orow(cn+i0,q+i0,B[5]-B[4],&sm2[c],&ct2[c]); } }
      for(int c=0;c<16;c++){ if(ct2[c]){ double o=sm2[c]/ct2[c]; if(o<-0.5)o=-0.5; if(o>0.5)o=0.5; off[c]=(float)o; } } }
    MK(15) RSCR=(float*)q; recon(q,buf,r,nz,ny,nx,Ls,bx,nb,w,(float)d,off,vmin,vmax,&h,v,FK); RSCR=NULL;
    /* fused error statistics: sum of squares (exact in int64 per row) and max |e| */
    long long s=0; int mx=0; for(long i0=0;i0<n;i0+=nx){ long long sr=0; int mr=0; for(long i=i0;i<i0+nx;i++){ int e=(int)v[i]-r[i]; sr+=(long long)e*e; int ae=e<0?-e:e; mr= ae>mr? ae:mr; } s+=sr; if(mr>mx)mx=mr; } mse=(double)s/n; MAXE=mx; }
  int tau; if(linf>0) tau=(int)floor(linf); else if(kt>0) tau=(int)ceil(kt*sqrt(mse>1e-9?mse:1e-9)); else tau=(int)MAXE; /* PSNR mode: guarantee = achieved max error, no correction layer */ if(tau<0)tau=0; int st=2*tau+1; long nout=0;
  if((int)MAXE>tau) for(long i=0;i<n;i++){ int e=(int)v[i]-r[i]; int kk=0; if(e>tau) kk=(e+tau)/st; else if(e<-tau) kk=-((-e+tau)/st); r[i]=kk; if(kk) nout++; }
  h.pad2=GZW; h.pad=PKZ; h.magic=0x434c5731; h.nz=nz;h.ny=ny;h.nx=nx;h.Lz=Ls[0];h.Ly=Ls[1];h.Lx=Ls[2];h.vmin=vmin;h.vmax=vmax;h.tau=tau;h.d=(float)d; memcpy(h.off,off,sizeof(off));
  TM[2]=now()-t0; TM[3]=iters; t0=now(); long p=sizeof(Hdr);
  h.l1=L1; p+=h.l1;
  if(nout){ uint8_t* ec=(uint8_t*)balloc(n); make_ec(buf,nz,ny,nx,tau,ec);
  Enc e; e_init(&e,out+p,cap-p); out_code(r,ec,nz,ny,nx,&e,NULL); free(ec); h.l2=e_finish(&e); p+=h.l2; } else h.l2=0;
  TM[4]=now()-t0; memcpy(out,&h,sizeof(h));
  if(stats){ stats[0]=mse; stats[1]=tau; stats[2]=nout; stats[3]=h.l1; stats[4]=h.l2; stats[5]=d; stats[6]=RHO; stats[7]=kt; }
  free(cn); free(q); return p;
}
WLC_API int wlc_shape(const uint8_t* in,int* sh){ Hdr h; memcpy(&h,in,sizeof(h)); sh[0]=h.nz; sh[1]=h.ny; sh[2]=h.nx; sh[3]=h.tau; return 0; }
WLC_API int wlc_decode(const uint8_t* in,long len,uint16_t* outv){
  Hdr h; memcpy(&h,in,sizeof(h)); int nz=h.nz,ny=h.ny,nx=h.nx; long n=(long)nz*ny*nx; int Ls[3]={h.Lz,h.Ly,h.Lx}; PKZ=h.pad;
  int32_t* q=(int32_t*)balloc(sizeof(int32_t)*n); float* buf=(float*)balloc(sizeof(float)*n);
  long p=sizeof(Hdr); GZW= h.pad2? h.pad2 : 2; dec_coefs(in+p,h.l1,nz,ny,nx,Ls[0],Ls[1],Ls[2],q); p+=h.l1;
  int bx[256*9]; int nb=boxes(nz,ny,nx,Ls[0],Ls[1],Ls[2],bx); float w[256]; box_w(bx,nb,Ls,w);
  RSCR=(float*)q; recon(q,buf,NULL,nz,ny,nx,Ls,bx,nb,w,h.d,h.off,h.vmin,h.vmax,&h,NULL,0); RSCR=NULL;
  memset(q,0,sizeof(int32_t)*n); if(h.l2>0){ Dec d; d_init(&d,in+p,h.l2); uint8_t* ec=(uint8_t*)balloc(n); make_ec(buf,nz,ny,nx,h.tau,ec); out_code(q,ec,nz,ny,nx,NULL,&d); free(ec); } int st=2*h.tau+1;
  /* fused: r = clamp(rint(buf+c)) (same rounding as recon's r path), then + q*(2tau+1), clamp, store uint16 - no n-int32 r array */
  float c=0.5f*(h.vmin+h.vmax); long i=0;
#if WLC_AVX2
  { __m256 C=_mm256_set1_ps(c); __m256i LO=_mm256_set1_epi32(h.vmin), HI=_mm256_set1_epi32(h.vmax), ST=_mm256_set1_epi32(st);
  for(;i+8<=n;i+=8){ __m256i t=_mm256_cvtps_epi32(_mm256_add_ps(_mm256_loadu_ps(buf+i),C)); t=_mm256_min_epi32(_mm256_max_epi32(t,LO),HI);
    t=_mm256_add_epi32(t,_mm256_mullo_epi32(_mm256_loadu_si256((const __m256i*)(q+i)),ST)); t=_mm256_min_epi32(_mm256_max_epi32(t,LO),HI);
    __m128i lo=_mm256_castsi256_si128(t), hi=_mm256_extracti128_si256(t,1); _mm_storeu_si128((__m128i*)(outv+i),_mm_packus_epi32(lo,hi)); } }
#endif
  for(;i<n;i++){ long r=lrintf(buf[i]+c); if(r<h.vmin)r=h.vmin; if(r>h.vmax)r=h.vmax; long t=r+(long)q[i]*st; if(t<h.vmin)t=h.vmin; if(t>h.vmax)t=h.vmax; outv[i]=(uint16_t)t; }
  free(q); free(buf); return 0;
}

/* ---------------- float32 path: same pipeline, quantized L-inf correction k*2tau ---------------- */
static inline float ffin(float xr,int k,float st){ return xr+(float)k*st; }
static float fstep(float tau,float mn,float mx){ float am=fabsf(mn)>fabsf(mx)?fabsf(mn):fabsf(mx); float u=am*(1.0f/8388608.0f); float st=2*tau*(1.0f-1.0f/1024)-6*u; if(st<tau*0.25f) st=tau*0.25f; return st; }
float* DBGO=NULL; WLC_API void wlc_dbg(float* p){ DBGO=p; }
#ifndef DMED
#define DMED 1
#endif
static int fcmp_(const void* a,const void* b){ float x=*(const float*)a, y=*(const float*)b; return (x>y)-(x<y); }
/* direct mode codes x-c0 with no transform: centre on the median, not the mid-range (one outlier would otherwise shift every index) */
static float med_sample(const float* v,long n){ long st=n/65536; if(st<1)st=1; long m=(n+st-1)/st; float* s=(float*)malloc(sizeof(float)*m); long j=0;
  for(long i=0;i<n&&j<m;i+=st) s[j++]=v[i]; qsort(s,j,sizeof(float),fcmp_); float r=s[j/2]; free(s); return r; }
static int FDIR=-1; static const uint8_t* MKM=NULL; static double MKSSE=0;
#ifndef TRIALC
#define TRIALC 0.8
#endif
static long wlc_encode_f_core(const float* v,int nz,int ny,int nx,double psnr,double linf,double kt,double rnd,int maxl,uint8_t* out,long cap,double* stats,double lam){
  double pt_=cpun();
  long n=(long)nz*ny*nx, sz=(long)ny*nx; int sh[3]={nz,ny,nx}; int Ls[3]; levels_for(sh,Ls,maxl); int direct=0; { const char* ev=getenv("WLC_DIRECT"); double rc_=lag1corr_f(v,nz,ny,nx); direct= FDIR>=0? FDIR : ev? atoi(ev) : (rc_<0.3);
    if(FDIR<0 && !ev && rc_>=0.3 && rc_<TRIALC){ /* ambiguous correlation: trial-code a leading slab both ways, keep the smaller */
      int tz=nz, ty=ny; if(nz>=32){ tz=nz/8; if(tz<16) tz=16; } else { ty=ny/8; if(ty<64) ty=ny<64?ny:64; }
      long tn=(long)tz*ty*nx; long tc=tn*5+(1<<20); uint8_t* tb=(uint8_t*)malloc(tc); long sz2[2];
      for(int dm=0;dm<2;dm++){ FDIR=dm; const float* src=v; float* tmp=NULL;
        if(ty<ny){ tmp=(float*)malloc(sizeof(float)*tn); for(int z=0;z<tz;z++) memcpy(tmp+(long)z*ty*nx,v+(long)z*ny*nx,sizeof(float)*ty*nx); src=tmp; }
        { long eab_=EABORT; EABORT=0; sz2[dm]=wlc_encode_f(src,tz,ty,nx,psnr,linf,kt,rnd,maxl,tb,tc,NULL,lam); EABORT=eab_; } free(tmp); }
      FDIR=-1; free(tb); direct= sz2[1]<sz2[0]; }
    if(direct){ Ls[0]=Ls[1]=Ls[2]=0; } } int L=Lmax3(Ls[0],Ls[1],Ls[2]); int dl=(direct && linf>0); if(dl){ rnd=0.5; lam=0; }
  float vmn,vmx; fminmax(v,n,&vmn,&vmx);
  double R=(double)vmx-vmn; if(R<=0)R=1e-30;
  PKZ=0; if(!direct){ int m2=ztrial(v,1,nz,ny,nx,psnr,linf,kt,rnd,maxl,lam,Ls,R); if(m2!=maxl){ maxl=m2; levels_for(sh,Ls,maxl); L=Lmax3(Ls[0],Ls[1],Ls[2]); PKZ=pkz_of(nz,maxl); } }
  int autok= (psnr<=0 && kt<=0); double tm=0;
  float c0=0.5f*(vmn+vmx); if(direct && DMED) c0=med_sample(v,n); float* cn;
  if(ZCACHE && ZCM==maxl && !direct){ cn=ZCACHE; ZCACHE=NULL; } /* transform already computed by the z-depth estimate */
  else { free(ZCACHE); ZCACHE=NULL; cn=(float*)balloc(sizeof(float)*n); for(long i=0;i<n;i++) cn[i]=v[i]-c0; dwt3(cn,nz,ny,nx,Ls[0],Ls[1],Ls[2],0); } MK(8) 
  int bx[256*9]; int nb=boxes(nz,ny,nx,Ls[0],Ls[1],Ls[2],bx); float w[256]; box_w(bx,nb,Ls,w);
  for(int k=0;k<nb;k++){ int* B=bx+9*k; for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ float* p=cn+(long)z*sz+(long)y*nx; for(int x=B[4];x<B[5];x++) p[x]*=w[k]; } }
  long step = n/600000; if(step<1) step=1; long ns=(n+step-1)/step;
  uint8_t* sc=(uint8_t*)malloc(ns); float* sv=(float*)malloc(sizeof(float)*ns); int32_t* q=(int32_t*)balloc(sizeof(int32_t)*n);
  { uint8_t* cm=(uint8_t*)q; for(int k=0;k<nb;k++){ int* B=bx+9*k; uint8_t c=(uint8_t)clsof(B[6],B[7],L);
      for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++) memset(cm+(long)z*sz+(long)y*nx+B[4],c,B[5]-B[4]); }
    for(long i=0,j=0;i<n;i+=step,j++){ sc[j]=cm[i]; sv[j]=fabsf(cn[i]); } csort(sc,sv,ns); }
  double lo=R*1e-9, hi=4*R; float off[16];
  if(dl){ lo=2.0*nextafterf((float)linf,0.0f)*(1.0-1.0/4096); kt=0; }
  else { if(autok) kt=auto_kt(sc,sv,ns,rnd,0,linf,lo,hi,0.65);
  tm = psnr>0? R*R/pow(10.0,psnr/10.0)/1.1 : (linf/kt)*(linf/kt);
  lo=rc_solve(sc,sv,ns,rnd,0,tm,lo,hi,off); }
  { float am=0; for(long i=0;i<n;i++){ float t=fabsf(cn[i]); am=t>am?t:am; } double fl=(double)am/1073741824.0; if(lo<fl) lo=fl; } /* r9 robustness: keep |q| < 2^30 (sub-ulp tau overflowed int32 in r8pg) */
  double rho_; { long c_=0; float id2=(float)(1.0/lo); for(long j=0;j<ns;j++) if(floorf(sv[j]*id2+(float)rnd)>0) c_++; rho_=(double)c_/ns; } RHO=rho_;
  free(sc); free(sv);
  MK(9) double d=lo; float* buf=cn; /* cn is dead after RDO/centroids */
  double mse=0; long L1=0; Hdr h; memset(&h,0,sizeof(h));
  { float id=(float)(1.0/d); double sm[16]={0}; long cnt[16]={0};
    for(int k=0;k<nb;k++){ int* B=bx+9*k; int c=clsof(B[6],B[7],L);
      for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ long i0=(long)z*sz+(long)y*nx+B[4]; qrow(cn+i0,q+i0,B[5]-B[4],id,(float)rnd,&sm[c],&cnt[c]); } }
    for(int c=0;c<16;c++){ double o= cnt[c]? sm[c]/cnt[c] : 0.5-rnd; if(o<-0.5)o=-0.5; if(o>0.5)o=0.5; off[c]= dl? 0.f : (float)o; }
    if(!direct && !dl) for(int it_=he_on();it_>0;it_--) he_step(cn,q,NULL,v,c0,nz,ny,nx,Ls,bx,nb,w,L,(float)d,off,(float)rnd,1.0f,99);
    GZW= getenv("WLC_ZW")? atoi(getenv("WLC_ZW")) : (direct? 2 : zw_choose(v,1,nz,ny,nx)); MK(10) L1=enc_coefs(q,nz,ny,nx,Ls[0],Ls[1],Ls[2],out+sizeof(Hdr),cap-sizeof(Hdr), lam>0?cn:NULL, 1.0f, off, (float)lam);
    if(EABORT && L1>EABORT){ free(cn); free(q); if(stats) memset(stats,0,8*sizeof(double)); return (long)sizeof(Hdr)+L1; }
    { double sm2[16]={0}; long ct2[16]={0};
      for(int k=0;k<nb;k++){ int* B=bx+9*k; int c=clsof(B[6],B[7],L);
        for(int z=B[0];z<B[1];z++) for(int y=B[2];y<B[3];y++){ long i0=(long)z*sz+(long)y*nx+B[4]; orow(cn+i0,q+i0,B[5]-B[4],&sm2[c],&ct2[c]); } }
      if(!dl) for(int c=0;c<16;c++){ if(ct2[c]){ double o=sm2[c]/ct2[c]; if(o<-0.5)o=-0.5; if(o>0.5)o=0.5; off[c]=(float)o; } } }
    MK(11) RSCR=(float*)q; recon2(q,buf,NULL,nz,ny,nx,Ls,bx,nb,w,(float)d,off,0,0,&h,NULL,v,dl?0:FK,c0); RSCR=NULL;
    double s=0; finrow(v,buf,n,c0,vmn,vmx,&s,&MAXE); mse=s/n; if(MKM){ double s2=0; for(long i=0;i<n;i++) if(!MKM[i]){ double e_=(double)v[i]-buf[i]; s2+=e_*e_; } MKSSE=s2; } }MK(12) 
  float tau; if(linf>0) tau=nextafterf((float)linf,0.0f); else if(kt>0) tau=(float)(kt*sqrt(mse>0?mse:1e-30)); else { tau=nextafterf(nextafterf(MAXE,INFINITY),INFINITY); /* float |err| max, 2 ulp safety */ } if(!(tau>0)) tau=(float)(R*1e-7); float st=fstep(tau,vmn,vmx); long nout=0;
  int32_t* kk=q; /* reuse */
  float tfast=tau*(1.0f-1.0f/65536);
#if WLC_AVX2
  const __m256 TF_=_mm256_set1_ps(tfast), AB_=_mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
#endif
  for(long i=0;i<n;i++){
#if WLC_AVX2
    if(!DBGO && (i&7)==0 && i+8<=n){ /* 8 voxels all within tfast -> k=0 (identical test, vectorised) */
      __m256 d_=_mm256_and_ps(_mm256_sub_ps(_mm256_loadu_ps(v+i),_mm256_loadu_ps(buf+i)),AB_);
      if(_mm256_movemask_ps(_mm256_cmp_ps(d_,TF_,_CMP_LE_OQ))==0xff){ _mm256_storeu_si256((__m256i*)(kk+i),_mm256_setzero_si256()); i+=7; continue; } }
#endif
    float xr=buf[i]; if(!DBGO && (fabsf(v[i]-xr)<=tfast || (MKM && MKM[i]))){ kk[i]=0; continue; } /* r9: floor voxels are restored exactly from the mask */ double e=(double)v[i]-(double)xr; int k=0;
    if(fabs(e)>tau){ k=(int)lrint(e/st); if(k==0) k=e>0?1:-1;
      for(int g=0;g<8 && fabs((double)v[i]-(double)ffin(xr,k,st))>tau;g++){ int k2= k+(((double)v[i]>(double)ffin(xr,k,st))?1:-1); if(k2==0) break; k=k2; }
      nout++; }
    kk[i]=k; if(DBGO) DBGO[i]= k? ffin(xr,k,st) : xr; }
  MK(13) h.pad2=GZW; h.pad=PKZ; h.magic=0x434c5731; h.isf=1; h.nz=nz;h.ny=ny;h.nx=nx;h.Lz=Ls[0];h.Ly=Ls[1];h.Lx=Ls[2]; h.fmn=vmn; h.fmx=vmx; h.ftau=tau; h.fc0=c0; h.d=(float)d; memcpy(h.off,off,sizeof(off));
  long p=sizeof(Hdr); h.l1=L1; p+=h.l1;
  if(nout){ uint8_t* ec=(uint8_t*)balloc(n); make_ec(buf,nz,ny,nx,tau,ec);
  Enc e; e_init(&e,out+p,cap-p); out_code(kk,ec,nz,ny,nx,&e,NULL); free(ec); h.l2=e_finish(&e); p+=h.l2; } else h.l2=0;
  MK(14) memcpy(out,&h,sizeof(h));
  if(stats){ stats[0]=mse; stats[1]=tau; stats[2]=nout; stats[3]=h.l1; stats[4]=h.l2; stats[5]=d; stats[6]=RHO; stats[7]=kt; }
  free(cn); free(q); return p;
}
WLC_API int wlc_info(const uint8_t* in,int* sh,float* tf){ uint32_t mg_; memcpy(&mg_,in,4); if(mg_==0x434c5733){ int32_t t[8]; memcpy(t,in,32); sh[0]=t[1]; sh[1]=t[2]; sh[2]=t[3]; sh[3]=0; sh[4]=1; memcpy(tf,in+24,4); return 0; }
  Hdr h; memcpy(&h,in,sizeof(h)); sh[0]=h.nz; sh[1]=h.ny; sh[2]=h.nx; sh[3]=h.tau; sh[4]=h.isf; tf[0]=h.ftau; return 0; }
static int wlc_decode_f_core(const uint8_t* in,long len,float* outv){
  Hdr h; memcpy(&h,in,sizeof(h)); int nz=h.nz,ny=h.ny,nx=h.nx; long n=(long)nz*ny*nx; int Ls[3]={h.Lz,h.Ly,h.Lx}; PKZ=h.pad;
  int32_t* q=(int32_t*)balloc(sizeof(int32_t)*n); float* buf=(float*)balloc(sizeof(float)*n);
  long p=sizeof(Hdr); GZW= h.pad2? h.pad2 : 2; dec_coefs(in+p,h.l1,nz,ny,nx,Ls[0],Ls[1],Ls[2],q); p+=h.l1;
  int bx[256*9]; int nb=boxes(nz,ny,nx,Ls[0],Ls[1],Ls[2],bx); float w[256]; box_w(bx,nb,Ls,w);
  RSCR=(float*)q; recon2(q,buf,NULL,nz,ny,nx,Ls,bx,nb,w,h.d,h.off,0,0,&h,NULL,NULL,0,h.fc0); RSCR=NULL;
  { /* r10: AVX clamp; without a correction layer it writes the output directly (identical values) */
    float* dst= h.l2>0? buf : outv; long i=0;
#if WLC_AVX2
    { const __m256 C8=_mm256_set1_ps(h.fc0), MN8=_mm256_set1_ps(h.fmn), MX8=_mm256_set1_ps(h.fmx);
    for(;i+8<=n;i+=8) _mm256_storeu_ps(dst+i,_mm256_min_ps(_mm256_max_ps(_mm256_add_ps(_mm256_loadu_ps(buf+i),C8),MN8),MX8)); }
#endif
    for(;i<n;i++){ float xr=buf[i]+h.fc0; if(xr<h.fmn)xr=h.fmn; if(xr>h.fmx)xr=h.fmx; dst[i]=xr; } }
  if(h.l2>0){ memset(q,0,sizeof(int32_t)*n); Dec d; d_init(&d,in+p,h.l2); uint8_t* ec=(uint8_t*)balloc(n); make_ec(buf,nz,ny,nx,h.ftau,ec); out_code(q,ec,nz,ny,nx,NULL,&d); free(ec);
    float st=fstep(h.ftau,h.fmn,h.fmx); long i=0;
#if WLC_AVX2
    { const __m256 ST8=_mm256_set1_ps(st); const __m256i Z8=_mm256_setzero_si256();
    for(;i+8<=n;i+=8){ __m256i k8=_mm256_loadu_si256((const __m256i*)(q+i)); __m256 xr=_mm256_loadu_ps(buf+i); __m256 y=_mm256_add_ps(xr,_mm256_mul_ps(_mm256_cvtepi32_ps(k8),ST8));
      _mm256_storeu_ps(outv+i,_mm256_blendv_ps(y,xr,_mm256_castsi256_ps(_mm256_cmpeq_epi32(k8,Z8)))); } }
#endif
    for(;i<n;i++) outv[i]= q[i]? ffin(buf[i],q[i],st) : buf[i]; }
  free(q); free(buf); return 0;
}

/* ---------------- r9: floor-mask mode for floor-dominated fields (engineering; prior art: MPEG-4 shape-adaptive padding,
   push-pull inpainting, JBIG-style context-coded binary masks). If >= MKF of the voxels equal the minimum ("floor", e.g. a
   clipped log mixing ratio), the floor mask is coded losslessly (3D causal-context binary coder), the floor region is filled
   smoothly from the non-floor data (push-pull pyramid) so the wavelet sees no cliff, the filled volume is coded as usual, and
   the decoder writes the floor value back exactly under the mask. L-inf: masked voxels exact, others within tau.
   Stream: inner WLC stream (magic CLW2) + mask stream + footer {int64 mask bytes, float floor}. ---------------- */
#ifndef MKBOTH
#define MKBOTH 1
#endif
#ifndef MKREF
#define MKREF 1
#endif
#ifndef MKG
#define MKG 0.58
#endif
#ifndef MKTOL
#define MKTOL 1.5
#endif
#ifndef MKF
#define MKF 0.5
#endif
#ifndef MKON
#define MKON 1
#endif
#define MAGIC2 0x434c5732
static int MKBUSY=0;
static void pushpull(float* a,const uint8_t* m,int nz,int ny,int nx){ /* m=1: unknown (floor) -> filled */
  int cz=(nz+1)/2, cy=(ny+1)/2, cx=(nx+1)/2; long sz=(long)ny*nx, csz=(long)cy*cx, cn=(long)cz*csz;
  float* s=(float*)calloc(cn,sizeof(float)); float* w=(float*)calloc(cn,sizeof(float));
  long known=0;
  for(int z=0;z<nz;z++) for(int y=0;y<ny;y++){ const float* r=a+(long)z*sz+(long)y*nx; const uint8_t* mr=m+(long)z*sz+(long)y*nx; long o=(long)(z>>1)*csz+(long)(y>>1)*cx;
    for(int x=0;x<nx;x++) if(!mr[x]){ s[o+(x>>1)]+=r[x]; w[o+(x>>1)]+=1; known++; } }
  if(!known){ free(s); free(w); for(long i=0;i<(long)nz*sz;i++) a[i]=0; return; }
  uint8_t* cm=(uint8_t*)malloc(cn); long unk=0; for(long i=0;i<cn;i++){ if(w[i]>0){ s[i]/=w[i]; cm[i]=0; } else { cm[i]=1; unk++; } }
  if(unk && (cz>1||cy>1||cx>1)) pushpull(s,cm,cz,cy,cx); else if(unk){ s[0]=0; }
  for(int z=0;z<nz;z++) for(int y=0;y<ny;y++){ float* r=a+(long)z*sz+(long)y*nx; const uint8_t* mr=m+(long)z*sz+(long)y*nx; const float* cr=s+(long)(z>>1)*csz+(long)(y>>1)*cx;
    for(int x=0;x<nx;x++) if(mr[x]) r[x]=cr[x>>1]; }
  free(s); free(w); free(cm); }
static inline int mctx(const uint8_t* m,int z,int y,int x,int ny,int nx){ long sz=(long)ny*nx; long i=(long)z*sz+(long)y*nx+x; int c=0;
  if(x>0) c|=m[i-1]; if(y>0) c|=m[i-nx]<<1; if(z>0) c|=m[i-sz]<<2; if(x>0&&y>0) c|=m[i-nx-1]<<3; if(z>0&&y>0) c|=m[i-sz-nx]<<4; if(z>0&&x>0) c|=m[i-sz-1]<<5;
  if(x>1) c|=m[i-2]<<6; if(y>1) c|=m[i-2*nx]<<7; return c; }
static long mask_code(uint8_t* m,int nz,int ny,int nx,Enc* e,Dec* d){ Prob P[256]; pinit(P,256);
  for(int z=0;z<nz;z++) for(int y=0;y<ny;y++) for(int x=0;x<nx;x++){ long i=((long)z*ny+y)*nx+x; int c=mctx(m,z,y,x,ny,nx);
    if(e) e_bit(e,&P[c],m[i]); else m[i]=(uint8_t)d_bit(d,&P[c]); }
  return 0; }
static long wlc_encode_f_mk(const float* v,int nz,int ny,int nx,double psnr,double linf,double kt,double rnd,int maxl,uint8_t* out,long cap,double* stats,double lam){
  long n=(long)nz*ny*nx; const char* ev=getenv("WLC_MASK"); int on= ev? atoi(ev) : MKON;
  if(MKBUSY || !on || n<4096) { int b=MKBUSY; MKBUSY=1; long r=wlc_encode_f_core(v,nz,ny,nx,psnr,linf,kt,rnd,maxl,out,cap,stats,lam); MKBUSY=b; return r; }
  float fl,fx; fminmax(v,n,&fl,&fx); long cnt=0; for(long i=0;i<n;i++) cnt+= (v[i]==fl);
  if(!(cnt>=MKF*n) || cnt==n){ MKBUSY=1; long r=wlc_encode_f_core(v,nz,ny,nx,psnr,linf,kt,rnd,maxl,out,cap,stats,lam); MKBUSY=0; return r; }
  uint8_t* m=(uint8_t*)malloc(n); float* f=(float*)balloc(sizeof(float)*n);
  for(long i=0;i<n;i++){ m[i]=(v[i]==fl); f[i]=v[i]; }
  pushpull(f,m,nz,ny,nx); float bl,bx_; fminmax(f,n,&bl,&bx_);
  double R=(double)fx-fl, Rb=(double)bx_-bl; if(!(Rb>0)) Rb=R;
  double pi=0; if(psnr>0){ double T=R*R/pow(10.0,psnr/10.0); pi=10*log10(Rb*Rb/T) - MKG*10*log10((double)n/(n-cnt)); } /* only the n-cnt unmasked voxels count; error is denser there than in the smooth fill (MKG<1, set on dev data) */
  MKBUSY=1; MKM=m; long L=wlc_encode_f_core(f,nz,ny,nx,pi,linf,kt,rnd,maxl,out,cap,stats,lam);
  if(psnr>0 && MKREF){ /* one rate-control refinement on the true (unmasked) error, measured inside the encoder */
    double ach= MKSSE>0? 10*log10(R*R/(MKSSE/n)) : 999; if(getenv("WLC_MKDBG")) fprintf(stderr,"MK floor %.4f frac %.3f psnr target %.2f inner %.2f -> achieved %.2f\n",fl,(double)cnt/n,psnr,pi,ach);
    if(ach<psnr-0.5 || ach>psnr+MKTOL){ /* undershoot is fixed; overshoot tolerated up to MKTOL dB (CPU) */ pi+= psnr-ach; L=wlc_encode_f_core(f,nz,ny,nx,pi,linf,kt,rnd,maxl,out,cap,stats,lam); if(getenv("WLC_MKDBG")) fprintf(stderr,"MK refined -> %.2f\n",10*log10(R*R/(MKSSE/n))); } }
  MKBUSY=0; MKM=NULL; Hdr h; memcpy(&h,out,sizeof(h)); h.magic=MAGIC2; memcpy(out,&h,sizeof(h));
  Enc e; e_init(&e,out+L,cap-L-12); mask_code(m,nz,ny,nx,&e,NULL); int64_t ml=e_finish(&e); L+=ml; memcpy(out+L,&ml,8); memcpy(out+L+8,&fl,4); L+=12;
  if(getenv("WLC_MKDBG")) fprintf(stderr,"MK mask %ld bytes (%.4f bpv)\n",(long)ml,ml*8.0/n);
  free(m); free(f);
  if(MKBOTH){ /* r9b: the smooth fill only pays at coarse/mid rates (at tight tau the non-constant fill itself costs bits):
       also code the volume the plain r8pg way and keep the better stream. Linf: fewer bytes (both meet tau).
       PSNR: of the streams reaching target-0.5 dB the smaller; if neither, the higher PSNR. */
    double mse_m= MKSSE/n; double st2[8]; long cp=cap; uint8_t* o2=(uint8_t*)malloc(cp);
    MKBUSY=1; long Lp=wlc_encode_f_core(v,nz,ny,nx,psnr,linf,kt,rnd,maxl,o2,cp,st2,lam); MKBUSY=0;
    int plain;
    if(psnr>0){ double pm=10*log10(R*R/(mse_m>0?mse_m:1e-300)), pp=10*log10(R*R/(st2[0]>0?st2[0]:1e-300)); int okm=pm>=psnr-0.5, okp=pp>=psnr-0.5;
      plain= (okm&&okp)? (Lp<L) : (okm||okp)? okp : (pp>pm); if(getenv("WLC_MKDBG")) fprintf(stderr,"MK both: mask %ld@%.2f plain %ld@%.2f -> %s\n",L,pm,Lp,pp,plain?"plain":"mask"); }
    else { plain= Lp<L; if(getenv("WLC_MKDBG")) fprintf(stderr,"MK both: mask %ld plain %ld -> %s\n",L,Lp,plain?"plain":"mask"); }
    if(plain){ memcpy(out,o2,Lp); if(stats) memcpy(stats,st2,sizeof(st2)); L=Lp; }
    free(o2); }
  return L; }
static int wlc_decode_f_mk(const uint8_t* in,long len,float* outv){ Hdr h; memcpy(&h,in,sizeof(h)); if(h.magic!=MAGIC2) return wlc_decode_f_core(in,len,outv);
  int64_t ml; float fl; memcpy(&ml,in+len-12,8); memcpy(&fl,in+len-4,4); long Li=len-12-ml; wlc_decode_f_core(in,Li,outv);
  long n=(long)h.nz*h.ny*h.nx; uint8_t* m=(uint8_t*)calloc(n,1); Dec d; d_init(&d,in+Li,ml); mask_code(m,h.nz,h.ny,h.nx,NULL,&d);
  for(long i=0;i<n;i++) if(m[i]) outv[i]=fl; free(m); return 0; }

WLC_API int wlc_zpick(const void* v,int isf,int nz,int ny,int nx,double psnr,double linf,int maxl,int skip){ int sh[3]={nz,ny,nx}, Ls[3]; levels_for(sh,Ls,maxl);
  double R; if(isf){ float a,b; fminmax((const float*)v,(long)nz*ny*nx,&a,&b); R=(double)b-a; if(R<=0)R=1e-30; } else { int a=65535,b=0; long n=(long)nz*ny*nx; for(long i=0;i<n;i++){ int t=((const uint16_t*)v)[i]; a=t<a?t:a; b=t>b?t:b; } R=b-a; if(R<1)R=1; }
  int s0=ZTSKIP_; ZTSKIP_=skip; int m=ztrial(v,isf,nz,ny,nx,psnr,linf,(psnr>0?0.0:-1.0),0.38,maxl,0.08,Ls,R); ZTSKIP_=s0; free(ZCACHE); ZCACHE=NULL; return m; }

WLC_API int wlc_zpick2(const void* v,int isf,int nz,int ny,int nx,double psnr,double linf,int maxl,int share){ int s0=ZSHARE_; ZSHARE_=share; int m=wlc_zpick(v,isf,nz,ny,nx,psnr,linf,maxl,0); ZSHARE_=s0; return m; }

WLC_API int wlc_zpick3(const void* v,int isf,int nz,int ny,int nx,double psnr,double linf,int maxl,int dec){ ZLASTMG=-1; ZLASTD=-1; int s0=ZDEC_; ZDEC_=dec; int m=wlc_zpick(v,isf,nz,ny,nx,psnr,linf,maxl,0); ZDEC_=s0; return m; }

/* ---------------- r10: PALETTE / PREDICTIVE mode for low-cardinality volumes (engineering; prior art: palette coding
   (PNG, HEVC-SCC palette mode), near-lossless predictive coding with context-modelled residuals (JPEG-LS/LOCO-I, CALIC,
   SZ Lorenzo predictor), value clustering for an L-inf bound). Exact distinct values -> sorted palette -> uint16 index volume.
   L-inf mode: palette entries are merged into clusters of width <= 2*tau (each value within tau of its cluster
   representative, verified per entry); PSNR mode: lossless. Indices coded in raster order: flag "equals left neighbour"
   (16 equality-pattern contexts), else the residual to a per-volume predictor (3D Lorenzo / left / 2D plane / z-linear,
   chosen on a sample) with activity contexts: zero flag, sign, Elias-gamma magnitude. Used only if a cheap sample test
   passes (<= PLK distinct values in a 1/64 sample and estimated size < PLR x the wavelet stream); the encoder keeps
   whichever stream is actually smaller. Stream: PHdr (magic CLW3) + palette floats + range-coded indices. ---------------- */
#define MAGIC3 0x434c5733
#ifndef PLON
#define PLON 1
#endif
#ifndef PLK
#define PLK 8192
#endif
#ifndef PLR
#define PLR 1.15
#endif
#ifndef PLPS
#define PLPS 62
#endif
#ifndef PLPQ
#define PLPQ 1
#endif
#ifndef PLSK
#define PLSK 1.08
#endif
#ifndef PLEQT
#define PLEQT 0.5
#endif
typedef struct { uint32_t magic; int32_t nz,ny,nx,K,pm; float tau; int32_t pad; } PHdr;
static int u32cmp(const void* a,const void* b){ float x=*(const float*)a, y=*(const float*)b; return x<y?-1:(x>y?1:0); }
/* one pass: distinct float bit patterns (<= maxk, else -1) and a raw id per voxel (if I != NULL). Open-addressing hash
   sized from hs_log2; returns -2 if the table would pass half full (caller retries bigger). Runs of equal voxels skip the hash. */
static int pal_scan(const float* v,long n,long step,uint32_t* pal,int maxk,uint16_t* I,int hlog){
  int HS=1<<hlog; uint32_t* hk=(uint32_t*)malloc(sizeof(uint32_t)*HS); int32_t* hv=(int32_t*)malloc(sizeof(int32_t)*HS); for(int t=0;t<HS;t++) hv[t]=-1;
  int K=0; uint32_t lastk=0; int lastid=-1; const uint32_t* u=(const uint32_t*)v;
  for(long i=0;i<n;i+=step){ uint32_t k=u[i]; if(k==lastk && lastid>=0){ if(I) I[i]=(uint16_t)lastid; continue; }
    if((k&0x7f800000u)==0x7f800000u && (k&0x007fffffu)){ K=-1; break; } /* NaN */
    uint32_t h=(k*2654435761u)>>(32-hlog); while(hv[h]>=0 && hk[h]!=k) h=(h+1)&(HS-1);
    if(hv[h]<0){ if(K>=maxk){ K=-1; break; } if(2*(K+1)>HS){ K=-2; break; } hv[h]=K; hk[h]=k; pal[K++]=k; }
    lastk=k; lastid=hv[h]; if(I) I[i]=(uint16_t)lastid; }
  free(hk); free(hv); return K; }
static const uint32_t* PSRT; static int idcmp(const void* a,const void* b){ float x,y; memcpy(&x,PSRT+*(const int*)a,4); memcpy(&y,PSRT+*(const int*)b,4); return x<y?-1:(x>y?1:(*(const int*)a-*(const int*)b)); }
static inline int pl_pred(const uint16_t* I,long i,int x,int y,int z,int nx,long sz,int pm,int K,int* act,int* ea){
  int a,b,c,d,e,f,g;
  if(x&&y&&z){ a=I[i-1]; b=I[i-nx]; c=I[i-sz]; d=I[i-nx-1]; e=I[i-sz-1]; f=I[i-sz-nx]; g=I[i-sz-nx-1]; }
  else { a= x? I[i-1] : (y? I[i-nx] : (z? I[i-sz] : 0)); b= y? I[i-nx] : a; c= z? I[i-sz] : a;
    d= (x&&y)? I[i-nx-1] : b; e= (x&&z)? I[i-sz-1] : c; f= (y&&z)? I[i-sz-nx] : c; g=(x&&y&&z)? I[i-sz-nx-1] : d; }
  int p; if(pm==0) p=a+b+c-d-e-f+g; else if(pm==1) p=a; else if(pm==2) p=a+b-d; else p= z>1? 2*c-I[i-2*sz] : c;
  if(p<0)p=0; if(p>K-1)p=K-1; unsigned t=(unsigned)(abs(a-d)+abs(b-d)+abs(c-g)+abs(a-e)); int q= t? 32-__builtin_clz(t) : 0; if(q>11) q=11; *act=q;
  *ea=(a==b)+2*(a==c)+4*(a==d); return a | (p<<16); }
typedef struct { Prob eq[64][2]; Prob eb[16][2]; Prob ec[16][2]; Prob zf[12][2]; Prob sg[12]; Prob eg[12][24]; Prob mb[24][2]; } PModel;
/* PV bit0: richer equal-to-left context (adds z-edge and y-edge equalities); bit1: equal-to-above / equal-to-previous-slice
   flags before the residual. Chosen per volume (blocky -> 3, natural 8-bit seismic -> 0), stored in PHdr.pad. */
static void pal_code(uint16_t* I,int nz,int ny,int nx,int K,int pm,int PV,Enc* e,Dec* d){
  PModel* M=(PModel*)malloc(sizeof(PModel)); pinit((Prob*)M,sizeof(PModel)/sizeof(Prob)); long sz=(long)ny*nx; int peq=0;
  for(int z=0;z<nz;z++) for(int y=0;y<ny;y++){ if(ABT) break; long i=(long)z*sz+(long)y*nx; for(int x=0;x<nx;x++,i++){
    int act,ea; int r_=pl_pred(I,i,x,y,z,nx,sz,pm,K,&act,&ea); int a=r_&0xffff, p=r_>>16;
    int bb= y? I[i-nx] : a, cc= z? I[i-sz] : a; int ee=(x&&z)? I[i-sz-1] : cc, dd=(x&&y)? I[i-nx-1] : bb;
    int ctx= (PV&1)? (ea|((cc==ee)<<3)|((bb==dd)<<4)|(peq<<5)) : (ea|(peq<<3));
    int isq; if(e){ isq= I[i]==a; e_bit(e,&M->eq[ctx][0],isq); } else isq=d_bit(d,&M->eq[ctx][0]);
    peq=isq; if(isq){ if(d) I[i]=(uint16_t)a; continue; }
    if(PV&2){ int c2=ea|((cc==ee)<<3);
      if(bb!=a){ int f; if(e){ f=I[i]==bb; e_bit(e,&M->eb[c2][0],f); } else f=d_bit(d,&M->eb[c2][0]); if(f){ if(d) I[i]=(uint16_t)bb; continue; } }
      if(cc!=a && cc!=bb){ int f; if(e){ f=I[i]==cc; e_bit(e,&M->ec[c2][0],f); } else f=d_bit(d,&M->ec[c2][0]); if(f){ if(d) I[i]=(uint16_t)cc; continue; } } }
    int r= e? (int)I[i]-p : 0; int nzr; if(e){ nzr=r!=0; e_bit(e,&M->zf[act][0],nzr); } else nzr=d_bit(d,&M->zf[act][0]);
    if(!nzr){ if(d) I[i]=(uint16_t)p; continue; }
    int neg; if(e){ neg=r<0; e_bit(e,&M->sg[act],neg); } else neg=d_bit(d,&M->sg[act]);
    unsigned u= e? (unsigned)(neg?-r:r) : 0; int nb=0;
    if(e){ while((u>>nb)>1) nb++; for(int t=0;t<nb;t++) e_bit(e,&M->eg[act][t],1); e_bit(e,&M->eg[act][nb],0);
      for(int t=nb-1;t>=0;t--){ int bit=(u>>t)&1; if(t==nb-1) e_bit(e,&M->mb[nb][0],bit); else if(t==nb-2) e_bit(e,&M->mb[nb][1],bit); else e_bitp(e,bit,32768); } }
    else { while(nb<23 && d_bit(d,&M->eg[act][nb])) nb++; u=1; for(int t=nb-1;t>=0;t--){ int bit; if(t==nb-1) bit=d_bit(d,&M->mb[nb][0]); else if(t==nb-2) bit=d_bit(d,&M->mb[nb][1]); else bit=d_bitp(d,32768); u=(u<<1)|bit; }
      int rr= neg? -(int)u : (int)u; I[i]=(uint16_t)(p+rr); } } }
  free(M); }
/* cheap estimate of the coded size: zeroth-order cost of the same symbols on every 8th row (bits) */
static double PLEQF=0; /* fraction of voxels equal to the left neighbour (last pal_est) */
static double pal_est(const uint16_t* I,int nz,int ny,int nx,int K,int pm){
  long sz=(long)ny*nx; double c_eq[16][2]={{0}}, c_z[12][2]={{0}}, mag=0, rows=0;
  for(int z=0;z<nz;z++) for(int y=(z*5)&15;y<ny;y+=16){ rows++; long i=(long)z*sz+(long)y*nx; int peq=0; for(int x=0;x<nx;x++,i++){
    int act,ea; int r_=pl_pred(I,i,x,y,z,nx,sz,pm,K,&act,&ea); int a=r_&0xffff,p=r_>>16; int isq=I[i]==a; c_eq[ea|(peq<<3)][isq]++; peq=isq; if(isq) continue;
    int r=(int)I[i]-p; c_z[act][r!=0]++; if(r){ unsigned u=abs(r); int nb=0; while((u>>nb)>1) nb++; mag+= 1+ (nb? 2*nb : 1) + 0.0; } } }
  double b=mag; for(int k=0;k<16;k++){ double t=c_eq[k][0]+c_eq[k][1]; for(int j=0;j<2;j++) if(c_eq[k][j]) b-=c_eq[k][j]*log2(c_eq[k][j]/t); }
  for(int k=0;k<12;k++){ double t=c_z[k][0]+c_z[k][1]; for(int j=0;j<2;j++) if(c_z[k][j]) b-=c_z[k][j]*log2(c_z[k][j]/t); }
  { double q1=0,tt=0; for(int k=0;k<16;k++){ q1+=c_eq[k][1]; tt+=c_eq[k][0]+c_eq[k][1]; } PLEQF= tt>0? q1/tt : 0; }
  return rows>0? b*((double)nz*ny/rows) : 1e30; }
static int PLDBG=-1;
/* returns stream length or -1 (not applicable / not worth it); lw = size of the competing wavelet stream */
/* lim>0: stop coding once the index stream passes lim bytes (the wavelet stream is already smaller); result then > lim */
static long pal_encode(const float* v,int nz,int ny,int nx,double psnr,double linf,uint8_t* out,long cap,long lw,long lim){
  long n=(long)nz*ny*nx; if(n<4096 || nz<2 || ny<2 || nx<2) return -1; if(PLDBG<0) PLDBG= getenv("WLC_PLDBG")!=NULL;
  double t0_=(double)clock()/CLOCKS_PER_SEC;
  uint32_t* pb=(uint32_t*)malloc(sizeof(uint32_t)*65536); long step=n/(1L<<20); if(step<1) step=1;
  int Ks=pal_scan(v,n,step,pb,PLK,NULL,15); if(PLDBG) fprintf(stderr,"PL sample distinct %d\n",Ks); if(Ks<0){ free(pb); return -1; }
  uint16_t* I=(uint16_t*)balloc(sizeof(uint16_t)*n); int hl=12; while((1<<hl)<8*Ks && hl<17) hl++;
  int K; for(;;){ K=pal_scan(v,n,1,pb,65535,I,hl); if(K!=-2) break; hl++; } if(K<1){ free(pb); free(I); return -1; }
  int* ord=(int*)malloc(sizeof(int)*K); for(int k=0;k<K;k++) ord[k]=k; PSRT=pb; qsort(ord,K,sizeof(int),idcmp);
  float* pal=(float*)malloc(sizeof(float)*K); for(int k=0;k<K;k++) memcpy(&pal[k],&pb[ord[k]],4);
  /* clusters over the sorted values: L-inf merges entries within 2*tau (checked: every member within tau of the representative); PSNR: identity */
  uint16_t* cid=(uint16_t*)malloc(sizeof(uint16_t)*K); float* rep=(float*)malloc(sizeof(float)*K); int Kc=0; double tau= linf>0? linf : 0;
  double* hc=NULL, tmse=0;
  if(psnr>0 && linf<=0 && K>1 && PLPQ){ /* PSNR mode: widest cluster tau whose exact MSE (from the value histogram, midpoint reps) meets the target */
    long* cnt=(long*)calloc(K,sizeof(long)); for(long i=0;i<n;i++) cnt[I[i]]++; hc=(double*)malloc(sizeof(double)*K); for(int k=0;k<K;k++) hc[k]=(double)cnt[ord[k]]; free(cnt);
    double R=(double)pal[K-1]-(double)pal[0], tm=R*R/pow(10.0,psnr/10.0)*n, lo=0, hi=R;
    for(int it=0;it<40;it++){ double t=0.5*(lo+hi); double se=0;
      for(int k=0;k<K;){ int j=k; while(j+1<K && (double)pal[j+1]-(double)pal[k] <= 2*t*(1-1e-6)) j++; double r=(float)(0.5*((double)pal[k]+(double)pal[j]));
        for(int u=k;u<=j;u++){ double e=(double)pal[u]-r; se+=hc[u]*e*e; } k=j+1; if(se>tm) break; }
      if(se<=tm) lo=t; else hi=t; }
    tmse=tm; tau=lo; if(PLDBG) fprintf(stderr,"PL psnr-mode cluster tau %.4g\n",tau); }
  for(;;){ Kc=0;
  for(int k=0;k<K;){ int j=k; if(tau>0){ while(j+1<K && (double)pal[j+1]-(double)pal[k] <= 2*tau*(1-1e-6)) j++; }
    for(;;){ float r=(float)(0.5*((double)pal[k]+(double)pal[j])); int ok=1; for(int t=k;t<=j;t++) if(fabs((double)pal[t]-(double)r)>tau){ ok=0; break; }
      if(ok || j==k){ rep[Kc]= (j==k)? pal[k] : r; break; } j--; }
    for(int t=k;t<=j;t++) cid[t]=(uint16_t)Kc; Kc++; k=j+1; }
  if(!hc || tau<=0) break; { double se=0; for(int k=0;k<K;k++){ double e=(double)pal[k]-(double)rep[cid[k]]; se+=hc[k]*e*e; } if(se<=tmse) break; tau*=0.97; } }
  free(hc);
  uint16_t* map=(uint16_t*)malloc(sizeof(uint16_t)*K); for(int k=0;k<K;k++) map[ord[k]]=cid[k];
  for(long i=0;i<n;i++) I[i]=map[I[i]]; free(map); free(ord); free(cid); free(pal); free(pb);
  if(PLDBG) fprintf(stderr,"PL t_index %.2f\n",(double)clock()/CLOCKS_PER_SEC-t0_);
  int pm=0; double be=1e300; for(int m=0;m<4;m++){ double e_=pal_est(I,nz,ny,nx,Kc,m); if(PLDBG) fprintf(stderr,"PL K %d clusters %d pred %d est %.0f bytes (wavelet %ld)\n",K,Kc,m,e_/8,lw); if(e_<be){ be=e_; pm=m; } }
  if(PLDBG) fprintf(stderr,"PL t_est %.2f\n",(double)clock()/CLOCKS_PER_SEC-t0_);
  long L=-1;
  if(lw<=0 || be/8 + 4.0*Kc + sizeof(PHdr) < PLR*lw){
    int pv= PLEQF>PLEQT? 3 : 0; if(getenv("WLC_PLPV")) pv=atoi(getenv("WLC_PLPV")); if(PLDBG) fprintf(stderr,"PL eqfrac %.3f pv %d\n",PLEQF,pv);
    if(lim>0 && nz>=64){ /* cheap test (wavelet already finished): real coder on 4-slice slabs every 32 slices; skip if clearly larger */
      int ns=0; for(int z=0;z+4<=nz;z+=32) ns+=4; long ssz=(long)ny*nx; uint16_t* J=(uint16_t*)malloc(sizeof(uint16_t)*ns*ssz); int o_=0;
      for(int z=0;z+4<=nz;z+=32){ memcpy(J+(long)o_*ssz,I+(long)z*ssz,sizeof(uint16_t)*4*ssz); o_+=4; }
      long cs=(long)ns*ssz*2+(1<<16); uint8_t* sb=(uint8_t*)malloc(cs); Enc se; e_init(&se,sb,cs); long eab=EABORT; EABORT=0;
      /* slabs are coded as one stack: the 4-slice blocks keep most of the z context */
      pal_code(J,ns,ny,nx,Kc,pm,pv,&se,NULL); EABORT=eab; double est=(double)e_finish(&se)*nz/ns + sizeof(PHdr)+4.0*Kc; free(sb); free(J);
      if(PLDBG) fprintf(stderr,"PL subest %.0f vs wavelet %ld\n",est,lim);
      if(est > PLSK*lim){ free(I); free(rep); return -1; } }
    PHdr h={MAGIC3,nz,ny,nx,Kc,pm,(float)tau,pv}; long p=sizeof(PHdr)+4L*Kc;
    if(p<cap){ memcpy(out,&h,sizeof(h)); memcpy(out+sizeof(h),rep,4L*Kc); Enc e; e_init(&e,out+p,cap-p); long eab=EABORT; EABORT= lim>0? (lim-p>0? lim-p : 1) : 0; pal_code(I,nz,ny,nx,Kc,pm,pv,&e,NULL); EABORT=eab; long l2=e_finish(&e); L=p+l2; if(L>cap) L=cap+1; } }
  if(PLDBG) fprintf(stderr,"PL t_code %.2f\n",(double)clock()/CLOCKS_PER_SEC-t0_);
  free(I); free(rep); return L; }
static int pal_decode(const uint8_t* in,long len,float* outv){ PHdr h; memcpy(&h,in,sizeof(h)); long n=(long)h.nz*h.ny*h.nx;
  const float* rep=(const float*)(in+sizeof(h)); long p=sizeof(h)+4L*h.K; uint16_t* I=(uint16_t*)balloc(sizeof(uint16_t)*n);
  Dec d; d_init(&d,in+p,len-p); pal_code(I,h.nz,h.ny,h.nx,h.K,h.pm,h.pad,NULL,&d); for(long i=0;i<n;i++) outv[i]=rep[I[i]]; free(I); return 0; }
/* r10 wrapper: palette candidate first (only if the 1M-voxel sample has <= PLK distinct values), then the wavelet/mask encoder
   with an exact early abort once its coefficient stream alone exceeds the palette stream. Keeps the smaller stream. */
WLC_API long wlc_encode_f(const float* v,int nz,int ny,int nx,double psnr,double linf,double kt,double rnd,int maxl,uint8_t* out,long cap,double* stats,double lam){
  const char* ev=getenv("WLC_PAL"); if(MKBUSY || !(ev? atoi(ev) : PLON)) return wlc_encode_f_mk(v,nz,ny,nx,psnr,linf,kt,rnd,maxl,out,cap,stats,lam);
  long n=(long)nz*ny*nx; long c2= n*2+65536*4+(1<<16); uint8_t* o2=(uint8_t*)malloc(c2); long L2=-1, L, eab=EABORT;
  /* order only affects speed (each side aborts exactly once it is larger than the finished other side): low PSNR targets,
     where a (lossless) palette rarely wins, run the wavelet first */
  int wfirst= (psnr>0 && psnr<PLPS); if(getenv("WLC_PLORD")) wfirst=atoi(getenv("WLC_PLORD"));
  if(wfirst){ L=wlc_encode_f_mk(v,nz,ny,nx,psnr,linf,kt,rnd,maxl,out,cap,stats,lam); if(L>0) L2=pal_encode(v,nz,ny,nx,psnr,linf,o2,c2,0,L); }
  else { L2=pal_encode(v,nz,ny,nx,psnr,linf,o2,c2,0,0); if(L2>0) EABORT= L2-(long)sizeof(Hdr)>0 ? L2-(long)sizeof(Hdr) : 1;
    L=wlc_encode_f_mk(v,nz,ny,nx,psnr,linf,kt,rnd,maxl,out,cap,stats,lam); EABORT=eab; }
  if(PLDBG>0) fprintf(stderr,"PL wavelet %ld palette %ld order %s -> %s\n",L,L2,wfirst?"W":"P",(L2>0&&L2<L)?"palette":"wavelet");
  if(L2>0 && L2<L && L2<=cap){ memcpy(out,o2,L2); if(stats) memset(stats,0,8*sizeof(double)); L=L2; } free(o2); return L; }
WLC_API int wlc_decode_f(const uint8_t* in,long len,float* outv){ uint32_t m; memcpy(&m,in,4); if(m==MAGIC3) return pal_decode(in,len,outv); return wlc_decode_f_mk(in,len,outv); }
