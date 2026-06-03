#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 3-panel Polignac figure. Panels 1-2 rebuild a small-X comet scatter
standalone (raw pi_d/base coloured by omega_odd(d), and collapsed C(d)); panel 3
reads the cross-d quartile constancy from ../data/polignac_quartile.csv.
Self-contained: scatter rebuilt at X_FIG; statistics in the paper are at X=1e8.
"""
import numpy as np, math, csv
from sympy import factorint
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
X_FIG=20_000_000; DMAX=1000
def sieve_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return s
sieve=sieve_upto(X_FIG+DMAX+2)
Pi2=1.0
for p in np.nonzero(sieve[:10**6+1])[0]:
    if p>2: Pi2*=(1-1/(p-1)**2)
def Sfac(d):
    s=2*Pi2
    for p in factorint(int(d)):
        if p>2: s*=(p-1)/(p-2)
    return s
pm=sieve[:X_FIG+1]; base=X_FIG/math.log(X_FIG)**2
ds=[];raw=[];coll=[];om=[]
for d in range(2,DMAX+1,2):
    nh=X_FIG-d; cnt=int(np.count_nonzero(pm[2:nh+1]&sieve[2+d:nh+1+d]))
    S=Sfac(d); ds.append(d); raw.append(cnt/base); coll.append(cnt/(base*S))
    om.append(sum(1 for p in factorint(int(d)) if p>2))
ds=np.array(ds);raw=np.array(raw);coll=np.array(coll);om=np.array(om)
fig,axes=plt.subplots(1,3,figsize=(17,4.8))
sc=axes[0].scatter(ds,raw,c=om,cmap='viridis',s=10,alpha=.7)
axes[0].set_xlabel('gap $d$',fontsize=11); axes[0].set_ylabel(r'$\pi_d(X)/\mathrm{base}$',fontsize=11)
axes[0].set_title('Polignac comet: prime-pair count / base\n(bands = $\\omega_{\\mathrm{odd}}(d)$, i.e. $\\mathfrak{S}_d$)',fontsize=11)
plt.colorbar(sc,ax=axes[0],label=r'$\omega_{\mathrm{odd}}(d)$'); axes[0].grid(alpha=.2)
axes[1].scatter(ds,coll,c=om,cmap='viridis',s=10,alpha=.7)
axes[1].axhline(np.median(coll),color='gray',ls=':',lw=1.4)
axes[1].set_xlabel('gap $d$',fontsize=11); axes[1].set_ylabel(r'$C(d)=\pi_d/(\mathrm{base}\cdot\mathfrak{S}_d)$',fontsize=11)
axes[1].set_title('Divided by $\\mathfrak{S}_d$: all gaps collapse to one $C$\n(CV $0.068\\%$ at $X{=}10^8$, $d$ up to 1000)',fontsize=11)
axes[1].set_ylim(np.median(coll)*0.96,np.median(coll)*1.04); axes[1].grid(alpha=.2)
qd=list(csv.DictReader(open('../data/polignac_quartile.csv')))
mids=[(int(r['d_lo'])+int(r['d_hi']))/2 for r in qd]; cq=[float(r['meanC']) for r in qd]
axes[2].plot(mids,cq,'o-',color='#c0392b',lw=2,ms=10)
axes[2].axhline(np.mean(cq),color='gray',ls=':',lw=1.2,label=f'mean {np.mean(cq):.5f}')
axes[2].set_ylim(1.1305,1.1320)
axes[2].set_xlabel('gap $d$ (quartile midpoint)',fontsize=11); axes[2].set_ylabel(r'mean $C(d)$',fontsize=11)
axes[2].set_title('Cross-$d$ constancy: no trend\n($X{=}10^8$, quartiles differ $<0.015\\%$)',fontsize=11)
axes[2].legend(fontsize=9); axes[2].grid(alpha=.25)
plt.suptitle('Polignac comet vs the Hardy--Littlewood singular series $\\mathfrak{S}_d=2\\Pi_2\\prod_{p\\mid d,p>2}(p-1)/(p-2)$: every gap collapses to one constant',fontsize=12,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper14_polignac.pdf',bbox_inches='tight')
plt.savefig('fig_paper14_polignac.png',dpi=160,bbox_inches='tight')
print("figure saved")
