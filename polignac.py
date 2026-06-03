#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polignac comet collapse at large X. For each even gap d, pi_d(X) = #{n<=X :
n, n+d both prime}. Hardy-Littlewood predicts
    pi_d(X) ~ S(d) * X/ln^2(X),  S(d) = 2*Pi_2 * prod_{p|d, p>2} (p-1)/(p-2),
with Pi_2 = prod_{p>2}(1-1/(p-1)^2) the twin constant. We test whether
    C(d) = pi_d(X) / (base(X) * S(d))
collapses to a single constant across ALL even d (the difference-side unification,
dual to the Goldbach sum-side collapse of Part XIII), report the pointwise CV, the
within-omega_odd(d) band CVs, and whether C(d) has any residual trend in d.

ATTRIBUTION: S(d) is the classical Hardy-Littlewood singular series, not a new
formula. We verify the collapse and the cross-d constancy; we make no statement
about Polignac's conjecture (existence of infinitely many prime pairs at each
gap). Default X=1e8, DMAX=1000. Requires: numpy, sympy.
"""
import numpy as np, math, os, time
from sympy import factorint
def sieve_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return s
X=int(os.environ.get("X",100_000_000))
DMAX=int(os.environ.get("DMAX",1000))
t0=time.time()
sieve=sieve_upto(X+DMAX+2)
print(f"sieved to {X+DMAX:,} ({time.time()-t0:.0f}s)")
# Pi_2
Pi2=1.0
for p in np.nonzero(sieve[:10**6+1])[0]:
    if p>2: Pi2*=(1-1/(p-1)**2)
print(f"Pi_2 = {Pi2:.8f}  (literature 0.66016182)")
def Sfac(d):
    s=2*Pi2
    for p in factorint(int(d)):
        if p>2: s*=(p-1)/(p-2)
    return s
def omega_odd(d):
    return sum(1 for p in factorint(int(d)) if p>2)
primemask=sieve[:X+1]
base=X/math.log(X)**2
rows=[]
for d in range(2,DMAX+1,2):
    n_hi=X-d
    cnt=int(np.count_nonzero(primemask[2:n_hi+1] & sieve[2+d:n_hi+1+d]))
    S=Sfac(d); rows.append((d,omega_odd(d),cnt,S/(2*Pi2),cnt/(base*S)))
print(f"counted {len(rows)} gaps ({time.time()-t0:.0f}s)")
Cs=np.array([r[4] for r in rows])
print(f"\n(a) Collapse C=pi_d/(base*S_d), d=2..{DMAX}: mean={Cs.mean():.5f}, CV={100*Cs.std()/Cs.mean():.3f}%")
from collections import defaultdict
byom=defaultdict(list)
for d,om,cnt,sf,C in rows: byom[om].append(C)
print(f"\n(b) by omega_odd(d):")
print(f"{'omega_odd(d)':>13}{'#gaps':>7}{'meanC':>10}{'CV%':>7}")
for om in sorted(byom):
    v=np.array(byom[om]); print(f"{om:>13}{len(v):>7}{v.mean():>10.5f}{100*v.std()/v.mean():>7.3f}")
print(f"\n(c) C(d) trend in d (quartiles of gap size; flat => d-independent):")
ds=np.array([r[0] for r in rows])
for qi in range(4):
    lo=2+qi*DMAX//4; hi=2+(qi+1)*DMAX//4
    sel=(ds>=lo)&(ds<hi)
    if sel.any(): print(f"   d in [{lo},{hi}): meanC={Cs[sel].mean():.5f}  (n={sel.sum()})")
print(f"\nsample gaps:")
print(f"{'d':>5}{'om':>4}{'count':>12}{'S/2Pi2':>9}{'C':>10}")
for d,om,cnt,sf,C in rows:
    if d in [2,4,6,8,12,30,42,210,420,630,DMAX] or d<=12:
        print(f"{d:>5}{om:>4}{cnt:>12}{sf:>9.3f}{C:>10.5f}")
print(f"\ntotal {time.time()-t0:.0f}s")

# ---- emit CSV (polignac_data_X{X}.csv by-omega, polignac_quartile_X{X}.csv, polignac_sample_X{X}.csv) ----
import csv as _csv
with open(f'polignac_data_X{X}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['omega_odd','n_gaps','meanC','CV_pct'])
    for om in sorted(byom):
        v=np.array(byom[om]); _w.writerow([om,len(v),f'{v.mean():.5f}',f'{100*v.std()/v.mean():.3f}'])
with open(f'polignac_quartile_X{X}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['d_lo','d_hi','meanC','n'])
    for qi in range(4):
        lo=2+qi*DMAX//4; hi=2+(qi+1)*DMAX//4; sel=(ds>=lo)&(ds<hi)
        if sel.any(): _w.writerow([lo,hi,f'{Cs[sel].mean():.5f}',int(sel.sum())])
with open(f'polignac_sample_X{X}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['d','omega_odd','count','S_over_2Pi2','C'])
    for d,om,cnt,sf,C in rows: _w.writerow([d,om,cnt,f'{sf:.3f}',f'{C:.5f}'])
print(f"\n[ok] wrote polignac_data_X{X}.csv, polignac_quartile_X{X}.csv, polignac_sample_X{X}.csv")
