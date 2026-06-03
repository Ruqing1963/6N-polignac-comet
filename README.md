# 6N Polignac Comet under omega-Stratification (Part XIV)

Every prime gap's comet collapses to **one** constant after dividing by the
Hardy-Littlewood singular series. The difference-side dual of the Goldbach
collapse (Part XIII).

**The Polignac comet.** For an even gap d, pi_d(X) = #{n<=X : n, n+d both prime}.
By Hardy-Littlewood,

```
    pi_d(X) ~ S(d) * X/ln^2(X),
    S(d) = 2*Pi_2 * prod_{p|d, p>2} (p-1)/(p-2),
    Pi_2 = prod_{p>2}(1 - 1/(p-1)^2)   (twin constant, 0.66016182).
```

For d=2 (twins) S = 2*Pi_2; for d=6=2·3 the odd factor 3 gives (3-1)/(3-2)=2, so
pi_6 ≈ 2·pi_2 — the familiar excess of gap-6 pairs. Plotting pi_d(X) against d
gives the "Polignac comet": horizontal bands indexed by omega_odd(d), each band a
value of S(d). This is the exact dual of the Goldbach comet (sum 2N, bands by
omega(N)): here the variable is the gap d, bands by omega_odd(d).

**Universal collapse (X = 10⁸, all 500 even d ≤ 1000).** With base(X)=X/ln²X,

```
    C(d) = pi_d(X) / (base(X) * S(d))
```

is a single constant: **mean 1.1313, CV 0.068%**, with no trend in d (quartile
means agree to 0.015%) and within-omega_odd(d) band CV ≤ 0.11%. Dividing by S(d)
collapses the comet to a flat line. The collapse is stronger than the sum side in
one respect: **the constant is shared by every gap.** Gaps with wildly different
absolute counts — pi_2=440,312, pi_6=879,908, pi_210=1,409,148 at X=10⁸ — all
reduce to C≈1.131 once S(d) is removed.

| omega_odd(d) | # gaps | mean C | band CV |   | d quartile | mean C |
|-------------:|-------:|-------:|--------:|---|-----------:|-------:|
| 0 |   9 | 1.13147 | 0.114% |   | [2,252)    | 1.13142 |
| 1 | 247 | 1.13141 | 0.076% |   | [252,502)  | 1.13130 |
| 2 | 221 | 1.13124 | 0.056% |   | [502,752)  | 1.13136 |
| 3 |  23 | 1.13150 | 0.033% |   | [752,1002) | 1.13127 |
| all | 500 | 1.13134 | 0.068% | |            |         |

**Sum–difference symmetry.** Parts I–XII resolved the prime difference at fixed gap
(twins d=2; cousins/sexy d=4,6). Part XIII collapsed the prime sum (Goldbach). This
completes the pair: across all even gaps the difference-side comet collapses onto
S(d) with a gap-independent constant. The sum-side enrichment prod(q-1)/(q-2) in N
and the difference-side prod(q-1)/(q-2) in d are the same singular-series mechanism
on the two linear constraints p1+p2=2N and p1−p2=d.

> **Attribution.** S(d) is the **classical Hardy-Littlewood prime-pair singular
> series**, NOT a new formula. The contribution is the omega-stratified reading of
> the comet's bands as S(d), the pointwise collapse to CV 0.068% over 500 gaps at
> X=10⁸, and the empirical **cross-gap constancy** of C(d) — a falsifiable
> statement (the quartile test) the data confirm.
>
> **Scope.** Only the conditional count's shape and constant. NO statement about
> Polignac's conjecture (that each even gap recurs infinitely often).
>
> **Log correction.** C carries the logarithmic correction of base=X/ln²X (it
> shifts with X: 1.147 at X=2e7, 1.131 at X=1e8), but — unlike the Goldbach decile
> drift — this does NOT break the cross-gap constancy, since X is fixed while d
> varies.

Part I: doi:10.5281/zenodo.20470367 · IX: doi:10.5281/zenodo.20520492 ·
XIII: doi:10.5281/zenodo.20530812

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   ├── polignac_data.csv       omega_odd, n_gaps, meanC, CV_pct  (X=1e8)
│   ├── polignac_quartile.csv   d_lo, d_hi, meanC, n  (cross-d constancy)
│   └── polignac_sample.csv     d, omega_odd, count, S/2Pi2, C  (sample gaps)
├── code/
│   ├── polignac.py         sieves, counts pi_d(X) for all even d by shifted
│   │                       intersection, evaluates S(d), reports collapse CV,
│   │                       omega bands, cross-d quartile trend; emits CSVs
│   └── make_polignac_fig.py self-contained 3-panel figure (rebuilds a small-X
│                            comet scatter; reads ../data for the quartile means)
├── figures/                fig_paper14_polignac.{pdf,png}
└── paper/                  Chen_6N_Paper14.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `sympy`, `matplotlib`.

```bash
pip install numpy sympy matplotlib

# 1. Collapse test. Default X=1e8, DMAX=1000 (~30 s, ~100-200 MB).
#    Emits polignac_data_X{X}.csv, polignac_quartile_X{X}.csv, polignac_sample_X{X}.csv.
python code/polignac.py
X=20000000 DMAX=300 python code/polignac.py     # quick

# 2. Figure (rebuilds a small-X scatter standalone; reads ../data quartiles).
cd code && python make_polignac_fig.py
```

### Conventions

- pi_d(X) = #{n<=X : n, n+d both prime}, counted by intersecting the prime
  indicator with its shift by d.
- S(d) = 2*Pi_2*prod_{p|d,p>2}(p-1)/(p-2); Pi_2 = prod_{p>2}(1-1/(p-1)^2).
- omega_odd(d) = number of distinct odd prime factors of d.
- base(X) = X/ln²X; C(d) = pi_d(X)/(base*S(d)).

## License

MIT — see `LICENSE`.
