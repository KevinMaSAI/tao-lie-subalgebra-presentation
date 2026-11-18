# Code for "The relations of a Lie subalgebra considered by Tao: a complete presentation"

Author: Qiyuan Ma (Shanghai Jiao Tong University), version 1.1, November 16, 2026.
Code: MIT License (see LICENSE).

These scripts are an independent computational check of the proof in the paper (Section 12).
The proof itself does not depend on them.

Requirements: Python 3 (tested with 3.13.5); `sympy` (tested with 1.14.0) for
`pbw_count.py`, `egf_check.py`, `lformula.py` and `prior_work_check.py`. Exact computations use `fractions`.
Modular computations use p = 2^31 - 1.

Notation: n is the size of the index set of a multilinear component (content [n]),
and k is the E-degree (the number of letters G_T in the presentation of the ideal k).

| Script | Command | What it does | Output file | Runtime (recorded) |
|---|---|---|---|---|
| `kmodel_exact.py` | `python3 kmodel_exact.py 8` | dims of k, l, h and of ker(theta) in every component (n,k), exact over Q | `out_exact_7.txt`, `out_exact_8.txt` | about 110 s for n = 8 |
| `kmodel.py` | `python3 kmodel.py 8` | same, mod p | `out_kmodel_7.txt`, `out_kmodel_8.txt` | about 74 s for n = 8 |
| `gmodel.py` | `python3 gmodel.py 6` | independent model of g in the letters D_i, E_i (mod p): dim g and dim h | `out_gmodel_6.txt` | about 3 s |
| `ktop.py` | `python3 ktop.py 8` | dims of k^top (mod p) | `out_ktop_8.txt` | about 35 s |
| `pbw_check.py` | `python3 pbw_check.py` | degree-3 Groebner (PBW) check for U(g) and its quadratic dual, all 8 natural letter orders, exact | `out_pbw_check.txt` | < 1 s |
| `sanity_random.py` | `python3 sanity_random.py` | the degree-3 test of `pbw_check.py` applied to random quadratic relations (it must fail; shows the test is not vacuous) | `out_sanity_random.txt` | < 1 s |
| `ambiguities.py` | `python3 ambiguities.py` | explicit resolution of the 4 overlap ambiguities (Proposition 7.1), exact | `out_ambiguities.txt` | < 1 s |
| `pbw_count.py` | `python3 pbw_count.py 10` | counts multilinear normal words of U(g), takes log, compares with x(y+z) - log(1-Y) | `out_pbw_count_10.txt` | about 3 s |
| `egf_check.py` | `python3 egf_check.py` | identity (5) of the paper and Eulerian numbers, n <= 10 | `out_egf_check.txt` | about 3 s |
| `lformula.py` | `python3 lformula.py 10` | dims of k and l from the proven formulas (Cor. 6.2, Prop. 9.1) | `out_lformula_10.txt` | about 1 s |
| `kproj.py` | `python3 kproj.py 9` | direct injectivity check in degree 9 (mod p) via projection to big letters | `out_kproj_9.txt` | 554 s, peak RSS 403 MB |
| `prior_work_check.py` | `python3 prior_work_check.py` | dictionary of Section 11: the relators r_{i,j,S} equal the relators of tHyperCom [DST, Prop. 2.4.2] (n <= 6); dims of U(k) in content [n] by D-degree, exact over Q (n <= 6), compared with Losev-Manin's Poincare polynomials [LM, Thm. 2.3]; sympy check that EGF(U(g)) e^{-xy} at z=1 equals (y-1)/(y-e^{(y-1)x}) | `out_prior_work_check.txt` | about 41 s |
| `character.py`, `decompose.py` | `python3 character.py N` | exploratory S_n-character computations; not used in the paper | (none recorded) | |

Notes.
* `out_kmodel_9.txt` records that a run of `kmodel.py 9` was stopped at about 1.67 GB RSS to protect
  the machine's memory; `kproj.py` replaces it.
* In `out_exact_*.txt` and `out_kmodel_*.txt` the last column ("kernel") is dim ker(theta); it is 0
  in every component.
* The output files were produced on 2026-11-15 (Asia/Shanghai); `out_prior_work_check.txt` on 2026-11-16.
  On 2026-11-16 the scripts were re-run (kmodel_exact 7 and 8, ktop 8, gmodel 6, pbw_check, sanity_random,
  ambiguities, pbw_count 10, egf_check, lformula 10, kproj 9) and reproduced the recorded outputs.
* Version 1.1 changes: added `prior_work_check.py` and its output; corrected the docstring of
  `sanity_random.py` (for random relations dim I_3 is *larger* than the number of non-normal words).
