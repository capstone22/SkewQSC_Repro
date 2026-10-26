# SkewQSC_Repro

Reproducibility package for the manuscript:

**Synchronization Recovery via Ore Polynomial Orbits: A Skew-Cyclic Algebraic
Framework for Quantum Codes**

Authors:

- Gwang-Sik Kim
- Woosang Im
- Young-Sik Kim

Affiliation:

Department of Electrical Engineering and Computer Science, DGIST (Daegu
Gyeongbuk Institute of Science and Technology), Daegu, South Korea.

Corresponding author:

Young-Sik Kim, `ysk@dgist.ac.kr`

Funding:

This work was partly supported by Institute of Information & Communications
Technology Planning & Evaluation (IITP) grant funded by the Korea government
(MSIT) (RS-2026-25530181, Development of Inter-Domain Quantum Security System
Combining PQC and QKD, 50%) and by the InnoCORE
program of the Ministry of Science and ICT (26-InnoCORE-01, 50%).

## Review access

Editors and reviewers need access to these files during peer review. Cite a
fixed commit or release of this repository and confirm that reviewers can open
it, or submit the matching `skewqsc_repro_package.zip` as supplementary
material. A link to a private repository grants no access by itself. The
manuscript and response letter should name the route actually provided.

The code uses only the Python standard library.

## Reproduce the reported computations

Run these commands from this directory. They cover the numerical results
requested in referee comment 6.

| Manuscript result | Command | Where to check |
|-------------------|---------|----------------|
| Table 2, especially the last two rows | `python3 reproduce_search.py` and `python3 verify.py` | Fixed-parameter maxima in the printed output; code dimensions, distances, and decoders in the 46 checks. |
| Table 3, all 14 rows | `python3 reproduce_search.py` | Printed assertions and `landscape_data.csv`. |
| Figure 4, all 24 plotted points | `python3 reproduce_search.py` | `figure4_points.csv`, derived from the 80 degree-level rows in `landscape_data.csv`. |
| Positive-rate search behind Table 2 | `python3 search2.py` | 4,382 hits in `positive_hits.csv` and `hits.pkl`. |

The exact `(n,k1,k2)` values for the four Table 2 rows appear below.

## Files

| File | Role |
|------|------|
| `ffield.py` | Finite fields, row reduction, dual spaces, and brute-force minimum distance. |
| `ore.py` | Ore-polynomial arithmetic, right division, skew index, ordinary order, and skew-cyclic generator matrices. |
| `mindist.py` | Faster minimum-distance computation from parity-check column dependencies. |
| `decoder.py` | Bounded-distance q-ary syndrome decoder. |
| `performance.py` | Monte Carlo bit/phase decoding performance for the flagship positive-rate codes. |
| `search2.py` | Exhaustive search for positive-rate strict-gain skew CSS families; regenerates both hit files. |
| `verify.py` | Focused construction, decoder, and selected landscape checks. |
| `reproduce_search.py` | Exhaustive fixed-parameter checks for Table 2 and degree-level checks for Table 3 and Figure 4. |
| `landscape_data.csv` | The 80 degree-level search results from which the Figure 4 series are obtained. |
| `figure4_points.csv` | The 24 length-level maxima plotted in Figure 4. |
| `hits.pkl` | Pickled search hits from the positive-rate strict-gain search. |
| `positive_hits.csv` | Human-readable version of the positive-rate search hits, including distances and generator polynomials. |

## Quick start

Run the construction, decoder, and selected landscape checks:

```bash
python3 verify.py
```

The expected final line is:

```text
RESULT:  46/46 checks passed   ALL CHECKS PASSED
```

Run the performance simulation:

```bash
python3 performance.py
```

The script reports bit-decoder, phase-decoder, and combined success estimates
for four positive-rate flagship CSS codes at several independent symbol-error
probabilities.

Recompute all four matched maxima in Table 2, all 14 rows of Table 3, and all
24 plotted points of Figure 4 (about 20 seconds on the test machine):

```bash
python3 reproduce_search.py
```

The script asserts the printed maxima against the manuscript values and
regenerates `landscape_data.csv` and `figure4_points.csv`. The last two skew
code witnesses in Table 2, their dimensions, dual containment, distances, and
weight-one decoders are checked separately by `verify.py`.

Regenerate the positive-rate search hits (about one minute on the test machine):

```bash
python3 search2.py
```

The expected total is 4382 strict-gain hits in the stated `F4`, `F8`, and
`F9` positive-rate search ranges. This regenerates both `hits.pkl` and the
readable `positive_hits.csv`; running either script changes generated files in
the working directory. Do not unpickle data from untrusted sources.

`search2.py` tests lengths divisible by the Frobenius order up to `n=16`
over `F4`, `n=12` over `F8`, and `n=10` over `F9`. It enumerates all monic
right-divisor generators `g1` with `k1>n/2`, keeps the Euclidean
dual-containing inner codes, and enumerates all nested `g2` and their
quotients. A hit has a positive logical dimension and a skew index strictly
larger than the cyclic optimum at the same `(n,k1,k2)`. Distances are then
computed for every hit; the four Table 2 witnesses can be found in the CSV.

## What is verified

`verify.py` checks the following items from scratch.

1. The length-8 marker construction over `F4`, including `ind_theta(f_s)=8`,
   divisibility `ind_theta(f_s) | n`, cyclic optimum `4`, matched gain `2`,
   m-step linearization, weight distributions, and the `[[8,0]]` CSS condition.
2. Divisibility `ind_theta(f) | n` for the enumerated right divisors in its
   selected finite-field search ranges. The separate Table 2 search asserts
   right-divisibility and index divisibility for every retained nested quotient.
3. Selected matched-gain landscape points over `F4`, `F8`, `F9`, `F16`, `F25`, and `F27`.
4. Four positive-logical-dimension strict-gain CSS codes with `d1,d2 >= 3`.
5. Actual bounded-distance syndrome decoders for the four positive-rate codes,
   including exhaustive recovery of every weight-1 bit and phase error.

## Field and ring conventions

Elements of `GF(p^r)` are encoded as integers `0..p^r-1` using base-`p`
coefficients of the residue-class polynomial. Polynomials are coefficient lists
in increasing degree order.

The Ore multiplication rule is `x a = theta(a) x`. In all manuscript
computations, `theta` is the Frobenius automorphism `a -> a^p`.

## Table 2 positive-rate codes

| Row | Field | `(n,k1,k2)` | CSS code | `ind` | cyclic optimum | gain | `d1` | `d2` |
|-----|-------|--------------|----------|-------|----------------|------|------|------|
| 1 (A) | `F4` | `(16,10,11)` | `[[16,4]]` | 2 | 1 | 2 | 4 | 4 |
| 2 (B) | `F8` | `(9,6,7)` | `[[9,3]]` | 3 | 1 | 3 | 4 | 3 |
| 3 (A2) | `F4` | `(16,9,11)` | `[[16,2]]` | 4 | 2 | 2 | 4 | 4 |
| 4 (B2) | `F8` | `(12,7,8)` | `[[12,2]]` | 3 | 1 | 3 | 5 | 4 |

## Search ranges

For Figure 4, the exhaustive searches use the following exact finite boxes.
The figure takes the maximum over the listed degrees at each listed length;
it makes no claim about unsearched lengths or degrees.

| Field | Plotted lengths `n` | Degree cap `s_max` |
|-------|---------------------|--------------------|
| `F4` | `2,4,6,8,10,12,14,16` | 5 |
| `F8` | `3,6,9,12` | 5 |
| `F9` | `2,4,6,8,10` | 4 |
| `F16` | `4,8,12` | 2 |
| `F25` | `2,4` | 2 |
| `F27` | `3,6` | 2 |

At each `n`, the tested degrees are `s=1,...,min(s_max,n-1)`.
`reproduce_search.py` enumerates all `(q-1)q^(s-1)` monic degree-`s`
polynomials with nonzero constant term. For the skew numerator it retains a
polynomial only when it right-divides `x^n-1` in `GF(q)[x;theta]`; for the
cyclic denominator it uses ordinary polynomial divisibility. It then takes
the largest right-remainder return time or ordinary polynomial order,
respectively. The Frobenius is `theta(a)=a^p` and the field moduli are those
in `verify.py`: `F4=(1,1,1)`, `F8=(1,1,0,1)`, `F9=(1,0,1)`,
`F16=(1,1,0,0,1)`, `F25=(2,0,1)`, `F27=(1,2,0,1)`, in increasing
coefficient order over their prime fields.

For Table 2, the code fixes each `(n,k1,k2)` and enumerates monic outer
generators `g2` of degree `n-k2` that right-divide `x^n-1`, then every monic
quotient `fs` of degree `k2-k1` with nonzero constant term. It retains the
pair only if `g1=fs*g2` also right-divides `x^n-1`. The cyclic denominator is
computed independently by enumerating all nested ordinary cyclic generator
pairs at the same dimensions, without dual-containment or distance filters.

## Archival

Keep the cited commit or release identical to the data and code examined by
reviewers. If a separate archive is cited, confirm that its files match this
version and can be downloaded by the intended readers.
