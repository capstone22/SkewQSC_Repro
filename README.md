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
(MSIT) (RS-2026-25528270, DevOps-based Post Quantum Cryptography Switching
Automation Technology and Open Platform Technology, 50%) and by the InnoCORE
program of the Ministry of Science and ICT (26-InnoCORE-01, 50%).

Access policy:

This repository is intended to support peer review and publication of the
manuscript. During peer review, access may be restricted to editors and
reviewers. Upon acceptance, the repository will be made publicly available and
archived with a Zenodo DOI.

The code uses only the Python standard library.

## Files

| File | Role |
|------|------|
| `ffield.py` | Finite fields, row reduction, dual spaces, and brute-force minimum distance. |
| `ore.py` | Ore-polynomial arithmetic, right division, skew index, ordinary order, and skew-cyclic generator matrices. |
| `mindist.py` | Faster minimum-distance computation from parity-check column dependencies. |
| `decoder.py` | Bounded-distance q-ary syndrome decoder. |
| `performance.py` | Monte Carlo bit/phase decoding performance for the flagship positive-rate codes. |
| `search2.py` | Exhaustive search for positive-rate strict-gain skew CSS families. |
| `verify.py` | Reproduction driver for the manuscript tables, figures, examples, and decoder checks. |
| `hits.pkl` | Pickled search hits from the positive-rate strict-gain search. |

## Quick start

Run the full verification suite:

```bash
python3 verify.py
```

The expected final line is:

```text
RESULT:  42/42 checks passed   ALL CHECKS PASSED
```

Run the performance simulation:

```bash
python3 performance.py
```

The script reports bit-decoder, phase-decoder, and combined success estimates
for four positive-rate flagship CSS codes at several independent symbol-error
probabilities.

## What is verified

`verify.py` checks the following items from scratch.

1. The length-8 marker construction over `F4`, including `ind_theta(f_s)=8`,
   divisibility `ind_theta(f_s) | n`, cyclic optimum `4`, matched gain `2`,
   m-step linearization, weight distributions, and the `[[8,0]]` CSS condition.
2. Divisibility `ind_theta(f_s) | n` for every realizable quotient in the listed
   finite-field search ranges.
3. The matched-gain landscape over `F4`, `F8`, `F9`, `F16`, `F25`, and `F27`.
4. Four positive-logical-dimension strict-gain CSS codes with `d1,d2 >= 3`.
5. Actual bounded-distance syndrome decoders for the four positive-rate codes,
   including exhaustive recovery of every weight-1 bit and phase error.

## Field and ring conventions

Elements of `GF(p^r)` are encoded as integers `0..p^r-1` using base-`p`
coefficients of the residue-class polynomial. Polynomials are coefficient lists
in increasing degree order.

The Ore multiplication rule is `x a = theta(a) x`. In all manuscript
computations, `theta` is the Frobenius automorphism `a -> a^p`.

## Flagship positive-rate codes

| Code | Field | CSS code | `ind` | cyclic optimum | gain | `d1` | `d2` |
|------|-------|----------|-------|----------------|------|------|------|
| A | `F4` | `[[16,4]]` | 2 | 1 | 2 | 4 | 4 |
| A2 | `F4` | `[[16,2]]` | 4 | 2 | 2 | 4 | 4 |
| B | `F8` | `[[9,3]]` | 3 | 1 | 3 | 4 | 3 |
| B2 | `F8` | `[[12,2]]` | 3 | 1 | 3 | 5 | 4 |

## Search ranges

The exhaustive searches enumerate monic right divisors of `x^n-1` for lengths
`n` divisible by the automorphism order.

| Field | Range |
|-------|-------|
| `F4` | `n <= 16`, synchronization quotient degree `s <= 5` |
| `F8` | `n <= 12`, `s <= 5` |
| `F9` | `n <= 10`, `s <= 4` |
| `F16`, `F27` | landscape searches up to `n <= 12` and `n <= 9` respectively |
| `F25` | landscape searches up to `n <= 10` |

## Archival

For publication, archive this repository using GitHub and Zenodo, then cite the
Zenodo DOI in the manuscript data and code availability statements.
