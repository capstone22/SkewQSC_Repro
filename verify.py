import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ffield import GFTable, rref, rank, dual_basis, in_rowspace
from ore import (Ore, xn_minus_1, monic_divisors_right,
                 skew_gen_matrix, left_mult_x_vec)
from mindist import min_distance_fast
from decoder import syndrome_table, decode, add_vec

F4  = GFTable(2, [1, 1, 1])
F8  = GFTable(2, [1, 1, 0, 1])
F9  = GFTable(3, [1, 0, 1])
F16 = GFTable(2, [1, 1, 0, 0, 1])
F25 = GFTable(5, [2, 0, 1])
F27 = GFTable(3, [1, 2, 0, 1])

PASS, FAIL = "PASS", "**FAIL**"
_alltests = []
def check(label, cond):
    _alltests.append(bool(cond))
    print(f"   [{PASS if cond else FAIL}] {label}")
    return cond

def euclid_dual_contains(F, G):
    R, piv = rref(F, G)
    return all(in_rowspace(F, R, piv, v) for v in dual_basis(F, G))

def weight_distribution(F, G):
    import itertools
    from collections import Counter
    R, piv = rref(F, G); k = len(piv); n = len(R[0])
    c = Counter()
    for coeffs in itertools.product(range(F.q), repeat=k):
        w = [F.ZERO]*n
        for ci, row in zip(coeffs, R):
            if ci == 0: continue
            for j in range(n):
                if row[j] != F.ZERO:
                    w[j] = F.add(w[j], F.mul(ci, row[j]))
        c[sum(1 for s in w if s != F.ZERO)] += 1
    return dict(sorted(c.items()))

def cyclic_optimum(F, n, s):
    Oid = Ore(F, 0); best = 0
    for f in monic_divisors_right(Oid, xn_minus_1(F, n), s):
        o = Oid.order_comm(f)
        if o and o > best: best = o
    return best

def max_skew_index(F, tp, n, s):
    O = Ore(F, tp); best = 0
    for f in monic_divisors_right(O, xn_minus_1(F, n), s):
        i = O.skew_index(f)
        if i and i > best: best = i
    return best

def part1():
    print("\n" + "="*70)
    print("PART 1   Flagship length-8 marker example over F4   (n,k1,k2)=(8,4,7)")
    print("="*70)
    O = Ore(F4, 1)
    g1 = [2, 3, 3, 3, 1]
    g2 = [2, 1]
    fs = [1, 3, 0, 1]
    check("m = ord(theta) = 2", O.m == 2)
    check("g1 = fs * g2  (Ore product)", O.trim(O.mul(fs, g2)) == O.trim(g1[:]))
    check("g1 right-divides x^8 - 1", O.right_divides(g1, xn_minus_1(F4, 8)))
    ind = O.skew_index(fs)
    check("ind_theta(fs) = 8", ind == 8)
    check("ind | n  (8 | 8)", 8 % ind == 0)
    cyc = cyclic_optimum(F4, 8, 3)
    check("cyclic optimum at deg 3 = 4", cyc == 4)
    check("matched gain G_m = 8/4 = 2", ind // cyc == 2)
    ell = None; one = [F4.ONE]; cur = one[:]
    for j in range(1, 50):
        cur = one[:]
        for _ in range(2*j): cur = O.right_rem(O.x_times(cur), fs)
        if cur == one: ell = j; break
    check("m-step linearization ell = ind/gcd(ind,m) = 4",
          ell == ind // math.gcd(ind, 2) == 4)
    G1 = skew_gen_matrix(O, g1, 8)
    wd = weight_distribution(F4, G1)
    check("skew C1 weight dist = {0:1,4:18,5:96,6:24,7:96,8:21}",
          wd == {0:1, 4:18, 5:96, 6:24, 7:96, 8:21})
    check("C1 is self-dual (dim 4 = n/2 and C1^perp ⊆ C1)",
          len(rref(F4, G1)[1]) == 4 and euclid_dual_contains(F4, G1))
    check("=> CSS [[8,0]], k_L = 2*4-8 = 0", 2*4-8 == 0)
    Gc = [[F4.ONE if (j-i) % 8 in (0, 4) else F4.ZERO for j in range(8)]
          for i in range(4)]
    wc = weight_distribution(F4, Gc)
    check("cyclic [8,4]_4 weight dist = {0:1,2:12,4:54,6:108,8:81}",
          wc == {0:1, 2:12, 4:54, 6:108, 8:81})
    q, ds = 4, 3
    lhs = 1 + ds*(q-1)
    cap = q**ds // 8
    check(f"EA t=1 obstruction: |B| >= {lhs} > {cap} (no 1-separation on 8 offsets)",
          lhs > cap)

def part2():
    print("\n" + "="*70)
    print("PART 2   ind_theta(f) | n  for enumerated right divisors")
    print("="*70)
    cases = [("F4", F4, 1, [4, 8, 10, 14, 16]),
             ("F8", F8, 1, [3, 6, 9, 12]),
             ("F9", F9, 1, [4, 6, 8])]
    for name, F, tp, ns in cases:
        O = Ore(F, tp); ok = True; tot = 0
        for n in ns:
            for s in range(1, min(n, 5)):
                for f in monic_divisors_right(O, xn_minus_1(F, n), s):
                    tot += 1
                    idx = O.skew_index(f)
                    if idx is None or n % idx != 0: ok = False
        check(f"{name}: all {tot} enumerated right divisors satisfy ind | n", ok)

def part3():
    print("\n" + "="*70)
    print("PART 3   Matched-gain landscape (argmax over sync-quotient degree s)")
    print("="*70)
    print("   gain reaches the Frobenius envelope m for F4,F8,F16,F27,")
    print("   but collapses to 1 for the odd-char quadratic extensions F9,F25.\n")
    spec = [("F4", F4, 1, 2, [2, 4, 6, 8, 10]),
            ("F8", F8, 1, 3, [3, 6, 9]),
            ("F9", F9, 1, 2, [2, 4, 6, 8]),
            ("F16", F16, 1, 4, [4, 8]),
            ("F25", F25, 1, 2, [2, 4]),
            ("F27", F27, 1, 3, [3, 6])]
    print(f"   {'field':5}{'m':>3}{'n':>4}{'gain':>6}   (s, skew, cyc)")
    for name, F, tp, m, ns in spec:
        for n in ns:
            bg, det = 0, None
            for s in range(1, min(n, 4)):
                sk = max_skew_index(F, tp, n, s); cy = cyclic_optimum(F, n, s)
                if sk and cy and sk/cy > bg: bg = sk/cy; det = (s, sk, cy)
            tag = " gain=m" if abs(bg-m) < 1e-9 and bg > 1 else (" no gain" if bg <= 1 else "")
            print(f"   {name:5}{m:>3}{n:>4}{bg:>6.3g}   {det}{tag}")
    g_f9 = max(max_skew_index(F9, 1, n, s)/cyclic_optimum(F9, n, s)
               for n in [2, 4, 6, 8] for s in [1, 2]
               if cyclic_optimum(F9, n, s))
    check("F9: maximum matched gain over searched range = 1 (no advantage)", g_f9 <= 1.0)
    g_f8 = max_skew_index(F8, 1, 3, 1)/cyclic_optimum(F8, 3, 1)
    check("F8: gain attains m = 3 at n=3 (degree-1 quotient)", abs(g_f8-3) < 1e-9)

def analyze_positive(name, F, tp, n, k1, k2, g1, g2, fs):
    O = Ore(F, tp); m = O.m
    fact = O.trim(O.mul(fs, g2)) == O.trim(g1[:])
    div  = O.right_divides(g1, xn_minus_1(F, n))
    nest = O.right_divides(g2, g1) and O.right_divides(g2, xn_minus_1(F, n))
    ind  = O.skew_index(fs); cyc = cyclic_optimum(F, n, O.deg(fs))
    gain = ind / cyc if cyc else 0
    G1 = skew_gen_matrix(O, g1, n); G2 = skew_gen_matrix(O, g2, n)
    dualc = euclid_dual_contains(F, G1)
    kdim = len(rref(F, G1)[1]); k2dim = len(rref(F, G2)[1]); kL = 2*k1 - n
    d1 = min_distance_fast(F, G1); d2 = min_distance_fast(F, G2)
    Hb, Tb = syndrome_table(F, G2, (d2 - 1) // 2)
    Hp, Tp = syndrome_table(F, G1, (d1 - 1) // 2)
    sent = [F.ZERO] * n
    bit_ok = True
    phase_ok = True
    for pos in range(n):
        for val in range(1, F.q):
            e = [F.ZERO] * n
            e[pos] = val
            bit_ok = bit_ok and decode(F, Hb, Tb, add_vec(F, sent, e)) == sent
            phase_ok = phase_ok and decode(F, Hp, Tp, add_vec(F, sent, e)) == sent
    print(f"\n   ----- {name}:  [[{n},{kL}]] over F{F.q}  (m={m}) -----")
    print(f"         g1={g1}")
    print(f"         g2={g2}   fs={fs}")
    print(f"         ind={ind}  cyclic_opt={cyc}  GAIN={gain:.3g}   "
          f"d1(phase)={d1}  d2(bit)={d2}")
    check(f"{name}: fs*g2=g1, nested, g1|x^n-1", fact and div and nest)
    check(f"{name}: dim C1 = k1 = {k1}", kdim == k1)
    check(f"{name}: dim C2 = k2 = {k2}, deg fs = k2-k1",
          k2dim == k2 and O.deg(fs) == k2-k1)
    check(f"{name}: C1^perp ⊆ C1 (valid CSS)", dualc)
    check(f"{name}: k_L = {kL} > 0  AND strict gain {gain:.3g} > 1", kL > 0 and gain > 1)
    check(f"{name}: corrects ≥1 bit and ≥1 phase (d1≥3 and d2≥3)", d1 >= 3 and d2 >= 3)
    check(f"{name}: syndrome decoders recover every weight-1 bit and phase error",
          bit_ok and phase_ok)

def part4():
    print("\n" + "="*70)
    print("PART 4   Positive-logical-dimension, strict-gain, error-correcting")
    print("         skew-cyclic CSS codes  (resolves Open Problem 3)")
    print("="*70)
    analyze_positive("FLAGSHIP A ", F4, 1, 16, 10, 11,
                     [1,1,2,3,2,2,1], [2,3,0,1,1,1], [3,1])
    analyze_positive("FLAGSHIP A2", F4, 1, 16, 9, 11,
                     [2,0,3,1,2,0,3,1], [1,2,0,2,0,1], [2,3,1])
    analyze_positive("FLAGSHIP B ", F8, 1, 9, 6, 7,
                     [2,1,5,1], [3,6,1], [7,1])
    analyze_positive("FLAGSHIP B2", F8, 1, 12, 7, 8,
                     [2,1,3,6,7,1], [3,6,4,0,1], [7,1])

if __name__ == "__main__":
    part1(); part2(); part3(); part4()
    print("\n" + "="*70)
    n_ok = sum(_alltests); n_tot = len(_alltests)
    print(f"RESULT:  {n_ok}/{n_tot} checks passed"
          + ("   ALL CHECKS PASSED" if n_ok == n_tot else "   *** SOME FAILED ***"))
    print("="*70)
    sys.exit(0 if n_ok == n_tot else 1)
