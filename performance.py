from ffield import GFTable
from ore import Ore, skew_gen_matrix
from decoder import trial_success_rate


def cases():
    F4 = GFTable(2, [1, 1, 1])
    F8 = GFTable(2, [1, 1, 0, 1])
    return [
        ("A", F4, 1, 16, [1, 1, 2, 3, 2, 2, 1], [2, 3, 0, 1, 1, 1], 4, 4),
        ("A2", F4, 1, 16, [2, 0, 3, 1, 2, 0, 3, 1], [1, 2, 0, 2, 0, 1], 4, 4),
        ("B", F8, 1, 9, [2, 1, 5, 1], [3, 6, 1], 4, 3),
        ("B2", F8, 1, 12, [2, 1, 3, 6, 7, 1], [3, 6, 4, 0, 1], 5, 4),
    ]


def main():
    trials = 2000
    probs = [0.005, 0.01, 0.02, 0.05]
    print("case field n p bit_success phase_success combined_success trials")
    for name, F, tp, n, g1, g2, d1, d2 in cases():
        O = Ore(F, tp)
        G1 = skew_gen_matrix(O, g1, n)
        G2 = skew_gen_matrix(O, g2, n)
        t_phase = (d1 - 1) // 2
        t_bit = (d2 - 1) // 2
        for p in probs:
            bit = trial_success_rate(F, G2, t_bit, p, trials, 1000 + n)
            phase = trial_success_rate(F, G1, t_phase, p, trials, 2000 + n)
            print(f"{name} F{F.q} {n} {p:.3f} {bit:.4f} {phase:.4f} {bit * phase:.4f} {trials}")


if __name__ == "__main__":
    main()
