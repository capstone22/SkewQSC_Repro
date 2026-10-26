"""Recompute the search-backed rows of Tables 2-3 and Figure 4."""

import csv
import itertools
from pathlib import Path

from search2 import cyclic_opt
from ore import Ore, monic_divisors_right, xn_minus_1
from verify import (
    F4, F8, F9, F16, F25, F27, cyclic_optimum, max_skew_index,
)


# Every n and s in these finite search boxes is enumerated, not sampled.
FIELDS = {
    "F4": (F4, 2, (2, 4, 6, 8, 10, 12, 14, 16), 5),
    "F8": (F8, 3, (3, 6, 9, 12), 5),
    "F9": (F9, 2, (2, 4, 6, 8, 10), 4),
    "F16": (F16, 4, (4, 8, 12), 2),
    "F25": (F25, 2, (2, 4), 2),
    "F27": (F27, 3, (3, 6), 2),
}

TABLE2 = (
    ("F4", 16, 10, 11, 2, 1),
    ("F8", 9, 6, 7, 3, 1),
    ("F4", 16, 9, 11, 4, 2),
    ("F8", 12, 7, 8, 3, 1),
)

TABLE3 = (
    ("F4", 8, 3, 8, 4), ("F4", 8, 1, 2, 1),
    ("F4", 6, 2, 6, 6), ("F4", 16, 5, 16, 8),
    ("F8", 3, 1, 3, 1), ("F8", 12, 3, 12, 4),
    ("F16", 4, 1, 4, 1), ("F16", 8, 1, 4, 1),
    ("F16", 12, 2, 12, 6), ("F27", 3, 1, 3, 1),
    ("F27", 6, 1, 6, 2), ("F9", 6, 2, 6, 6),
    ("F9", 8, 2, 8, 8), ("F25", 4, 2, 4, 4),
)

FIGURE4 = {
    "F4": (2, 2, 1, 2, 2, 1, 2, 2),
    "F8": (3, 3, 3, 3),
    "F9": (1, 1, 1, 1, 1),
    "F16": (4, 4, 2),
    "F25": (1, 1),
    "F27": (3, 3),
}


def skew_nested_opt(field, n, k1, k2):
    ring = Ore(field, 1)
    target = xn_minus_1(field, n)
    degree2 = n - k2
    quotient_degree = k2 - k1
    best = 0
    for g2 in monic_divisors_right(ring, target, degree2):
        for constant in range(1, field.q):
            for tail in itertools.product(range(field.q), repeat=quotient_degree - 1):
                fs = [constant, *tail, field.ONE]
                g1 = ring.mul(fs, g2)
                if ring.right_divides(g1, target):
                    index = ring.skew_index(fs)
                    assert ring.right_divides(fs, target) and n % index == 0
                    best = max(best, index)
    return best


def main():
    for name, n, k1, k2, expected_skew, expected_cyclic in TABLE2:
        field = FIELDS[name][0]
        skew = skew_nested_opt(field, n, k1, k2)
        cyclic = cyclic_opt(field, n, k1, k2)
        assert (skew, cyclic) == (expected_skew, expected_cyclic), (
            name, n, k1, k2, skew, cyclic, expected_skew, expected_cyclic
        )
        print(f"Table 2: {name} ({n},{k1},{k2}) maxima = ({skew}, {cyclic})", flush=True)

    results = {}
    rows = []
    figure_rows = []
    for name, (field, m, lengths, cap) in FIELDS.items():
        for n, expected_gain in zip(lengths, FIGURE4[name]):
            best = 0
            for s in range(1, min(cap, n - 1) + 1):
                skew = max_skew_index(field, 1, n, s)
                cyclic = cyclic_optimum(field, n, s)
                gain = skew / cyclic if cyclic else 0
                results[name, n, s] = (skew, cyclic)
                rows.append((name, m, n, s, skew, cyclic, gain))
                best = max(best, gain)
            assert abs(best - expected_gain) < 1e-12, (name, n, best, expected_gain)
            figure_rows.append((name, m, n, best))
            print(f"Figure 4: {name} n={n} max gain={best:g}", flush=True)

    for name, n, s, expected_skew, expected_cyclic in TABLE3:
        actual = results[name, n, s]
        assert actual == (expected_skew, expected_cyclic), (
            name, n, s, actual, (expected_skew, expected_cyclic)
        )
        print(f"Table 3: {name} n={n} s={s} maxima={actual}", flush=True)

    path = Path(__file__).with_name("landscape_data.csv")
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("field", "m", "n", "s", "max_skew_index", "max_cyclic_order", "gain"))
        writer.writerows(rows)
    print(f"Wrote {len(rows)} exhaustive degree-level results to {path.name}")

    figure_path = Path(__file__).with_name("figure4_points.csv")
    with figure_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("field", "m", "n", "max_gain"))
        writer.writerows(figure_rows)
    print(f"Wrote {len(figure_rows)} plotted points to {figure_path.name}")


if __name__ == "__main__":
    main()
