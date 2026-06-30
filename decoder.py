import itertools
import random
from ffield import dual_basis


def zero_vec(F, n):
    return [F.ZERO] * n


def add_vec(F, a, b):
    return [F.add(x, y) for x, y in zip(a, b)]


def syndrome(F, H, v):
    out = []
    for row in H:
        s = F.ZERO
        for a, b in zip(row, v):
            if a != F.ZERO and b != F.ZERO:
                s = F.add(s, F.mul(a, b))
        out.append(s)
    return tuple(out)


def error_patterns(F, n, t):
    yield zero_vec(F, n)
    values = [a for a in range(F.q) if a != F.ZERO]
    for w in range(1, t + 1):
        for positions in itertools.combinations(range(n), w):
            for coeffs in itertools.product(values, repeat=w):
                e = zero_vec(F, n)
                for i, c in zip(positions, coeffs):
                    e[i] = c
                yield e


def syndrome_table(F, G, t):
    H = dual_basis(F, G)
    table = {}
    ambiguous = set()
    for e in error_patterns(F, len(G[0]), t):
        s = syndrome(F, H, e)
        if s in table and table[s] != e:
            ambiguous.add(s)
        else:
            table[s] = e
    for s in ambiguous:
        table.pop(s, None)
    return H, table


def decode(F, H, table, received):
    s = syndrome(F, H, received)
    e = table.get(s)
    if e is None:
        return None
    return add_vec(F, received, e)


def random_error(F, n, probability, rng):
    values = [a for a in range(F.q) if a != F.ZERO]
    e = zero_vec(F, n)
    for i in range(n):
        if rng.random() < probability:
            e[i] = rng.choice(values)
    return e


def trial_success_rate(F, G, t, probability, trials, seed):
    rng = random.Random(seed)
    H, table = syndrome_table(F, G, t)
    ok = 0
    n = len(G[0])
    sent = zero_vec(F, n)
    for _ in range(trials):
        e = random_error(F, n, probability, rng)
        received = add_vec(F, sent, e)
        decoded = decode(F, H, table, received)
        if decoded == sent:
            ok += 1
    return ok / trials
