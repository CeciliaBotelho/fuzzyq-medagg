"""Reproduces Tables III and V of the main paper (analytical, no quantum simulation).
Table V: baselines (IFWA, IFWG, mean pairwise Szmidt-Kacprzyk similarity) vs XNOR consensus.
Table III: multi-source stress scenarios."""
from itertools import combinations
from math import prod

def SP(a, b): return a + b - a * b
def xnor(x, y):
    (x1, x2), (y1, y2) = x, y
    return SP(x2, y1) * SP(x1, y2), SP(x1 * y2, x2 * y1)

def ifwa(xs, w=None):
    n = len(xs); w = w or [1 / n] * n
    return 1 - prod((1 - m) ** wi for (m, _), wi in zip(xs, w)), prod(v ** wi for (_, v), wi in zip(xs, w))

def ifwg(xs, w=None):
    n = len(xs); w = w or [1 / n] * n
    return prod(m ** wi for (m, _), wi in zip(xs, w)), 1 - prod((1 - v) ** wi for (_, v), wi in zip(xs, w))

def sk_similarity(x, y):
    """Szmidt-Kacprzyk: 1 - d_H, d_H = 1/2(|mu1-mu2|+|nu1-nu2|+|pi1-pi2|)."""
    (a, b), (c, d) = x, y
    return 1 - 0.5 * (abs(a - c) + abs(b - d) + abs((1 - a - b) - (1 - c - d)))

def mean_similarity(xs):
    p = list(combinations(xs, 2)); return sum(sk_similarity(a, b) for a, b in p) / len(p)

def consensus(xs, omega=None):
    n = len(xs); omega = omega or [1] * n
    pairs = list(combinations(range(n), 2))
    Z = sum(omega[i] * omega[j] for i, j in pairs)
    mu = sum(omega[i] * omega[j] / Z * xnor(xs[i], xs[j])[0] for i, j in pairs)
    nu = sum(omega[i] * omega[j] / Z * xnor(xs[i], xs[j])[1] for i, j in pairs)
    return mu, nu, 1 - mu - nu

def action(xs, omega=None):
    mu, nu, pi = consensus(xs, omega)
    s = sum(m for m, _ in xs) / len(xs) - sum(v for _, v in xs) / len(xs)
    sc = {"Treat": mu * (1 + s) / 2, "DoNotTreat": mu * (1 - s) / 2, "Exams": nu + pi}
    best = max(sc.values())
    winners = [k for k, v in sc.items() if abs(v - best) < 1e-12]
    return s, sc, ("Exams" if len(winners) > 1 else winners[0])  # ties resolved in favour of Exams

print("== Table V: baselines vs XNOR consensus (uniform weights) ==")
for name, xs in [("Concordant support", [(0.85, 0.10), (0.80, 0.02)]),
                 ("Concordant rejection", [(0.10, 0.85), (0.02, 0.80)]),
                 ("Symmetric conflict", [(0.85, 0.10), (0.10, 0.85)]),
                 ("High hesitation", [(0.30, 0.20), (0.25, 0.25)]),
                 ("One dissonant of three", [(0.85, 0.10), (0.80, 0.05), (0.10, 0.85)])]:
    a, g, sim, c = ifwa(xs), ifwg(xs), mean_similarity(xs), consensus(xs)
    print(f"{name:24s} IFWA=({a[0]:.3f},{a[1]:.3f}) IFWG=({g[0]:.3f},{g[1]:.3f}) "
          f"1-dH={sim:.3f} C=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}) {action(xs)[2]}")

print("\n== Table III: multi-source stress scenarios ==")
A1, A2, A3 = (0.85, 0.10), (0.80, 0.05), (0.82, 0.08)
A, D1, D2 = [A1, A2, A3], (0.10, 0.85), (0.45, 0.45)
for name, xs, om in [("A only", A, None),
                     ("A + D1, uniform", A + [D1], None),
                     ("A + D1 duplicated", A + [D1, D1], None),
                     ("A + A1 duplicated + D1", A + [A1, D1], None),
                     ("A + D1, institution-balanced", A + [D1], [1, 1, 1, 3]),
                     ("A + D2, uniform", A + [D2], None),
                     ("A + D2, institution-balanced", A + [D2], [1, 1, 1, 3]),
                     ("Near-copies of A1 + D1", [A1, (0.84, 0.11), (0.86, 0.09), D1], None)]:
    c = consensus(xs, om); s, _, a = action(xs, om)
    print(f"{name:32s} n={len(xs)} C=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}) s={s:+.3f} {a}")
