"""
Validates weighted and hierarchical aggregation, closure in U~, order
independence and confidence-interval coverage.

Runs under pytest or directly:  python3 tests/test_weighting.py
"""

import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import fusion  # noqa: E402
from fusion import InvalidOpinion, aggregate, decide, pair_weights  # noqa: E402
from qif_xnor import xnor_analytic  # noqa: E402


def _truth(ops, w=None):
    pairs = list(itertools.combinations(range(len(ops)), 2))
    if w is None:
        w = {p: 1 / len(pairs) for p in pairs}
    mu = sum(w[p] * xnor_analytic(ops[p[0]], ops[p[1]])[0] for p in pairs)
    nu = sum(w[p] * xnor_analytic(ops[p[0]], ops[p[1]])[1] for p in pairs)
    return mu, nu


# =====================================================
# Weights
# =====================================================

def test_weights_sum_to_one():
    for mode, kw in (
        ("uniform", {}),
        ("weighted", {"weights": [1, 2, 3, 4]}),
        ("hierarchical", {"groups": ["a", "a", "b", "b"]}),
    ):
        w = pair_weights(4, combine=mode, **kw)
        assert abs(sum(w.values()) - 1.0) < 1e-12, mode
        assert all(v >= 0 for v in w.values())
        assert len(w) == 6


def test_uniform_is_weighted_with_equal_weights():
    a = pair_weights(5, "uniform")
    b = pair_weights(5, "weighted", weights=[1.0] * 5)
    for p in a:
        assert abs(a[p] - b[p]) < 1e-12


def test_hierarchical_equalises_institutions():
    """Each institution carries total source weight 1/K, whatever its size."""
    groups = ["A"] * 6 + ["B"] * 2
    w = pair_weights(8, "hierarchical", groups=groups)
    intra_a = sum(v for (i, j), v in w.items() if groups[i] == groups[j] == "A")
    intra_b = sum(v for (i, j), v in w.items() if groups[i] == groups[j] == "B")
    cross = sum(v for (i, j), v in w.items() if groups[i] != groups[j])
    # source weight 1/(K n_k) => intra_k = C(n_k,2)/(K n_k)^2, cross = 2/(K^2)*...
    assert abs(intra_a + intra_b + cross - 1.0) < 1e-12
    # the large hospital does not dominate: its intra weight is below the cross weight
    assert cross > intra_a


def test_bad_weight_inputs_are_rejected():
    for kw in (
        {"combine": "weighted"},                        # no weights
        {"combine": "weighted", "weights": [1, 2]},      # tamanho errado
        {"combine": "weighted", "weights": [1, -1, 1]},  # negativo
        {"combine": "hierarchical"},                     # no groups
        {"combine": "nope"},
    ):
        try:
            pair_weights(3, **kw)
        except InvalidOpinion:
            continue
        raise AssertionError(f"accepted invalid input: {kw}")


# =====================================================
# Fechamento em U~ e independencia de ordem
# =====================================================

def test_closed_in_U_for_every_scheme():
    ops = [(0.9, 0.05), (0.2, 0.7), (0.5, 0.4), (0.35, 0.35)]
    for mode, kw in (
        ("uniform", {}),
        ("weighted", {"weights": [3, 1, 4, 1]}),
        ("hierarchical", {"groups": ["A", "A", "B", "C"]}),
        ("tnorm", {}),
    ):
        c, _, _ = aggregate(ops, exact=True, combine=mode, **kw)
        assert c[0] >= -1e-12 and c[1] >= -1e-12, mode
        assert c[0] + c[1] <= 1.0 + 1e-9, (mode, c)


def test_order_independent_under_permutation():
    """The weights depend only on the unordered pair, so order is irrelevant."""
    ops = [(0.9, 0.05), (0.2, 0.7), (0.6, 0.3), (0.45, 0.4)]
    wts = [3.0, 1.0, 4.0, 1.0]
    grp = ["A", "A", "B", "C"]
    base = aggregate(ops, exact=True, combine="weighted", weights=wts)[0]
    baseh = aggregate(ops, exact=True, combine="hierarchical", groups=grp)[0]

    for perm in itertools.permutations(range(4)):
        o = [ops[i] for i in perm]
        got = aggregate(o, exact=True, combine="weighted",
                        weights=[wts[i] for i in perm])[0]
        assert abs(got[0] - base[0]) < 1e-12 and abs(got[1] - base[1]) < 1e-12
        goth = aggregate(o, exact=True, combine="hierarchical",
                         groups=[grp[i] for i in perm])[0]
        assert abs(goth[0] - baseh[0]) < 1e-12 and abs(goth[1] - baseh[1]) < 1e-12


def test_large_institution_no_longer_dominates():
    """30 sources from one hospital against 2 from another, in direct conflict."""
    ops = [(0.88, 0.07)] * 30 + [(0.08, 0.87)] * 2
    groups = ["A"] * 30 + ["B"] * 2

    flat = decide(ops, exact=True, combine="uniform")
    tier = decide(ops, exact=True, combine="hierarchical", groups=groups)

    assert flat["decision"] == "TREAT"          # the large hospital decides
    assert tier["decision"] == "REQUEST_EXAMS"     # the conflict surfaces
    assert tier["consensus"]["nu"] > flat["consensus"]["nu"]


def test_breakdown_separates_intra_from_inter():
    ops = [(0.88, 0.07)] * 3 + [(0.08, 0.87)] * 3
    groups = ["A"] * 3 + ["B"] * 3
    b = decide(ops, exact=True, combine="hierarchical", groups=groups)["breakdown"]
    assert b["intra"]["n_pairs"] == 6          # C(3,2) per institution
    assert b["inter"]["n_pairs"] == 9          # 3 x 3
    assert b["intra"]["mu"] > b["inter"]["mu"]  # internally coherent, externally in conflict


# =====================================================
# Intervalos de confianca
# =====================================================

def test_exact_mode_has_zero_width_intervals():
    u = decide([(0.9, 0.05), (0.2, 0.7)], exact=True)["uncertainty"]
    assert u["mu_ci95"] == 0.0 and u["nu_ci95"] == 0.0 and u["pi_ci95"] == 0.0


def test_interval_shrinks_as_one_over_sqrt_shots():
    ops = [(0.85, 0.1), (0.7, 0.2), (0.4, 0.5)]
    a = decide(ops, shots=4000, seed=1)["uncertainty"]["mu_ci95"]
    b = decide(ops, shots=64000, seed=1)["uncertainty"]["mu_ci95"]
    ratio = a / b
    assert 3.0 < ratio < 5.0, ratio      # expected 4x for 16x more shots


def test_interval_coverage_is_nominal():
    """
    The 95% interval must contain the true value in ~95% of the repetitions.
    With 120 repetitions the Monte Carlo error is ~2 pp, so 88-100% is
    accepted.
    """
    ops = [(0.90, 0.05), (0.85, 0.10), (0.20, 0.70), (0.55, 0.35)]
    tm, tn = _truth(ops)
    tp = 1 - tm - tn
    hits = [0, 0, 0]
    trials = 120
    for t in range(trials):
        r = decide(ops, shots=3000, seed=9000 + t)
        c, u = r["consensus"], r["uncertainty"]
        hits[0] += abs(c["mu"] - tm) <= u["mu_ci95"]
        hits[1] += abs(c["nu"] - tn) <= u["nu_ci95"]
        hits[2] += abs(c["pi"] - tp) <= u["pi_ci95"]
    for name, h in zip(("mu", "nu", "pi"), hits):
        assert 0.88 <= h / trials <= 1.0, (name, h / trials)


if __name__ == "__main__":
    fns = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_")]
    failed = 0
    for name, fn in fns:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {name}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} tests passed")
    sys.exit(1 if failed else 0)
