"""
Validates that the pairwise pipeline is SINGLE and covers every n >= 2.

There must be no special branch for n = 2: that is merely the case
P = C(2,2) = 1, in which the aggregate reduces to the single
boxplus_I(x~_1, x~_2).

Runs under pytest or directly:  python3 tests/test_pipeline.py
"""

import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fusion import InvalidOpinion, aggregate, decide, pair_weights  # noqa: E402
from qif_xnor import derive_seeds, run_xnor, xnor_analytic  # noqa: E402

POOL = [(0.90, 0.05), (0.85, 0.10), (0.20, 0.70), (0.55, 0.35), (0.40, 0.50)]


# =====================================================
# n < 2 e rejeitado
# =====================================================

def test_rejects_zero_opinions():
    for call in (lambda: aggregate([], exact=True), lambda: decide([], exact=True)):
        try:
            call()
        except InvalidOpinion as e:
            assert "two" in str(e)
            continue
        raise AssertionError("accepted n = 0")


def test_rejects_a_single_opinion():
    """
    n = 1 defines no agreement: there is no pair. The code previously
    returned the opinion itself through a special branch; it is now
    rejected.
    """
    for call in (
        lambda: aggregate([(0.7, 0.2)], exact=True),
        lambda: decide([(0.7, 0.2)], exact=True),
    ):
        try:
            call()
        except InvalidOpinion as e:
            assert "two" in str(e) and "got 1" in str(e)
            continue
        raise AssertionError("accepted n = 1")


# =====================================================
# P = C(n,2) for each n
# =====================================================

def test_pair_count_is_n_choose_2():
    for n, expected in ((2, 1), (3, 3), (4, 6), (5, 10)):
        ops = POOL[:n]
        _, steps, _ = aggregate(ops, exact=True)
        assert len(steps) == expected, (n, len(steps))
        assert len(pair_weights(n, "uniform")) == expected

        seen = {(s["pair_i"], s["pair_j"]) for s in steps}
        assert len(seen) == expected                 # no repetition
        assert all(i < j for i, j in seen)            # unordered pairs only
        assert seen == {(i + 1, j + 1)
                        for i, j in itertools.combinations(range(n), 2)}


def test_same_pipeline_shape_for_every_n():
    """The same keys and the same weight invariant, from n = 2 to n = 5."""
    keys = None
    for n in (2, 3, 4, 5):
        _, steps, stats = aggregate(POOL[:n], exact=True)
        assert abs(sum(s["weight"] for s in steps) - 1.0) < 1e-12
        if keys is None:
            keys = sorted(steps[0])
        assert sorted(steps[0]) == keys                # estrutura identica
        assert stats["n_pairs"] == n * (n - 1) // 2


# =====================================================
# n = 2 reduz ao XNOR unico
# =====================================================

def test_single_pair_weight_is_exactly_one():
    w = pair_weights(2, "uniform")
    assert list(w) == [(0, 1)]
    assert w[(0, 1)] == 1.0                            # exact, not approximate


def test_uniform_weights_sum_to_one():
    for n in (2, 3, 4, 5, 8):
        w = pair_weights(n, "uniform")
        assert abs(sum(w.values()) - 1.0) < 1e-12
        assert all(abs(v - 1.0 / len(w)) < 1e-12 for v in w.values())


def test_n2_consensus_equals_analytic_xnor():
    """C~ must be bit-for-bit equal to xnor_analytic for n = 2."""
    for x, y in itertools.combinations(POOL, 2):
        c, steps, _ = aggregate([x, y], exact=True)
        assert c == xnor_analytic(x, y), (x, y, c)
        assert len(steps) == 1
        assert steps[0]["weight"] == 1.0
        assert (steps[0]["mu"], steps[0]["nu"]) == c


def test_n2_consensus_matches_the_circuit():
    """
    Through the sampling path, C~ must coincide with run_xnor on the same
    pair.

    The seed of the pair is DERIVED from the base seed, also when there is a
    single pair -- the pipeline makes no exception for n = 2. That is why the
    comparison uses derive_seeds(base, 1)[0] and not the base itself.
    """
    x, y = (0.90, 0.05), (0.85, 0.10)
    c, _, _ = aggregate([x, y], shots=200_000, seed=17)
    mu, nu, _ = run_xnor(x, y, shots=200_000, seed=derive_seeds(17, 1)[0])
    assert c == (mu, nu)

    # and the result is reproducible from the base seed
    again, _, _ = aggregate([x, y], shots=200_000, seed=17)
    assert again == c


def test_n2_pi_is_the_complement():
    for x, y in itertools.combinations(POOL, 2):
        r = decide([x, y], exact=True)
        mu12, nu12 = xnor_analytic(x, y)
        c = r["consensus"]
        assert c["mu"] == mu12 and c["nu"] == nu12
        assert abs(c["pi"] - (1.0 - mu12 - nu12)) < 1e-12


# =====================================================
# Permutation invariance for n >= 3
# =====================================================

def test_permutation_invariance_for_n_at_least_3():
    for n in (3, 4, 5):
        ops = POOL[:n]
        base = aggregate(ops, exact=True)[0]
        for perm in itertools.permutations(range(n)):
            got = aggregate([ops[i] for i in perm], exact=True)[0]
            assert abs(got[0] - base[0]) < 1e-12, (n, perm, got, base)
            assert abs(got[1] - base[1]) < 1e-12, (n, perm, got, base)


def test_permutation_invariance_holds_for_n2_trivially():
    """D_I2: boxplus_I is symmetric, so the order of the single pair is irrelevant."""
    for x, y in itertools.combinations(POOL, 2):
        assert aggregate([x, y], exact=True)[0] == aggregate([y, x], exact=True)[0]


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
