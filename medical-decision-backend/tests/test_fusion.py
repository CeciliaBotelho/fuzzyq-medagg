"""
Validates the Evidence Fusion Core: score normalisation, direction via the
Xu-Yager order, aggregation of n opinions and input validation.

Runs under pytest or directly:  python3 tests/test_fusion.py
"""

import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import fusion  # noqa: E402
from fusion import InvalidOpinion, aggregate, decide, mean_ifv, score  # noqa: E402
from qif_xnor import derive_seeds, run_xnor, xnor_analytic  # noqa: E402


def _grid(step=0.1):
    n = int(round(1 / step))
    vals = [i * step for i in range(n + 1)]
    return [(mu, nu) for mu, nu in itertools.product(vals, repeat=2) if mu + nu <= 1 + 1e-12]


# =====================================================
# Input validation
# =====================================================

def test_rejects_outside_U():
    for bad in [(0.7, 0.5), (1.1, 0.0), (-0.1, 0.5), (0.5, 1.2)]:
        try:
            decide([bad], exact=True)
        except InvalidOpinion:
            continue
        raise AssertionError(f"accepted an invalid opinion: {bad}")


def test_rejects_empty():
    try:
        decide([], exact=True)
    except InvalidOpinion:
        return
    raise AssertionError("accepted an empty list")


def test_rejects_outside_U_before_counting():
    """U~ validation runs before the n check, so an invalid n=1 also fails."""
    try:
        decide([(0.7, 0.5)], exact=True)
    except InvalidOpinion as e:
        assert "intuitionistic" in str(e)
        return
    raise AssertionError("accepted an opinion outside U~")


# =====================================================
# Normalisation
# =====================================================

def test_scores_sum_to_one():
    """The agreement/conflict/hesitation partition sums to 1 without renormalising."""
    for a, b in itertools.product(_grid(0.2), repeat=2):
        r = decide([a, b], exact=True)
        total = sum(r["scores"].values())
        assert abs(total - 1.0) < 1e-9, (a, b, total)
        assert all(v >= -1e-12 for v in r["scores"].values())


def test_requestExams_is_conflict_plus_hesitation():
    for a, b in itertools.product(_grid(0.25), repeat=2):
        r = decide([a, b], exact=True)
        c = r["consensus"]
        assert abs(r["scores"]["REQUEST_EXAMS"] - (c["nu"] + c["pi"])) < 1e-9


# =====================================================
# Agreement and direction semantics
# =====================================================

def test_xnor_is_direction_agnostic():
    """Agreeing on treat and agreeing on do-not-treat yield the same agreement."""
    pos = decide([(0.95, 0.02), (0.90, 0.05)], exact=True)
    neg = decide([(0.02, 0.95), (0.05, 0.90)], exact=True)
    assert abs(pos["consensus"]["mu"] - neg["consensus"]["mu"]) < 1e-9
    assert pos["decision"] == "TREAT"
    assert neg["decision"] == "DO_NOT_TREAT"


def test_direction_follows_xu_yager_score():
    for a, b in itertools.product(_grid(0.2), repeat=2):
        r = decide([a, b], exact=True)
        s = score(mean_ifv([a, b]))
        if r["decision"] == "TREAT":
            assert s > 0, (a, b, s)
        elif r["decision"] == "DO_NOT_TREAT":
            assert s < 0, (a, b, s)


def test_conflict_forces_exams():
    r = decide([(0.9, 0.05), (0.05, 0.9)], exact=True)
    assert r["decision"] == "REQUEST_EXAMS"
    assert r["consensus"]["nu"] > 0.7


def test_hesitation_forces_exams():
    """A case where the previous version zeroed pi and lost the signal."""
    r = decide([(0.3, 0.3), (0.3, 0.3)], exact=True)
    assert r["decision"] == "REQUEST_EXAMS"
    assert r["consensus"]["pi"] > 0.5


def test_symmetric_opinions_tie_to_exams():
    """score = 0 cannot decide between treat and do-not-treat."""
    r = decide([(0.5, 0.5), (0.5, 0.5)], exact=True)
    assert r["decision"] == "REQUEST_EXAMS"


# =====================================================
# Aggregation of n opinions
# =====================================================

def test_single_opinion_is_rejected():
    """
    n = 1 defines no agreement: there is no pair. The pipeline requires
    n >= 2, with no special branch. Detailed coverage in
    tests/test_pipeline.py.
    """
    for x in _grid(0.2)[:5]:
        try:
            aggregate([x], exact=True)
        except InvalidOpinion:
            continue
        raise AssertionError(f"accepted n = 1 with {x}")


def test_n_opinions_produce_all_pairs():
    """C(n,2) applications of boxplus_I, one per pair of opinions."""
    for n in range(2, 7):
        ops = [(0.8, 0.1)] * n
        _, steps, _ = aggregate(ops, exact=True)
        assert len(steps) == n * (n - 1) // 2
        seen = {(s["pair_i"], s["pair_j"]) for s in steps}
        assert len(seen) == len(steps)
        assert all(i < j for i, j in seen)
        r = decide(ops, exact=True)
        assert r["n_opinions"] == n


def test_tnorm_combine_stays_in_U():
    """O modo conservador tambem e fechado em U~."""
    ops = [(0.85, 0.10), (0.80, 0.15), (0.90, 0.05), (0.75, 0.20)]
    for n in range(2, len(ops) + 1):
        c, _, _ = aggregate(ops[:n], exact=True, combine="tnorm")
        assert c[0] + c[1] <= 1.0 + 1e-9, (n, c)


def test_more_agreeing_opinions_do_not_break_U():
    ops = [(0.85, 0.10), (0.80, 0.15), (0.90, 0.05), (0.75, 0.20), (0.88, 0.08)]
    for n in range(2, len(ops) + 1):
        r = decide(ops[:n], exact=True)
        c = r["consensus"]
        assert c["mu"] + c["nu"] <= 1.0 + 1e-9, (n, c)
        assert c["pi"] >= -1e-12


def test_order_independent():
    """Pairwise agreement does not depend on the order the sources are given."""
    ops = [(0.9, 0.05), (0.2, 0.7), (0.6, 0.3), (0.45, 0.4)]
    base = aggregate(ops, exact=True)[0]
    for perm in itertools.permutations(ops):
        got = aggregate(list(perm), exact=True)[0]
        assert abs(got[0] - base[0]) < 1e-12, (perm, got, base)
        assert abs(got[1] - base[1]) < 1e-12, (perm, got, base)


def test_direction_agnostic_for_any_n():
    """
    Regression test for the sequential-fold defect: agreement among n
    sources must not depend on the direction they point to.
    """
    pos = [(0.90, 0.05), (0.85, 0.10), (0.88, 0.07), (0.92, 0.04)]
    neg = [(mu_, nu_) for nu_, mu_ in pos]  # espelha em torno da diagonal
    for n in range(2, len(pos) + 1):
        cp = aggregate(pos[:n], exact=True)[0]
        cn = aggregate(neg[:n], exact=True)[0]
        assert abs(cp[0] - cn[0]) < 1e-12, (n, cp, cn)
        assert abs(cp[1] - cn[1]) < 1e-12, (n, cp, cn)


def test_sequential_fold_would_be_direction_dependent():
    """
    Documents why the sequential fold was abandoned: at n = 3 it separates
    positive from negative agreement, which must be identical.
    """
    def fold(ops):
        acc = ops[0]
        for nxt in ops[1:]:
            acc = xnor_analytic(acc, nxt)
        return acc

    pos = [(0.90, 0.05), (0.85, 0.10), (0.88, 0.07)]
    neg = [(nu_, mu_) for mu_, nu_ in pos]
    assert abs(fold(pos)[0] - fold(neg)[0]) > 0.4


def test_circuit_and_analytic_agree():
    """Aggregation through the simulator reproduces the analytical path."""
    ops = [(0.85, 0.10), (0.80, 0.15), (0.70, 0.20)]
    exact = decide(ops, exact=True)
    sampled = decide(ops, shots=200_000)
    assert exact["decision"] == sampled["decision"]
    for k in ("mu", "nu"):
        assert abs(exact["consensus"][k] - sampled["consensus"][k]) < 8e-3


# =====================================================
# Reprodutibilidade
# =====================================================

def test_seed_makes_run_reproducible():
    ops = [(0.85, 0.10), (0.80, 0.15), (0.90, 0.05)]
    runs = [decide(ops, shots=4000, seed=2026) for _ in range(3)]
    for r in runs[1:]:
        assert r["consensus"] == runs[0]["consensus"]
        assert r["scores"] == runs[0]["scores"]
        assert r["steps"] == runs[0]["steps"]


def test_different_seeds_give_different_samples():
    ops = [(0.85, 0.10), (0.80, 0.15)]
    a = decide(ops, shots=4000, seed=1)["consensus"]
    b = decide(ops, shots=4000, seed=2)["consensus"]
    assert a != b


def test_derived_seeds_decorrelate_identical_pairs():
    """
    Reusing a single seed across the C(n,2) circuits would make pairs with
    identical inputs return identical counts, cancelling the benefit of
    averaging.
    """
    identical = [(0.7, 0.2)] * 4          # all 6 pairs are identical
    _, steps, _ = aggregate(identical, shots=4000, seed=11)
    assert len({s["mu"] for s in steps}) > 1

    # e o mesmo seed base continua devolvendo a mesma sequencia
    _, again, _ = aggregate(identical, shots=4000, seed=11)
    assert [s["mu"] for s in steps] == [s["mu"] for s in again]


def test_seed_sequence_streams_are_distinct_and_in_range():
    """SeedSequence.spawn substitui o espalhamento multiplicativo manual."""
    seeds = {v for base in range(1, 6) for v in derive_seeds(base, 20)}
    assert len(seeds) == 100
    assert all(0 <= v < 2**31 - 1 for v in seeds)
    assert derive_seeds(7, 5) == derive_seeds(7, 5)


def test_seedless_runs_stay_stochastic():
    """Fixing a seed by default would make a sampling method look exact."""
    x, y = (0.7, 0.2), (0.65, 0.25)
    vals = {run_xnor(x, y, shots=4000)[0] for _ in range(6)}
    assert len(vals) > 1


def test_provenance_records_the_run():
    ops = [(0.85, 0.10), (0.80, 0.15)]
    p = decide(ops, shots=4000, seed=5)["provenance"]
    assert p["shots"] == 4000 and p["seed"] == 5
    assert p["qiskit"] and p["qiskit_aer"]
    assert decide(ops, exact=True)["provenance"]["seed"] is None


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
