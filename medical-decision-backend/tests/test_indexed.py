"""
Validates the indexed register: it must converge to the pairwise mean.

Runs under pytest or directly:  python3 tests/test_indexed.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import indexed as ix  # noqa: E402
from fusion import aggregate  # noqa: E402

SHOTS = 120_000
TOL = 0.012


def test_converges_to_pairwise_mean():
    """P(mu=1 | valid branch) = mean of mu_boxplus over all pairs."""
    for ops in (
        [(0.90, 0.05), (0.85, 0.10)],
        [(0.90, 0.05), (0.85, 0.10), (0.20, 0.70)],
        [(0.90, 0.05), (0.85, 0.10), (0.20, 0.70), (0.55, 0.35)],
    ):
        mt, nt = ix.mean_pairwise_analytic(ops)
        mu, nu, _, _ = ix.run_indexed(ops, shots=SHOTS, seed=7)
        assert abs(mu - mt) < TOL, (ops, mu, mt)
        assert abs(nu - nt) < TOL, (ops, nu, nt)


def test_agrees_with_fusion_aggregate():
    """The indexed register and uniform aggregation compute the same quantity."""
    ops = [(0.85, 0.10), (0.80, 0.15), (0.70, 0.20)]
    exact, _, _ = aggregate(ops, exact=True, combine="uniform")
    mu, nu, _, _ = ix.run_indexed(ops, shots=SHOTS, seed=3)
    assert abs(mu - exact[0]) < TOL
    assert abs(nu - exact[1]) < TOL


def test_stays_in_U():
    ops = [(0.90, 0.05), (0.10, 0.85), (0.40, 0.40)]
    mu, nu, pi, _ = ix.run_indexed(ops, shots=SHOTS, seed=11)
    assert mu + nu <= 1.0 + TOL
    assert pi >= 0.0


def test_direction_agnostic():
    """As in the pairwise path, agreement does not depend on direction."""
    pos = [(0.90, 0.05), (0.85, 0.10), (0.88, 0.07)]
    neg = [(nu_, mu_) for mu_, nu_ in pos]
    a, _, _, _ = ix.run_indexed(pos, shots=SHOTS, seed=5)
    b, _, _, _ = ix.run_indexed(neg, shots=SHOTS, seed=5)
    assert abs(a - b) < TOL


def test_seed_reproduces():
    ops = [(0.85, 0.10), (0.80, 0.15), (0.70, 0.20)]
    a = ix.run_indexed(ops, shots=20_000, seed=99)
    b = ix.run_indexed(ops, shots=20_000, seed=99)
    assert a == b


def test_rejects_single_opinion():
    try:
        ix.build_indexed_circuit([(0.7, 0.2)])
    except ValueError:
        return
    raise AssertionError("accepted n = 1")


def test_width_is_linear_in_n():
    """qubits = 2n + 2*ceil(log2 n) + 12."""
    import math

    for n in (2, 3, 4, 8, 16):
        m = max(1, math.ceil(math.log2(n)))
        assert ix.metrics(n)["qubits"] == 2 * n + 2 * m + 12


def test_pairwise_circuit_depth_is_constant():
    """
    The counterpoint that justifies the design choice: the pairwise circuit
    has FIXED width and depth, whereas the indexed one grows with n.
    """
    from qif_xnor import build_xnor_circuit

    base = build_xnor_circuit((0.7, 0.2), (0.6, 0.3))
    assert base.num_qubits == 10
    # 9 active stages in the T0--T11 ordering of Fig. 2; fixed, independent of n
    assert base.depth() == 9
    assert base.count_ops()["ccx"] == 6

    deep = [ix.metrics(n)["depth"] for n in (2, 4, 8)]
    assert deep[0] < deep[1] < deep[2]


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
