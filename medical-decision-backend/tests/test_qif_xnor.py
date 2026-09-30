"""
Validates the QIF-XNOR Engine against the formulation of the manuscript.

Runs under pytest or directly:  python3 tests/test_qif_xnor.py
"""

import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qif_xnor import (  # noqa: E402
    IFV_ONE,
    IFV_ZERO,
    N_MEASURE_OPS,
    N_PREP_OPS,
    STAGE_OPS,
    build_xnor_circuit,
    figure_ops,
    run_xnor,
    stage_table,
    xnor_analytic,
    xor_analytic,
)

# Stages T1--T11 of Fig. 2, gate by gate. Does NOT include the amplitude
# encoding (T0 of the figure is the already-prepared state) nor the
# measurements, which the manuscript applies after the circuit.
FIGURE_SEQUENCE = [
    ("x", (1,)), ("x", (2,)),                                       # T1
    ("ccx", (1, 2, 4)),                                             # T2
    ("x", (0,)), ("x", (1,)), ("x", (2,)), ("x", (3,)), ("x", (4,)),  # T3
    ("ccx", (0, 3, 5)),                                             # T4
    ("x", (0,)), ("x", (3,)), ("x", (5,)),                          # T5
    ("ccx", (4, 5, 8)),                                             # T6
    ("ccx", (0, 3, 6)),                                             # T7
    ("ccx", (1, 2, 7)),                                             # T8
    ("x", (6,)), ("x", (7,)),                                       # T9
    ("ccx", (6, 7, 9)),                                             # T10
    ("x", (6,)), ("x", (7,)), ("x", (9,)),                          # T11
]

SHOTS = 200_000
MC_TOL = 6e-3


def _grid(step=1 / 12):
    n = int(round(1 / step))
    vals = [i * step for i in range(n + 1)]
    for mu, nu in itertools.product(vals, repeat=2):
        if mu + nu <= 1 + 1e-12:
            yield (mu, nu)


PAIRS = [(x, y) for x, y in itertools.product(list(_grid()), repeat=2)]


# =====================================================
# Correspondencia literal com a Fig. 2
# =====================================================

def test_figure_stages_match_gate_for_gate():
    """
    Locks the exact gate order of T1--T11. The circuit must reproduce the
    figure literally, not a mathematically equivalent realisation.
    """
    qc = build_xnor_circuit((0.9, 0.05), (0.85, 0.10))
    got = [
        (ci.operation.name, tuple(qc.find_bit(q).index for q in ci.qubits))
        for ci in figure_ops(qc)
    ]
    assert got == FIGURE_SEQUENCE, got


def test_figure_gate_counts():
    """15 Pauli-X and 6 Toffoli gates, counted inside T1--T11 only."""
    from collections import Counter

    qc = build_xnor_circuit((0.7, 0.2), (0.6, 0.3))
    counts = Counter(ci.operation.name for ci in figure_ops(qc))
    assert qc.num_qubits == 10
    assert counts["x"] == 15
    assert counts["ccx"] == 6
    assert set(counts) == {"x", "ccx"}      # nada mais pertence a figura


def test_state_preparation_is_separate_from_the_figure():
    """T0 of the figure is the already-prepared state; the u gates precede T1."""
    qc = build_xnor_circuit((0.7, 0.2), (0.6, 0.3))
    prep = list(qc.data[:N_PREP_OPS])
    assert len(prep) == 4
    assert all(ci.operation.name == "u" for ci in prep)
    assert [qc.find_bit(ci.qubits[0]).index for ci in prep] == [0, 1, 2, 3]
    assert not any(ci.operation.name == "u" for ci in figure_ops(qc))


def test_measurement_is_applied_after_the_figure():
    """The figure ends at T11; measurement comes after and is not a stage."""
    qc = build_xnor_circuit((0.9, 0.05), (0.85, 0.10))
    tail = list(qc.data[len(qc.data) - N_MEASURE_OPS:])
    assert [(qc.find_bit(m.qubits[0]).index, qc.find_bit(m.clbits[0]).index)
            for m in tail] == [(8, 0), (9, 1)]
    assert not any(ci.operation.name == "measure" for ci in figure_ops(qc))


def test_stage_table_covers_the_figure():
    """STAGE_OPS must not diverge from the constructed circuit."""
    qc = build_xnor_circuit((0.7, 0.2), (0.6, 0.3))
    assert sum(n for _, n in STAGE_OPS) == len(figure_ops(qc))
    rows = stage_table()
    assert [r["stage"] for r in rows] == ["T0"] + [f"T{i}" for i in range(1, 12)] + ["post"]
    assert [r["in_figure"] for r in rows] == [False] + [True] * 11 + [False]
    assert sum(len(r["ops"]) for r in rows) == len(qc.data)


def test_outputs_follow_the_figure_convention():
    """q8 -> c0 = mu_boxplus  e  q9 -> c1 = nu_boxplus."""
    mu_t, nu_t = xnor_analytic((0.9, 0.05), (0.85, 0.10))
    mu_c, nu_c, _ = run_xnor((0.9, 0.05), (0.85, 0.10), shots=SHOTS)
    assert abs(mu_c - mu_t) < MC_TOL and abs(nu_c - nu_t) < MC_TOL
    assert mu_t > nu_t          # case chosen so the reading is unambiguous


# =====================================================
# Axiomas da Definicao 1
# =====================================================

def test_DI1_boundary():
    """D_I1: XNOR(0~,0~) = XNOR(1~,1~) = 1~, e XNOR(0~,1~) = 0~."""
    assert xnor_analytic(IFV_ZERO, IFV_ZERO) == (1.0, 0.0)
    assert xnor_analytic(IFV_ONE, IFV_ONE) == (1.0, 0.0)
    assert xnor_analytic(IFV_ZERO, IFV_ONE) == (0.0, 1.0)
    assert xnor_analytic(IFV_ONE, IFV_ZERO) == (0.0, 1.0)


def test_DI2_symmetry():
    """D_I2: XNOR e simetrico."""
    for x, y in PAIRS:
        assert xnor_analytic(x, y) == xnor_analytic(y, x)


def test_DI3_identity():
    """D_I3: XNOR(1~, x~) = x~ -- o axioma violado pela versao anterior."""
    for x in _grid():
        mu, nu = xnor_analytic(IFV_ONE, x)
        assert abs(mu - x[0]) < 1e-12, (x, mu)
        assert abs(nu - x[1]) < 1e-12, (x, nu)


# =====================================================
# Condicao intuicionista e dualidade N_I
# =====================================================

def test_closed_in_U():
    """Eq.(1): the result stays in U~, hence pi >= 0 without clamping."""
    for x, y in PAIRS:
        mu, nu = xnor_analytic(x, y)
        assert 0.0 <= mu <= 1.0
        assert 0.0 <= nu <= 1.0
        assert mu + nu <= 1.0 + 1e-12, (x, y, mu + nu)


def test_NI_duality_is_a_swap():
    """XNOR = swap(XOR): the N_I-dual transposes components, it does not apply 1-x."""
    for x, y in PAIRS:
        mu_n, nu_n = xnor_analytic(x, y)
        mu_o, nu_o = xor_analytic(x, y)
        assert abs(mu_n - nu_o) < 1e-12
        assert abs(nu_n - mu_o) < 1e-12


def test_bitwise_negation_would_break_U():
    """
    Explicit regression test: the (1-mu, 1-nu) formulation leaves U~ over the
    vast majority of U~. This test fails if anyone reintroduces qc.x(8)/x(9).
    """
    violations = sum(
        1
        for x, y in PAIRS
        if (1 - xor_analytic(x, y)[0]) + (1 - xor_analytic(x, y)[1]) > 1 + 1e-9
    )
    assert violations > 0.9 * len(PAIRS)


# =====================================================
# Quantum circuit vs Proposition 2
# =====================================================

def test_circuit_matches_proposition_2():
    """As marginais de q9/q8 reproduzem mu_boxplus e nu_boxplus."""
    cases = [
        ((0.95, 0.02), (0.90, 0.05)),   # strong agreement, positive
        ((0.05, 0.90), (0.02, 0.95)),   # strong agreement, negative
        ((0.90, 0.05), (0.05, 0.90)),   # total disagreement
        ((0.70, 0.20), (0.65, 0.25)),   # moderate agreement
        ((0.30, 0.30), (0.30, 0.30)),   # high hesitation
        ((1.00, 0.00), (0.70, 0.20)),   # D_I3 on the circuit
        ((0.00, 1.00), (0.00, 1.00)),   # D_I1
    ]
    for x, y in cases:
        mu_t, nu_t = xnor_analytic(x, y)
        mu_c, nu_c, pi_c = run_xnor(x, y, shots=SHOTS)
        assert abs(mu_c - mu_t) < MC_TOL, (x, y, mu_c, mu_t)
        assert abs(nu_c - nu_t) < MC_TOL, (x, y, nu_c, nu_t)
        assert mu_c + nu_c <= 1.0 + MC_TOL
        assert abs(pi_c - max(0.0, 1 - mu_t - nu_t)) < 2 * MC_TOL


def test_circuit_is_deterministic_at_extremes():
    """At the vertices of U~ the circuit has no variance."""
    for x, y, exp in [
        (IFV_ONE, IFV_ONE, (1.0, 0.0)),
        (IFV_ZERO, IFV_ZERO, (1.0, 0.0)),
        (IFV_ONE, IFV_ZERO, (0.0, 1.0)),
    ]:
        mu, nu, _ = run_xnor(x, y, shots=2000)
        assert (mu, nu) == exp, (x, y, mu, nu)


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
    print(f"\n{len(fns) - failed}/{len(fns)} tests passed ({len(PAIRS)} pairs in U~)")
    sys.exit(1 if failed else 0)
