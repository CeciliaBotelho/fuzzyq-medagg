"""
Indexed Quantum Agreement Register
==================================

ADDITIONAL ANALYSIS -- design-space comparison, not one of the main
experiments of the manuscript. See experiments/README.md.

Evaluates the MEAN agreement among n opinions in ONE circuit, instead of
C(n,2) independent executions.

Idea
----
Two index registers a and b, each with m = ceil(log2 n) qubits, are placed in
uniform superposition. A multiplexed read brings the opinion selected by each
index into a 4-qubit working register, on which the same XNOR circuit of
qif_xnor is applied. Post-selecting the branches with a != b and both indices
in range,

    P(mu = 1 | valid) = 1/(n(n-1)) * SUM_{a != b} mu_boxplus(x_a, x_b)

and, since boxplus_I is symmetric (D_I2), the mean over ORDERED pairs
coincides with the mean over unordered pairs. This is exactly the uniform
aggregation of fusion.aggregate.

The shot budget then appears amortised across all pairs: a single execution
samples the mean, instead of C(n,2) executions each sampling one term.

Honest limitation
-----------------
The register grows as 2n + 2m + 12 qubits. Under classical SIMULATION the
cost is exponential in that width, so the apparent shot saving does not
translate into wall-clock time. Measurements reported in
experiments/compare_registers.py show no accuracy gain at an equal total shot
budget, while width and depth both increase. See metrics().
"""

import itertools
import math
from typing import Dict, List, Optional, Sequence, Tuple

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile
from qiskit_aer import Aer

from qif_xnor import BACKEND_NAME, IFV, encode_inputs, xnor_analytic

__all__ = ["build_indexed_circuit", "run_indexed", "metrics", "mean_pairwise_analytic"]


def _index_width(n: int) -> int:
    return max(1, math.ceil(math.log2(n)))


def mean_pairwise_analytic(opinions: Sequence[IFV]) -> IFV:
    """Closed-form reference: mean of boxplus_I over all pairs i < j."""
    pairs = list(itertools.combinations(range(len(opinions)), 2))
    mus, nus = [], []
    for i, j in pairs:
        mu, nu = xnor_analytic(opinions[i], opinions[j])
        mus.append(mu)
        nus.append(nu)
    return sum(mus) / len(pairs), sum(nus) / len(pairs)


# =====================================================
# Multiplexed read
# =====================================================

def _select_into(
    qc: QuantumCircuit,
    idx: QuantumRegister,
    ops: QuantumRegister,
    dest: List,
    sel,
    n: int,
) -> None:
    """
    For each value a in [0, n), computes sel = (idx == a) and, conditioned on
    sel, swaps the qubit pair of opinion a into dest.

    Swapping rather than copying respects the no-cloning theorem: in the
    branch where idx = a, opinion a migrates into the working register and its
    slot is left in |0>. Branches with a == b are discarded by post-selection,
    so the second read never needs an opinion already consumed in a valid
    branch.
    """
    m = len(idx)
    for a in range(n):
        neg = [idx[k] for k in range(m) if not (a >> k) & 1]

        for q in neg:
            qc.x(q)
        qc.mcx(list(idx), sel)
        for q in neg:
            qc.x(q)

        qc.cswap(sel, ops[2 * a], dest[0])
        qc.cswap(sel, ops[2 * a + 1], dest[1])

        for q in neg:                      # uncompute sel
            qc.x(q)
        qc.mcx(list(idx), sel)
        for q in neg:
            qc.x(q)


def _mark_invalid(
    qc: QuantumCircuit,
    idx_a: QuantumRegister,
    idx_b: QuantumRegister,
    bad,
    n: int,
) -> None:
    """
    Marks bad = 1 when either index is out of range or when a == b.

    Equality is detected without extra ancillas: idx_a is added bit by bit
    into idx_b (CNOTs), so idx_b comes to hold the difference; all-zero bits
    mean equal indices. The CNOTs are then undone.
    """
    m = len(idx_a)

    for reg in (idx_a, idx_b):             # out of range
        for v in range(n, 2 ** m):
            neg = [reg[k] for k in range(m) if not (v >> k) & 1]
            for q in neg:
                qc.x(q)
            qc.mcx(list(reg), bad)
            for q in neg:
                qc.x(q)

    for k in range(m):                     # idx_b <- idx_a XOR idx_b
        qc.cx(idx_a[k], idx_b[k])
    for k in range(m):
        qc.x(idx_b[k])
    qc.mcx(list(idx_b), bad)               # diferenca nula => iguais
    for k in range(m):
        qc.x(idx_b[k])
    for k in range(m):                     # desfaz
        qc.cx(idx_a[k], idx_b[k])


# =====================================================
# Circuit
# =====================================================

def build_indexed_circuit(opinions: Sequence[IFV]) -> QuantumCircuit:
    n = len(opinions)
    if n < 2:
        raise ValueError("the indexed register requires at least two opinions")

    m = _index_width(n)

    ops = QuantumRegister(2 * n, "op")
    idx_a = QuantumRegister(m, "ia")
    idx_b = QuantumRegister(m, "ib")
    work = QuantumRegister(4, "w")
    sel = QuantumRegister(1, "sel")
    anc = QuantumRegister(4, "anc")
    out = QuantumRegister(2, "out")
    bad = QuantumRegister(1, "bad")
    cr = ClassicalRegister(3, "c")

    qc = QuantumCircuit(ops, idx_a, idx_b, work, sel, anc, out, bad, cr,
                        name=f"indexed_xnor_n{n}")

    flat: List[float] = []
    for mu, nu in opinions:
        flat += [mu, nu]
    encode_inputs(qc, flat)                # opinions in amplitude

    qc.h(idx_a)                            # uniform superposition over pairs
    qc.h(idx_b)

    _select_into(qc, idx_a, ops, [work[0], work[1]], sel[0], n)
    _select_into(qc, idx_b, ops, [work[2], work[3]], sel[0], n)

    _mark_invalid(qc, idx_a, idx_b, bad[0], n)

    # ---- XNOR over the working register (same construction as qif_xnor) ----
    x1, x2, y1, y2 = work[0], work[1], work[2], work[3]

    qc.ccx(x2, y1, anc[0])                 # nu = S(T(x2,y1), T(x1,y2)) -> out[1]
    qc.ccx(x1, y2, anc[1])
    qc.x(anc[0])
    qc.x(anc[1])
    qc.ccx(anc[0], anc[1], out[1])
    qc.x(out[1])

    qc.x(x1)                               # mu = T(S(x1,y2), S(x2,y1)) -> out[0]
    qc.x(x2)
    qc.x(y1)
    qc.x(y2)
    qc.ccx(x1, y2, anc[2])
    qc.x(anc[2])
    qc.ccx(x2, y1, anc[3])
    qc.x(anc[3])
    qc.ccx(anc[2], anc[3], out[0])

    qc.measure(out[0], cr[0])              # mu
    qc.measure(out[1], cr[1])              # nu
    qc.measure(bad[0], cr[2])              # 0 = ramo valido
    return qc


def run_indexed(
    opinions: Sequence[IFV],
    shots: int = 20000,
    seed: Optional[int] = None,
    method: str = "automatic",
) -> Tuple[float, float, float, float]:
    """
    Returns (mean_mu, mean_nu, pi, acceptance_rate).

    The acceptance rate is the fraction of shots landing in valid branches;
    the effective shots for the estimate are shots * acceptance.
    """
    qc = build_indexed_circuit(opinions)
    backend = Aer.get_backend(BACKEND_NAME)
    if method != "automatic":
        backend.set_options(method=method)

    counts = (
        backend.run(transpile(qc, backend), shots=shots, seed_simulator=seed)
        .result()
        .get_counts()
    )

    total = sum(counts.values())
    valid = mu_hits = nu_hits = 0
    for bits, c in counts.items():         # bits little-endian: "c2c1c0"
        if bits[0] == "1":                 # bad
            continue
        valid += c
        if bits[2] == "1":
            mu_hits += c
        if bits[1] == "1":
            nu_hits += c

    if valid == 0:
        raise RuntimeError("nenhum ramo valido amostrado; aumente shots")

    mu = mu_hits / valid
    nu = nu_hits / valid
    return mu, nu, max(0.0, 1.0 - mu - nu), valid / total


def metrics(n: int, topology: str = "logical") -> Dict[str, object]:
    """
    Width, depth and resource estimate of the indexed circuit.

    topology="logical" measures the {cx, u} decomposition under free
    connectivity, which is a LOWER BOUND and not a hardware cost. Pass
    "heavy-hex" or "linear" to include routing. See resources.py.
    """
    import resources

    qc = build_indexed_circuit([(0.6, 0.3)] * n)
    est = resources.estimate(qc, topology)
    m = _index_width(n)
    return {
        "n": n,
        "pairs": n * (n - 1) // 2,
        "index_width": m,
        "qubits": qc.num_qubits,
        "depth": qc.depth(),
        "topology": topology,
        "hw_depth": est["depth"],
        "cx": est["cx"],
        "ops": dict(qc.count_ops()),
        "acceptance": n * (n - 1) / (2 ** (2 * m)),
    }
