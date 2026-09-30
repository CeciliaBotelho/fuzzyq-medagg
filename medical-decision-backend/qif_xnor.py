"""
QIF-XNOR Engine
===============

Implementation of the intuitionistic fuzzy XNOR operator (\\boxplus_I) and of
its realisation as a quantum circuit (Quantum Agreement Register).

Reference: the accompanying manuscript, "From Intuitionistic Fuzzy Evidence
to Clinical Actions: A Quantum-Fuzzy Approach for Agreement-Based Aggregation"
(under review; authors omitted for double-blind review).

Proposition 2 (Eqs. eq-ts / eq-st), with T = T_P (product) and S = S_P
(probabilistic sum):

    mu_boxplus(x,y)  = T( S(x2,y1), S(x1,y2) )
    nu_boxplus(x,y)  = S( T(x1,y2), T(x2,y1) )

REGISTER CONVENTION
-------------------
The qubit assignment follows Fig. 2 of the manuscript, stages T0--T11:

    q0..q3   encoded inputs x1, x2, y1, y2
    q4       S_P(x2,y1)        t-conorm ancilla
    q5       S_P(x1,y2)        t-conorm ancilla
    q6       T_P(x1,y2)        t-norm ancilla
    q7       T_P(x2,y1)        t-norm ancilla
    q8       T_P(q4,q5) = mu_boxplus     -> c0     [Eq. eq-ts]
    q9       S_P(q6,q7) = nu_boxplus     -> c1     [Eq. eq-st]

The gate sequence of stages T1--T11 reproduces the figure literally, not a
mathematically equivalent realisation: 15 Pauli-X and 6 Toffoli gates.

Two stages fall OUTSIDE that literal comparison:

  T0            in the figure this is the input state ALREADY prepared. In the
                code, amplitude encoding (the u gates) is a state-preparation
                step that precedes T1.
  measurement   the figure ends at the three Pauli-X gates of T11. The
                measurements of q8 and q9 are applied AFTER the circuit, as
                described in the manuscript, and are not a stage.

STAGE_OPS fixes the operation count per stage from T1 to T11; stage_table()
builds the table from the circuit that is actually constructed.

ON NEGATION
-----------
The XNOR is the N_I-dual of the XOR under the standard intuitionistic
negation N_I_S(x1,x2) = (x2,x1), which TRANSPOSES the components -- it is not
the fuzzy negation N_S(x) = 1-x. Applying Pauli-X to the outputs, that is
reading (1 - mu_oplus, 1 - nu_oplus), violates the intuitionistic condition
mu + nu <= 1 and breaks axiom D_I3.

The Pauli-X gates appearing in the circuit below are internal De Morgan steps,
S(a,b) = not T(not a, not b), and not the intuitionistic negation.
"""

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import qiskit
import qiskit_aer
from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer

__all__ = [
    "IFV",
    "T_P",
    "S_P",
    "xnor_analytic",
    "xor_analytic",
    "encode_inputs",
    "build_xnor_circuit",
    "run_xnor",
    "run_xnor_joint",
    "figure_ops",
    "stage_table",
    "derive_seeds",
    "provenance",
]

BACKEND_NAME = "qasm_simulator"

# Intuitionistic fuzzy value (mu, nu) with mu + nu <= 1.
IFV = Tuple[float, float]

# Extreme elements of U~ (Section III-A).
IFV_ONE: IFV = (1.0, 0.0)   # 1~ -- neutral element of boxplus_I (axiom D_I3)
IFV_ZERO: IFV = (0.0, 1.0)  # 0~

_TOL = 1e-9


# =====================================================
# 1) FUZZY AGGREGATIONS (Eqs. tp / sp)
# =====================================================

def T_P(a: float, b: float) -> float:
    """Product t-norm."""
    return a * b


def S_P(a: float, b: float) -> float:
    """Probabilistic-sum t-conorm (dual of T_P)."""
    return a + b - a * b


# =====================================================
# 2) ANALYTICAL REFERENCE (Proposition 2)
# =====================================================

def xnor_analytic(x: IFV, y: IFV) -> IFV:
    """Closed form of boxplus_I(x,y), used as the oracle in the tests."""
    x1, x2 = x
    y1, y2 = y
    mu = T_P(S_P(x2, y1), S_P(x1, y2))
    nu = S_P(T_P(x1, y2), T_P(x2, y1))
    return mu, nu


def xor_analytic(x: IFV, y: IFV) -> IFV:
    """oplus_I(x,y). Note that xnor_analytic(x,y) == swap(xor_analytic(x,y))."""
    mu, nu = xnor_analytic(x, y)
    return nu, mu


# =====================================================
# 3) QUBIT ENCODING (Section III-B)
# =====================================================

def encode_inputs(qc: QuantumCircuit, values: Sequence[float]) -> None:
    """
    Encodes each fuzzy value x as
        |psi> = sqrt(1-x)|0> + sqrt(x)|1>
    through U(theta, 0, pi), with theta = 2*arctan( sqrt(x) / sqrt(1-x) ),
    so that P(|1>) = sin^2(theta/2) = x.
    """
    for qubit_index, x in enumerate(values):
        if x <= 0.0:
            theta = 0.0
        elif x >= 1.0:
            theta = np.pi
        else:
            theta = 2.0 * np.arctan(np.sqrt(x) / np.sqrt(1.0 - x))
        qc.u(theta, 0, np.pi, qubit_index)


# =====================================================
# 4) QUANTUM AGREEMENT REGISTER
# =====================================================

# Operations per stage, in the order emitted by build_xnor_circuit.
# Covers T1--T11 only: state preparation and measurement are excluded.
N_PREP_OPS = 4          # amplitude-encoding u gates (before T1)
N_MEASURE_OPS = 2       # measurements of q8 and q9 (after T11)

STAGE_OPS = [
    ("T1", 2),    # X(q1), X(q2)
    ("T2", 1),    # CCX(q1,q2 -> q4)
    ("T3", 5),    # X(q0), X(q1), X(q2), X(q3), X(q4)
    ("T4", 1),    # CCX(q0,q3 -> q5)
    ("T5", 3),    # X(q0), X(q3), X(q5)
    ("T6", 1),    # CCX(q4,q5 -> q8)
    ("T7", 1),    # CCX(q0,q3 -> q6)
    ("T8", 1),    # CCX(q1,q2 -> q7)
    ("T9", 2),    # X(q6), X(q7)
    ("T10", 1),   # CCX(q6,q7 -> q9)
    ("T11", 3),   # X(q6), X(q7), X(q9)
]

STAGE_ROLE = {
    "T0": "prepared input state (amplitude encoding, state preparation)",
    "T1": "negate x2, y1",
    "T2": "q4 <- T_P(not x2, not y1)",
    "T3": "negate x1, y2; restore x2, y1; q4 = S_P(x2,y1)",
    "T4": "q5 <- T_P(not x1, not y2)",
    "T5": "restore x1, y2; q5 = S_P(x1,y2)",
    "T6": "q8 = T_P(q4,q5) = mu_boxplus",
    "T7": "q6 <- T_P(x1,y2)",
    "T8": "q7 <- T_P(x2,y1)",
    "T9": "negate q6, q7",
    "T10": "q9 <- T_P(not q6, not q7)",
    "T11": "restore q6, q7; q9 = S_P(q6,q7) = nu_boxplus",
}


def build_xnor_circuit(x: IFV, y: IFV) -> QuantumCircuit:
    """
    Intuitionistic fuzzy XNOR circuit.

    Stages T1--T11 follow Fig. 2 gate for gate; see STAGE_OPS. Amplitude
    encoding precedes T1 as state preparation and the measurements are applied
    after T11; neither belongs to the figure.

    Because every gate after the encoding is a permutation of the
    computational basis (X and Toffoli), the marginal probabilities of q8 and
    q9 reproduce Proposition 2 exactly under T_P/S_P.
    """
    x1, x2 = x
    y1, y2 = y

    qc = QuantumCircuit(10, 2, name="xnor_intuitionistic")

    # ---- state preparation (T0 of the figure: inputs already encoded) ----
    encode_inputs(qc, [x1, x2, y1, y2])

    # ---- T1 ----
    qc.x(1)
    qc.x(2)
    # ---- T2 ----
    qc.ccx(1, 2, 4)
    # ---- T3 ----
    qc.x(0)
    qc.x(1)
    qc.x(2)
    qc.x(3)
    qc.x(4)
    # ---- T4 ----
    qc.ccx(0, 3, 5)
    # ---- T5 ----
    qc.x(0)
    qc.x(3)
    qc.x(5)
    # ---- T6 ----
    qc.ccx(4, 5, 8)
    # ---- T7 ----
    qc.ccx(0, 3, 6)
    # ---- T8 ----
    qc.ccx(1, 2, 7)
    # ---- T9 ----
    qc.x(6)
    qc.x(7)
    # ---- T10 ----
    qc.ccx(6, 7, 9)
    # ---- T11 ----
    qc.x(6)
    qc.x(7)
    qc.x(9)

    # ---- measurement, applied after the figure's circuit ----
    qc.measure(8, 0)         # c0 = mu_boxplus
    qc.measure(9, 1)         # c1 = nu_boxplus
    return qc


def figure_ops(qc: QuantumCircuit) -> List:
    """The instructions belonging to stages T1--T11 of the figure."""
    return list(qc.data[N_PREP_OPS:len(qc.data) - N_MEASURE_OPS])


def stage_table(x: IFV = (0.90, 0.05), y: IFV = (0.85, 0.10)) -> List[Dict]:
    """
    Builds the stage table from the circuit that is actually constructed.

    There is no hand-written gate list: the operations are read from qc.data
    and merely sliced according to STAGE_OPS, whose total is checked against
    the number of operations in T1--T11. The T0 row describes the prepared
    state and the final row describes the measurement; neither belongs to the
    figure.
    """
    qc = build_xnor_circuit(x, y)
    fig = figure_ops(qc)
    total = sum(n for _, n in STAGE_OPS)
    if total != len(fig):
        raise AssertionError(
            f"STAGE_OPS soma {total} operacoes, T1--T11 tem {len(fig)}"
        )

    def describe(instrs):
        out = []
        for ci in instrs:
            out.append({
                "gate": ci.operation.name,
                "qubits": [qc.find_bit(q).index for q in ci.qubits],
                "clbits": [qc.find_bit(c).index for c in ci.clbits],
            })
        return out

    rows = [{
        "stage": "T0",
        "role": STAGE_ROLE["T0"],
        "in_figure": False,
        "ops": describe(qc.data[:N_PREP_OPS]),
    }]

    k = 0
    for label, count in STAGE_OPS:
        rows.append({
            "stage": label,
            "role": STAGE_ROLE[label],
            "in_figure": True,
            "ops": describe(fig[k:k + count]),
        })
        k += count

    rows.append({
        "stage": "post",
        "role": "measurement, applied after the figure's circuit",
        "in_figure": False,
        "ops": describe(qc.data[len(qc.data) - N_MEASURE_OPS:]),
    })
    return rows


def derive_seeds(base: int, count: int) -> List[int]:
    """
    Derives `count` independent, reproducible seeds from a base seed.

    One aggregation runs C(n,2) circuits. Reusing the same seed in all of them
    correlates the sampling errors, and in the degenerate case of two pairs
    with identical inputs the simulator would return identical counts,
    cancelling the benefit of averaging.

    The derivation uses numpy.random.SeedSequence.spawn, the standard
    mechanism for generating independent streams from a single entropy source,
    rather than an ad-hoc multiplicative spread. The values are reduced to the
    range accepted by Aer's seed_simulator.
    """
    children = np.random.SeedSequence(base).spawn(count)
    return [
        int(c.generate_state(1, dtype=np.uint32)[0]) % (2**31 - 1)
        for c in children
    ]


def provenance(shots: int, seed: Optional[int] = None) -> Dict[str, object]:
    """Metadata that makes a published result traceable."""
    return {
        "qiskit": qiskit.__version__,
        "qiskit_aer": qiskit_aer.__version__,
        "backend": BACKEND_NAME,
        "shots": shots,
        "seed": seed,
    }


def run_xnor(
    x: IFV,
    y: IFV,
    shots: int = 20000,
    seed: Optional[int] = None,
) -> Tuple[float, float, float]:
    """
    Runs the circuit and returns (mu_boxplus, nu_boxplus, pi_boxplus).

    pi_boxplus = 1 - mu - nu is the residual hesitation of the agreement.
    Proposition 2 guarantees mu + nu <= 1, so pi >= 0 up to sampling error and
    no clamping is required.

    seed=None keeps the run stochastic, which is the honest default: silently
    fixing a seed would make a sampling method look deterministic. Pass an
    explicit seed to reproduce published numbers.
    """
    j = run_xnor_joint(x, y, shots=shots, seed=seed)
    return j["mu"], j["nu"], max(0.0, 1.0 - j["mu"] - j["nu"])


def run_xnor_joint(
    x: IFV,
    y: IFV,
    shots: int = 20000,
    seed: Optional[int] = None,
) -> Dict[str, float]:
    """
    Like run_xnor, but also returns p11 = P(mu=1 and nu=1).

    mu and nu are read from the SAME measurement and are therefore
    correlated; p11 is what allows Cov(mu, nu) to be estimated, and with it
    the confidence interval of pi = 1 - mu - nu.
    """
    qc = build_xnor_circuit(x, y)
    backend = Aer.get_backend(BACKEND_NAME)
    job = backend.run(transpile(qc, backend), shots=shots, seed_simulator=seed)
    counts = job.result().get_counts()

    total = sum(counts.values())
    # Qiskit returns the bit string in little-endian order, "c1c0",
    # with c0 = q8 = mu_boxplus and c1 = q9 = nu_boxplus.
    mu = sum(c for b, c in counts.items() if b[1] == "1") / total
    nu = sum(c for b, c in counts.items() if b[0] == "1") / total
    p11 = sum(c for b, c in counts.items() if b == "11") / total

    return {"mu": mu, "nu": nu, "p11": p11, "shots": total}
