# FuzzyQ-MedAgg

Agreement-based aggregation of intuitionistic fuzzy opinions through a quantum
XNOR circuit.

This repository accompanies a manuscript under review. Author, affiliation and
institutional information have been removed from the artifact.

## Contents

1. [Overview](#overview)
2. [Method](#method)
3. [Repository layout](#repository-layout)
4. [Installation](#installation)
5. [Running](#running)
6. [Reproducing the experiments](#reproducing-the-experiments)
7. [Summary of results](#summary-of-results)
8. [Scope and limitations](#scope-and-limitations)

## Overview

Clinicians increasingly consult several AI models for the same case. Those
models disagree, and majority voting or averaging compresses the disagreement
away, so the output cannot distinguish "both sources are unsure" from "the
sources contradict each other".

FuzzyQ-MedAgg keeps the two apart. Each source reports an Atanassov
intuitionistic fuzzy value x = (mu, nu): belief for, belief
against, and the hesitation pi = 1 - mu - nu left over. Agreement between
two sources is measured by an intuitionistic fuzzy XNOR evaluated on a quantum
circuit, and the aggregate of n sources is mapped onto one of three clinical
actions: treat, do not treat, or request further examinations.

Diagnostic responsibility remains with the human specialist. The framework
provides a structured and auditable mechanism for handling uncertainty and
conflict among multiple AI predictions.

## Method

### Agreement operator

With the product t-norm T_P(a,b) = ab and the probabilistic sum
S_P(a,b) = a + b - ab:

```
mu_XNOR = T_P(S_P(x2, y1), S_P(x1, y2))
nu_XNOR = S_P(T_P(x1, y2), T_P(x2, y1))
```

The XNOR is the dual of the XOR under the intuitionistic negation
N_IS(x1, x2) = (x2, x1), which transposes components. It is not the
fuzzy negation 1 - x: that pair would leave U~, violating
mu + nu <= 1, and would break the identity axiom D_I3.

### Quantum circuit

The manuscript figure labels twelve columns, T0 to T11. T0 is the input
register already amplitude-encoded, so it carries no gate; T1 to T11 hold
the 15 Pauli-X and 6 Toffoli gates, and the code reproduces those eleven
columns gate for gate. Measurement of q8 and q9 is applied after T11 and
is not drawn in the figure.

| Property | Value |
| --- | --- |
| Qubits | 10 (4 input, 4 ancilla, 2 output) |
| Logical gates | 15 Pauli-X, 6 Toffoli |
| Depth | 7 for the logical block, 9 for the full circuit |
| Outputs | q8 = mu_XNOR, q9 = nu_XNOR |

Width and depth are fixed and do not depend on the number of sources.

### Aggregating n sources

XNOR is a binary equivalence and is non-associative, so no fold is
available. Agreement is formed from the P = C(n,2) pairs:

```
mu_C = sum over i<j of  w_ij * mu_ij
nu_C = sum over i<j of  w_ij * nu_ij      with  sum over i<j of w_ij = 1
```

A convex combination of points of the convex set U~ remains in
U~, so pi_C >= 0 requires no truncation. The weights depend only
on the unordered pair, which makes the result independent of presentation
order. A single pipeline covers every n >= 2; for n = 2 it reduces
exactly to XNOR(x1, x2) with w_12 = 1.

Three weighting schemes are available. The uniform scheme sets
w_ij = 1/P. The source-specific scheme uses per-source weights. The institution-balanced scheme sets omega_i = 1 / (K * n_k), balancing the source weights across the
K institutions before the pair weights are induced; this attenuates the
dominance of a larger site without making the final contributions equal.

### From degrees to a clinical action

The admissibility result gives mu_C + nu_C + pi_C = 1 exactly, so the three degrees
already partition the unit. The score s = mu_bar - nu_bar, the same
function used in the Xu-Yager ordering, supplies the direction of the mean
evidence:

```
treat        = mu_C * (1 + s) / 2
doNotTreat   = mu_C * (1 - s) / 2
requestExams = nu_C + pi_C
```

These sum to 1 by construction. The mapping is the decision layer of the
framework and is not a consequence of the Xu-Yager ordering.

## Repository layout

```
medical-decision-backend/
  qif_xnor.py            QIF-XNOR Engine: encoding, circuit T1-T11, closed form
  fusion.py              Evidence Fusion Core: pairwise aggregation, weighting,
                         confidence intervals, decision layer
  resources.py           circuit resource estimates under explicit topologies
  main.py                REST API: POST /decide, GET /health
  requirements.txt
  tests/                 68 tests, runnable without pytest
  experiments/
    tfs_evaluation/      every number reported in the manuscript
    compare_registers.py         additional analysis, design space
    compare_decision_rules.py    additional analysis, ablation

experiments/             manuscript tables and the gate-noise study
medical-decision-front/  demonstration prototype, produces no reported result
data/example_dataset.csv sample batch of cases
```

## Installation

Python 3.11 or later is required. Node.js 20 and pnpm are needed only for the
demonstration interface, on which no reported result depends.

```bash
cd medical-decision-backend
python3 -m venv .venv && source .venv/bin/activate   # optional
pip install -r requirements.txt
```

## Running

Test suite, 68 tests, no pytest required:

```bash
cd medical-decision-backend
for t in tests/*.py; do python3 "$t"; done
```

REST API on port 8000:

```bash
cd medical-decision-backend
uvicorn main:app --reload
curl -s localhost:8000/health
```

```bash
curl -X POST http://127.0.0.1:8000/decide \
  -H 'Content-Type: application/json' \
  -d '{"opinions":[{"mu":0.85,"nu":0.10},{"mu":0.80,"nu":0.02}], "exact":true}'
```

```json
{
  "decision": "TREAT",
  "scores":    { "treat": 0.617, "doNotTreat": 0.082, "requestExams": 0.301 },
  "consensus": { "mu": 0.699, "nu": 0.096, "pi": 0.205 },
  "provenance": { "backend": "closed form (Proposition 2)" }
}
```

Setting `"exact": false`, the default, samples the circuit instead. The
response then carries the 95% intervals together with the Qiskit and Aer
versions, the shot count and the seed.

Demonstration interface on port 3000, with the API already running:

```bash
cd medical-decision-front
cp .env.example .env.local     # only if the API is not on port 8000
pnpm install && pnpm dev
```

## Reproducing the experiments

All experiments are analytical or simulator-based. No quantum hardware and no
GPU are required. The full set runs in under ten minutes on a laptop.

### Requirements

Python 3.11, `qiskit==2.5.2`, `qiskit-aer==0.17.2`, `numpy`. Install into a
clean environment:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install "qiskit==2.5.2" "qiskit-aer==0.17.2" numpy
```

The exact versions matter only for bit-for-bit reproduction of the reported
figures. Any Qiskit 2.x with Aer 0.17 reproduces the same values to within
sampling error.

### Layout

```
experiments/
  qif_xnor_noise.py        circuit statistics and gate-noise study
  aggregation_tables.py    baselines and multi-source stress scenarios
  noise_results.csv        output of qif_xnor_noise.py, committed for reference
```

### 1. Circuit statistics and gate-noise sensitivity

```bash
python experiments/qif_xnor_noise.py
```

Runtime: about four minutes. The script writes `noise_results.csv` to the
working directory, so run it from the repository root to overwrite the
committed copy, or from `experiments/` to keep it in place.

The script rebuilds the ten-qubit QIF-XNOR circuit from the gate sequence of
the manuscript and measures how the operator degrades under depolarizing
noise. The circuit it builds carries 15 Pauli-X gates, 6 Toffoli gates, 4
state-preparation rotations and 2 measurements, with logical depth 9, and
transpiles to 36 CX gates at depth 31 in the {CX, U} basis:

```
logical_ops : {'x': 15, 'ccx': 6, 'ry': 4, 'measure': 2}
logical_depth: 9
cx_opt1      : 36      depth_opt1: 31
```

It then sweeps five noise levels. At each level the error against the
analytical operator is computed twice and printed side by side:

- exact, from the marginals of a density-matrix simulation, which involves no
  sampling and isolates the effect of gate noise;
- shot-based, from 30 independent runs of 20,000 shots with `seed_simulator`
  in 1..30, reported as mean plus or minus standard deviation over the 30 runs.

Expected output:

```
p1=0e+00: exact MAE mu=0.0000 nu=0.0000 | sampled (30x20000) mu=0.0018+-0.0006 nu=0.0018+-0.0008
p1=1e-04: exact MAE mu=0.0028 nu=0.0060 | sampled (30x20000) mu=0.0033+-0.0010 nu=0.0062+-0.0013
p1=5e-04: exact MAE mu=0.0136 nu=0.0290 | sampled (30x20000) mu=0.0139+-0.0012 nu=0.0290+-0.0017
p1=1e-03: exact MAE mu=0.0266 nu=0.0558 | sampled (30x20000) mu=0.0268+-0.0013 nu=0.0556+-0.0018
p1=5e-03: exact MAE mu=0.1079 nu=0.2053 | sampled (30x20000) mu=0.1082+-0.0011 nu=0.2051+-0.0020
```

The two estimates agree within one standard deviation at every level, which
attributes the degradation to gate noise rather than to finite sampling.

The noise model applies error probability p1 to every single-qubit gate and
p2 = 10 * p1 to every CX gate. The tenfold ratio is an experimental scenario
reflecting the typically dominant two-qubit error; it is not a hardware law.

### 2. Baselines and multi-source stress scenarios

```bash
python experiments/aggregation_tables.py
```

Runtime: instantaneous. This script is purely analytical and requires no
quantum simulation. It regenerates two blocks of results.

The first contrasts the XNOR-based consensus with two intuitionistic fuzzy
averaging operators, IFWA and IFWG, and with the mean pairwise
Szmidt-Kacprzyk similarity 1 - d_H, across five synthetic scenarios. The
similarity is included as a reference comparator, not as an aggregator: it
returns a reflexive scalar, whereas the proposed operator returns a triple
(mu_C, nu_C, pi_C).

```
Concordant support       IFWA=(0.827,0.045) IFWG=(0.825,0.061) 1-dH=0.870 C=(0.699,0.096,0.205) Treat
Concordant rejection     IFWA=(0.061,0.825) IFWG=(0.045,0.827) 1-dH=0.870 C=(0.699,0.096,0.205) DoNotTreat
Symmetric conflict       IFWA=(0.633,0.292) IFWG=(0.292,0.633) 1-dH=0.250 C=(0.186,0.725,0.089) Exams
High hesitation          IFWA=(0.275,0.224) IFWG=(0.274,0.225) 1-dH=0.950 C=(0.190,0.121,0.689) Exams
One dissonant of three   IFWA=(0.700,0.162) IFWG=(0.408,0.496) 1-dH=0.450 C=(0.343,0.509,0.148) Exams
```

Neither averaging operator represents pairwise concordance or conflict, and
both take a side the evidence does not support: for the symmetric conflict
IFWA returns a support-oriented aggregate and IFWG its mirror image. The
similarity, being reflexive, rates the two highly hesitant opinions as almost
identical, whereas the proposed consensus reports that their agreement is
itself hesitant.

The second block stresses the multi-source aggregation with three concordant
sources from one institution and a fourth from another.

```
A only                           n=3 C=(0.701,0.122,0.177) s=+0.747 Treat
A + D1, uniform                  n=4 C=(0.433,0.412,0.155) s=+0.372 Exams
A + D1 duplicated                n=5 C=(0.384,0.474,0.142) s=+0.148 Exams
A + A1 duplicated + D1           n=5 C=(0.496,0.362,0.142) s=+0.448 Exams
A + D1, institution-balanced     n=4 C=(0.299,0.557,0.144) s=+0.372 Exams
A + D2, uniform                  n=4 C=(0.573,0.257,0.170) s=+0.560 Treat
A + D2, institution-balanced     n=4 C=(0.509,0.325,0.167) s=+0.560 Exams
Near-copies of A1 + D1           n=4 C=(0.467,0.444,0.089) s=+0.375 Exams
```

Two points are worth noting. Duplicating a source shifts the consensus, which
exposes the sensitivity of uniform weights to source multiplicity. And for the
hesitant minority D2 the weighting scheme changes the recommendation, from
treat under uniform weights to request further examinations under
institution-balanced weights.

### 3. Full evaluation suite

The remaining reported numbers come from
`medical-decision-backend/experiments/tfs_evaluation/`. Those scripts import
`qif_xnor`, `fusion` and `resources` and only observe them, changing no
production behaviour. The generated CSV and JSON files are committed alongside
the scripts.

```bash
cd medical-decision-backend/experiments/tfs_evaluation
python3 run_all.py            # E1 through E6, about seven minutes
```

| Manuscript experiment | Script | Output |
| --- | --- | --- |
| Setup and circuit resources | `e1_setup.py` | `experimental_setup.json` |
| XNOR validation: D_I1, D_I2, D_I3, closure | `e2_properties.py` | `raw_property_validation.csv` |
| Shot analysis: closed form against circuit | `e3_analytical_vs_circuit.py` | `analytical_vs_circuit.csv`, `shot_summary.csv` |
| Finite-shot uncertainty: interval coverage | `e4_coverage.py` | `finite_shot_coverage.csv` and raw |
| Multi-source, permutation, weighting, action layer | `e5_multisource.py` | `multi_source_validation.csv`, `weighting_example.csv` |
| Scalability and resources | `e6_runtime.py` | `runtime_scalability.csv` |

The protocol, covering seeds, shot budgets, grid resolutions, repetition
counts and tolerances, is documented in
[`tfs_evaluation/README.md`](medical-decision-backend/experiments/tfs_evaluation/README.md).

### Determinism

Every reported value is reproducible. The analytical script is deterministic.
The noise study fixes `seed_transpiler=1` and `seed_simulator` in 1..30, so
rerunning it on the same versions returns the same numbers; on other Qiskit
versions the exact values move within the reported standard deviations.

## Summary of results

| Aspect | Result |
| --- | --- |
| Operator properties | No violation of D_I1, D_I2, D_I3 or closure over 1,758,276 pairs of U~; largest deviation 2.2e-16 |
| Circuit against closed form | Mean absolute error 1.8e-3 at 32,000 shots; sqrt(N) * MAE remains within 0.29 to 0.32 across a sixteenfold range of N |
| Interval coverage | Mean 0.956 over 27 estimates from 1,800 repetitions, against a nominal 0.95 |
| Aggregation | P = C(n,2) exactly; permutation invariance to 1e-16; decision scores sum to 1 to 1e-16 |
| Scalability | About 58 ms per pair, constant from n = 2 to n = 12; classical fusion accounts for 0.0014% of the total |

## Scope and limitations

No quantum computational advantage is claimed or measured. The closed form
computes the same quantity exactly and in constant time, and the circuit is
slower and less accurate than the arithmetic it realises.

All circuit results come from classical simulation with the Qiskit Aer
`qasm_simulator`, not from quantum hardware, and no noise model is applied.

Resource figures transpiled only to a `{cx, u}` basis assume free connectivity
and are logical lower bounds rather than hardware costs. `resources.py` also
reports topology-routed estimates.

The opinion values used throughout are synthetic and carry no clinical
validity.

Numerical verification over finite grids is evidence, not proof. The proofs
are given in the manuscript.
