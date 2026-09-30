# FuzzyQ-MedAgg

Agreement-based aggregation of intuitionistic fuzzy opinions through a quantum
XNOR circuit.

This repository accompanies a manuscript under double-blind review. Author,
affiliation and institutional information have been removed.

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
intuitionistic fuzzy value x̃ = (μ, ν): belief for, belief
against, and the hesitation π = 1 − μ − ν left over. Agreement between
two sources is measured by an intuitionistic fuzzy XNOR evaluated on a quantum
circuit, and the aggregate of n sources is mapped onto one of three clinical
actions: treat, do not treat, or request further examinations.

Diagnostic responsibility remains with the human specialist. The framework
provides a structured and auditable mechanism for handling uncertainty and
conflict among multiple AI predictions.

## Method

### Agreement operator

With the product t-norm T_P(a,b) = ab and the probabilistic sum
S_P(a,b) = a + b − ab:

$$
\mu_{\boxplus_{I_P}} = T_P\big(S_P(x_2,y_1),\; S_P(x_1,y_2)\big)
$$

$$
\nu_{\boxplus_{I_P}} = S_P\big(T_P(x_1,y_2),\; T_P(x_2,y_1)\big)
$$

The XNOR is the dual of the XOR under the intuitionistic negation
N_IS(x₁,x₂) = (x₂,x₁), which transposes components. It is not the
fuzzy negation 1 − x: that pair would leave Ũ, violating
μ + ν ≤ 1, and would break the identity axiom D_I3.

### Quantum circuit

Stages T₁ to T₁₁ reproduce the manuscript figure gate for gate.
Amplitude encoding precedes T₁ and measurement follows T₁₁; neither
belongs to the figure.

| Property | Value |
| --- | --- |
| Qubits | 10 (4 input, 4 ancilla, 2 output) |
| Logical gates | 15 Pauli-X, 6 Toffoli |
| Depth | 7 for the logical block, 9 for the full circuit |
| Outputs | q₈ = μ_⊞, q₉ = ν_⊞ |

Width and depth are fixed and do not depend on the number of sources.

### Aggregating n sources

⊞ is a binary equivalence and is non-associative, so no fold is
available. Agreement is formed from the P = C(n,2) pairs:

$$
\tilde{C} = \left( \sum_{i \lt j} w_{ij}\, \mu_{ij}, \;
\sum_{i \lt j} w_{ij}\, \nu_{ij} \right)
\qquad \text{with} \qquad \sum_{i \lt j} w_{ij} = 1
$$

A convex combination of points of the convex set Ũ remains in
Ũ, so π_C ≥ 0 requires no truncation. The weights depend only
on the unordered pair, which makes the result independent of presentation
order. A single pipeline covers every n ≥ 2; for n = 2 it reduces
exactly to ⊞(x̃₁, x̃₂) with w₁₂ = 1.

Three weighting schemes are available. The uniform scheme sets
w_ij = 1/P. The weighted scheme uses per-source weights. The hierarchical
scheme sets ω_i = 1/(K n_k), balancing the source weights across the
K institutions before the pair weights are induced; this attenuates the
dominance of a larger site without making the final contributions equal.

### From degrees to a clinical action

Proposition 2 gives μ_C + ν_C + π_C = 1 exactly, so the three degrees
already partition the unit. The score s(Ẽ) = μ_E − ν_E, the same
function used in the Xu-Yager ordering, supplies the direction of the mean
evidence:

$$
\text{treat} = \mu_C \frac{1+s}{2}
\qquad
\text{doNotTreat} = \mu_C \frac{1-s}{2}
\qquad
\text{requestExams} = \nu_C + \pi_C
$$

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

Every number reported in the manuscript comes from
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

## Summary of results

| Aspect | Result |
| --- | --- |
| Operator properties | No violation of D_I1, D_I2, D_I3 or closure over 1,758,276 pairs of Ũ; largest deviation 2.2 × 10⁻¹⁶ |
| Circuit against closed form | Mean absolute error 1.8 × 10⁻³ at 32,000 shots; √N · MAE remains within 0.29 to 0.32 across a sixteenfold range of N |
| Interval coverage | Mean 0.956 over 27 estimates from 1,800 repetitions, against a nominal 0.95 |
| Aggregation | P = C(n,2) exactly; permutation invariance to 10⁻¹⁶; decision scores sum to 1 to 10⁻¹⁶ |
| Scalability | About 58 ms per pair, constant from n = 2 to n = 12; classical fusion accounts for 0.0014% of the total |

## Scope and limitations

No quantum computational advantage is claimed or measured. The closed form of
Proposition 2 computes the same quantity exactly and in constant time, and the
circuit is slower and less accurate than the arithmetic it realises.

All circuit results come from classical simulation with the Qiskit Aer
`qasm_simulator`, not from quantum hardware, and no noise model is applied.

Resource figures transpiled only to a `{cx, u}` basis assume free connectivity
and are logical lower bounds rather than hardware costs. `resources.py` also
reports topology-routed estimates.

The opinion values used throughout are synthetic and carry no clinical
validity.

Numerical verification over finite grids is evidence, not proof. The proofs
are given in the manuscript.
