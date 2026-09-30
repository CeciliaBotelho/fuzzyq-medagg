# FuzzyQ-MedAgg, Experimental Evaluation

Material for the experimental section of the manuscript. Every number reported
there is produced by the scripts in this directory.

These scripts **observe** the production modules (`qif_xnor.py`, `fusion.py`,
`resources.py`) and never modify them. No production behaviour is changed by
running them.

## Requirements

Run from inside this directory, with the backend dependencies installed:

```bash
cd medical-decision-backend
pip install -r requirements.txt
cd experiments/tfs_evaluation
```

`_common.py` prepends the backend directory to `sys.path`, so the production
modules are importable without installation.

## Reproducing everything

```bash
python3 run_all.py
```

Approximate wall-clock on the reference machine (Apple M-series, 10 cores):
E1 <1 s, E2 ~1 s, E3 ~47 s, E4 ~6 min, E5 ~11 s, E6 ~46 s.

## Which script produces which experiment

Each manuscript experiment maps to exactly one script. Four of the eight are
blocks inside `e5_multisource.py`, which prints them under the headings
`A)`--`F)`; there is no duplicated code.

| Manuscript experiment | Script | Block | Outputs |
|---|---|---|---|
| XNOR validation | `e2_properties.py` | all | `raw_property_validation.csv` |
| Shot analysis | `e3_analytical_vs_circuit.py` | all | `analytical_vs_circuit.csv`, `shot_summary.csv` |
| Finite-shot uncertainty | `e4_coverage.py` | all | `finite_shot_coverage.csv`, `finite_shot_coverage_raw.csv` |
| Multi-source validation | `e5_multisource.py` | `A` pair counts, `B` reduction at n=2 | `multi_source_validation.csv` |
| Permutation invariance | `e5_multisource.py` | `C` | `multi_source_validation.csv` |
| Weighting validation | `e5_multisource.py` | `D` normalisation, `E` hierarchical example | `multi_source_validation.csv`, `weighting_example.csv` |
| Action-layer validation | `e5_multisource.py` | `F` | `multi_source_validation.csv` |
| Scalability | `e6_runtime.py` | all | `runtime_scalability.csv` |
| Setup and circuit resources | `e1_setup.py` | all | `experimental_setup.json` |

To obtain only one of the blocks inside `e5`, run the whole script and read the
corresponding section of its output; the blocks share the opinion pool and the
seed offsets, so they are not separable without changing the protocol.

## Reproducing one experiment at a time

```bash
python3 e1_setup.py                 # environment, backend, circuit resources
python3 e2_properties.py            # D_I1, D_I2, D_I3, closure in U~
python3 e3_analytical_vs_circuit.py # closed form vs shot-based estimation
python3 e4_coverage.py              # empirical coverage of the 95% intervals
python3 e5_multisource.py           # blocks A-F: pair counts, n=2 reduction,
                                    # permutation invariance, weight
                                    # normalisation, hierarchical example,
                                    # action-layer validation
python3 e6_runtime.py               # gate counts, topologies, runtime scaling
```

`e4_coverage.py` is the long one. It buffers stdout when redirected, so use
`python3 -u e4_coverage.py` if you want live progress.

## Outputs

| File | Produced by | Contents |
|---|---|---|
| `experimental_setup.json` | E1 | versions, backend, circuit resources, protocol |
| `raw_property_validation.csv` | E2 | one row per property/domain, violations, max error |
| `analytical_vs_circuit.csv` | E3 | per case and shot budget, analytical vs estimated |
| `shot_summary.csv` | E3 | MAE, RMSE, max error and sqrt(N)·MAE per budget |
| `finite_shot_coverage.csv` | E4 | coverage and interval half-widths per cell |
| `finite_shot_coverage_raw.csv` | E4 | one row per repetition (raw data) |
| `multi_source_validation.csv` | E5 | blocks A–F with verdicts |
| `weighting_example.csv` | E5 | the 30 + 2 institution example |
| `runtime_scalability.csv` | E6 | circuit resources and timings per n |

## Protocol

Defined in `_common.py` and recorded in `experimental_setup.json`.

- **Backend** `Aer.get_backend("qasm_simulator")`, invoked through
  `backend.run`. No `Sampler` / `AerSampler` primitive is used. Simulation
  method is `automatic`.
- **Transpilation** `transpile(qc, backend)` with no explicit
  `optimization_level`, `basis_gates` or `coupling_map`. The simulator declares
  no coupling map, so production transpilation assumes free connectivity. The
  topology-qualified estimates in E6 come from `resources.py`, which transpiles
  to `{cx, u}` against explicit synthetic coupling maps with a fixed
  `seed_transpiler`.
- **Shot budgets** 2 000, 8 000 and 32 000. The production default is 20 000
  and is additionally reported in E3.
- **Seeds** all in `_common.SEEDS`. A base seed is passed to
  `fusion.aggregate`, which derives one independent stream per pair using
  `numpy.random.SeedSequence(base).spawn(P)`, reduced modulo `2**31 - 1`, and
  applies each as `seed_simulator` on `backend.run`.
  - E3: `20260901 + 100*case_index + budget_index`, plus `10000*(rep+1)` for
    the repeated-measures summary.
  - E4: `70000 + 1000*config_index + repetition`.
  - E5: `4242` plus block-specific offsets.
  - E6: `909 + repetition`.
  - Random valid opinions: `numpy.random.default_rng(13337 + ...)`, rejection
    sampling on `mu + nu <= 1`.
- **Grids over U~** step `0.05` (231 valid points) for the symmetry sweep, and
  step `0.02` (1326 valid points) for identity and closure.
- **Repetitions** E3 summary 30 per case and budget; E4 coverage 200 per cell;
  E5 decision layer 150 per `n`; E6 timings 5 per `n` for the sampled path and
  more for the cheap ones, reported as medians.
- **Tolerances** `1e-12` for floating-point equality, `1e-9` slack on
  `mu + nu <= 1`.

## What these experiments do not establish

No quantum computational advantage is claimed or measured. All results come
from classical simulation of the circuits, not from quantum hardware, and no
noise model is applied. The opinion values are synthetic and carry no clinical
validity.
