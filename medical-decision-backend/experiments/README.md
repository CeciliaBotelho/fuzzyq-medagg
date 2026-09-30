# Experiments

This directory holds two categories of scripts. They are **not**
interchangeable, and only the first reproduces the results reported in the
manuscript.

## Main experiments of the manuscript

**`tfs_evaluation/`**, everything reported in the experimental section. See
[`tfs_evaluation/README.md`](tfs_evaluation/README.md) for the protocol
(seeds, shot budgets, grid resolutions, repetition counts, tolerances) and for
the mapping from each manuscript experiment to its script.

```bash
cd tfs_evaluation && python3 run_all.py
```

## Additional analyses (not part of the main experiments)

These two scripts support design decisions discussed in the text. They are
kept for completeness and are **not** required to reproduce any reported
result. Neither modifies production behaviour.

| Script | Kind | Purpose |
|---|---|---|
| `compare_registers.py` | design-space comparison | Compares the pairwise register against an indexed alternative that samples the mean over all pairs in a single execution. Shows no accuracy gain at an equal total shot budget while width and depth both increase. Exercises `indexed.py`. |
| `compare_decision_rules.py` | ablation | Places the previous decision rule (`q_conf` / hesitation heuristics) side by side with the current score-based decision layer, separating the effect of the circuit correction from the effect of the rule change. Reimplements the previous rule locally; `fusion.py` is untouched. |

```bash
python3 compare_registers.py         # design-space comparison
python3 compare_decision_rules.py    # decision-rule ablation
```

`compare_registers.py` writes a LaTeX table next to itself. That artifact is
derived output and is not tracked; regenerate it if needed.
