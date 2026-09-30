"""E1 -- Experimental setup: environment, circuit and gate counts."""

import json
from collections import Counter

from _common import (OUT_DIR, PRODUCTION_DEFAULT_SHOTS, SEEDS, SHOT_BUDGETS,
                     GRID_STEP_CLOSURE, GRID_STEP_PROPERTIES, banner, environment)

import os

from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer

from qif_xnor import (BACKEND_NAME, N_MEASURE_OPS, N_PREP_OPS, STAGE_OPS,
                      build_xnor_circuit, figure_ops)


def main() -> dict:
    banner("E1  EXPERIMENTAL SETUP")
    env = environment()
    for k, v in env.items():
        print(f"  {k:22} {v}")

    qc = build_xnor_circuit((0.70, 0.20), (0.65, 0.25))
    fig = figure_ops(qc)
    full = Counter(ci.operation.name for ci in qc.data)
    logical = Counter(ci.operation.name for ci in fig)

    # depth of the T1--T11 logical block only
    only = QuantumCircuit(qc.num_qubits, qc.num_clbits)
    for ci in fig:
        only.append(ci.operation, ci.qubits, ci.clbits)

    backend = Aer.get_backend(BACKEND_NAME)
    cfg = backend.configuration()
    tc = transpile(qc, backend)

    circuit = {
        "qubits_total": qc.num_qubits,
        "qubits_input": 4,          # q0..q3
        "qubits_ancilla": 4,        # q4..q7
        "qubits_output": 2,         # q8, q9
        "clbits": qc.num_clbits,
        "state_preparation_gates": {"u": full["u"]},
        "logical_xnor_gates": {"x": logical["x"], "ccx": logical["ccx"]},
        "measurement_gates": full["measure"],
        "total_operations": len(qc.data),
        "depth_logical_block_only": only.depth(),
        "depth_full_circuit": qc.depth(),
        "transpiled_ops": dict(tc.count_ops()),
        "transpiled_depth": tc.depth(),
        "stages": [{"stage": s, "n_ops": n} for s, n in STAGE_OPS],
        "n_prep_ops": N_PREP_OPS,
        "n_measure_ops": N_MEASURE_OPS,
    }

    backend_info = {
        "name": backend.name,
        "class": type(backend).__name__,
        "simulator": bool(cfg.simulator),
        "method": str(backend.options.method) + " (None means automatic)",
        "coupling_map": cfg.coupling_map,
        "max_qubits": cfg.n_qubits,
        "n_basis_gates": len(cfg.basis_gates),
        "primitives_used": "backend.run (no Sampler/AerSampler primitive)",
        "transpile_call": "transpile(qc, backend)",
        "transpile_optimization_level": "default (not set explicitly)",
        "transpile_basis_gates": "backend defaults (not set explicitly)",
        "transpile_coupling_map": "backend defaults (None for qasm_simulator)",
        "seed_application": "seed_simulator on backend.run, one derived seed per pair",
    }

    protocol = {
        "production_default_shots": PRODUCTION_DEFAULT_SHOTS,
        "shot_budgets_evaluated": list(SHOT_BUDGETS),
        "seeds": SEEDS,
        "grid_step_properties": GRID_STEP_PROPERTIES,
        "grid_step_closure": GRID_STEP_CLOSURE,
        "seed_derivation": "numpy.random.SeedSequence(base).spawn(P), reduced mod 2**31-1",
    }

    print("\n  one pairwise circuit:")
    print(f"    qubits            {circuit['qubits_total']}"
          f"  (input {circuit['qubits_input']},"
          f" ancilla {circuit['qubits_ancilla']},"
          f" output {circuit['qubits_output']})")
    print(f"    (a) state prep    {circuit['state_preparation_gates']}")
    print(f"    (b) logical XNOR  {circuit['logical_xnor_gates']}"
          f"  in {len(STAGE_OPS)} stages T1--T11")
    print(f"    (c) measurement   {circuit['measurement_gates']}")
    print(f"    depth             logical block {circuit['depth_logical_block_only']}"
          f" | full circuit {circuit['depth_full_circuit']}"
          f" | transpiled {circuit['transpiled_depth']}")
    print(f"    transpiled        {circuit['transpiled_ops']}")

    setup = {"environment": env, "backend": backend_info,
             "circuit": circuit, "protocol": protocol}
    path = os.path.join(OUT_DIR, "experimental_setup.json")
    with open(path, "w") as f:
        json.dump(setup, f, indent=2)
    print(f"\n  -> {os.path.basename(path)}")
    return setup


if __name__ == "__main__":
    main()
