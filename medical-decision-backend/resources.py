"""
Circuit resource estimation
===========================

Transpiling only to a gate basis ({cx, u}) measures the LOGICAL cost of the
decomposition and assumes full connectivity. Real superconducting hardware has
sparse connectivity, and routing inserts SWAPs that increase depth and CNOT
count -- often by a large factor.

This module separates the two explicitly:

  logical   {cx, u} basis, free connectivity. This is a LOWER BOUND, useful
            for comparing constructions, and must not be published as a
            hardware cost.
  linear    open path. Realistic worst case for connectivity.
  ring      cycle.
  grid      square lattice.
  heavy-hex topology of IBM superconducting devices.

Routing is stochastic; seed_transpiler fixes the outcome.
"""

import math
from typing import Dict, List, Optional

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

__all__ = ["TOPOLOGIES", "coupling_for", "estimate", "estimate_all"]

BASIS = ["cx", "u"]
TOPOLOGIES = ("logical", "linear", "ring", "grid", "heavy-hex")


def coupling_for(topology: str, n_qubits: int) -> Optional[CouplingMap]:
    """Smallest coupling map of the requested kind that fits n_qubits."""
    if topology == "logical":
        return None
    if topology == "linear":
        return CouplingMap.from_line(n_qubits)
    if topology == "ring":
        return CouplingMap.from_ring(n_qubits)
    if topology == "grid":
        side = math.ceil(math.sqrt(n_qubits))
        return CouplingMap.from_grid(side, side)
    if topology == "heavy-hex":
        d = 3
        while True:                      # odd distances; grow until it fits
            cm = CouplingMap.from_heavy_hex(d)
            if cm.size() >= n_qubits:
                return cm
            d += 2
    raise ValueError(f"unknown topology: {topology!r}")


def estimate(
    qc: QuantumCircuit,
    topology: str = "logical",
    seed: int = 2026,
    optimization_level: int = 1,
) -> Dict[str, object]:
    """Depth and CNOT count after decomposition and routing."""
    cm = coupling_for(topology, qc.num_qubits)
    tc = transpile(
        qc,
        basis_gates=BASIS,
        coupling_map=cm,
        optimization_level=optimization_level,
        seed_transpiler=seed,
    )
    return {
        "topology": topology,
        "device_qubits": None if cm is None else cm.size(),
        "circuit_qubits": qc.num_qubits,
        "depth": tc.depth(),
        "cx": tc.count_ops().get("cx", 0),
        "seed_transpiler": seed,
        "optimization_level": optimization_level,
    }


def estimate_all(
    qc: QuantumCircuit,
    topologies: Optional[List[str]] = None,
    seed: int = 2026,
) -> List[Dict[str, object]]:
    return [estimate(qc, t, seed=seed) for t in (topologies or list(TOPOLOGIES))]
