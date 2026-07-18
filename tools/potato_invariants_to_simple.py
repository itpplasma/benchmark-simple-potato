#!/usr/bin/env python3
"""Convert POTATO invariant rows into all matching SIMPLE initial states."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np


HANDOFF_VERSION = "# potato-simple-invariants-v1"
HANDOFF_COLUMNS = (
    "id",
    "R_cm",
    "Z_cm",
    "p",
    "xi",
    "sigma",
    "H0",
    "J_perp",
    "psi_star",
    "psi_axis",
    "psi_edge",
    "phi_elec",
    "v0_cm_s",
    "status",
)
INTEGER_COLUMNS = {"id", "sigma", "status"}
POTENTIAL_TOLERANCE = 1.0e-12
NORMALIZATION_TOLERANCE = 1.0e-12


def _parse_row(line: str, line_number: int) -> dict[str, float | int]:
    fields = line.split()
    if len(fields) != len(HANDOFF_COLUMNS):
        raise ValueError(f"line {line_number}: expected {len(HANDOFF_COLUMNS)} fields")
    values = [float(field) for field in fields]
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"line {line_number}: values must be finite")
    row = dict(zip(HANDOFF_COLUMNS, values))
    for name in INTEGER_COLUMNS:
        if not row[name].is_integer():
            raise ValueError(f"line {line_number}: {name} must be an integer")
        row[name] = int(row[name])
    _validate_row(row, line_number)
    return row


def _validate_row(row: dict[str, float | int], line_number: int) -> None:
    if row["status"] != 0:
        raise ValueError(f"line {line_number}: POTATO status is {row['status']}")
    if abs(row["phi_elec"]) > POTENTIAL_TOLERANCE:
        raise ValueError(f"line {line_number}: nonzero electrostatic potential")
    if row["p"] <= 0.0 or row["v0_cm_s"] <= 0.0:
        raise ValueError(f"line {line_number}: p and v0 must be positive")
    if abs(row["xi"]) > 1.0:
        raise ValueError(f"line {line_number}: xi is outside [-1, 1]")
    if row["psi_axis"] == row["psi_edge"]:
        raise ValueError(f"line {line_number}: axis and edge flux are equal")
    expected_sigma = int(np.sign(row["xi"]))
    if row["sigma"] != expected_sigma:
        raise ValueError(f"line {line_number}: sigma is inconsistent with xi")
    if not math.isclose(
        row["H0"], row["p"] ** 2, rel_tol=NORMALIZATION_TOLERANCE, abs_tol=1.0e-14
    ):
        raise ValueError(f"line {line_number}: H0 is not magnetic-only p^2")


def read_handoff(path: Path) -> list[dict[str, float | int]]:
    """Read and strictly validate a POTATO invariant handoff."""
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 3:
        raise ValueError("invariant handoff must contain a header and at least one row")
    if lines[0] != HANDOFF_VERSION:
        raise ValueError(f"unsupported invariant handoff version: {lines[0]}")
    columns = tuple(lines[1].removeprefix("# ").split())
    if columns != HANDOFF_COLUMNS:
        raise ValueError(f"unexpected invariant columns: {columns}")
    rows = [_parse_row(line, index) for index, line in enumerate(lines[2:], start=3)]
    identifiers = [row["id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("invariant identifiers must be unique")
    return rows


def _candidate(result: dict, index: int, inspector: Callable | None) -> dict:
    state = np.asarray(result["states"][index], dtype=float)
    cylindrical = np.asarray(result["cylindrical"][index], dtype=float)
    if state.shape != (5,) or cylindrical.shape != (3,):
        raise ValueError("SIMPLE returned an unexpected candidate shape")
    if not np.isfinite(state).all() or not np.isfinite(cylindrical).all():
        raise ValueError("SIMPLE returned a non-finite candidate")
    candidate = {
        "state": state.tolist(),
        "cylindrical_native": cylindrical.tolist(),
        "sigma": int(result["sigma"][index]),
        "section_branch": int(result["section_branch"][index]),
        "section_kind_code": int(result["section_kind_code"][index]),
        "section_kind": str(result["section_kind"][index]),
        "residual": float(result["residual"][index]),
    }
    if inspector is not None:
        reconstructed = inspector(state)
        expected = np.asarray(result["simple_invariants"], dtype=float)
        actual = np.asarray(
            [reconstructed["h0"], reconstructed["j_perp"], reconstructed["p_phi"]]
        )
        if reconstructed["status"] != 0:
            raise ValueError("SIMPLE failed to reconstruct a candidate")
        tolerances = ((1.0e-10, 1.0e-12), (1.0e-10, 1.0e-12), (1.0e-9, 1.0e-9))
        if any(
            not math.isclose(got, want, rel_tol=rtol, abs_tol=atol)
            for got, want, (rtol, atol) in zip(actual, expected, tolerances)
        ):
            raise ValueError("SIMPLE candidate does not reconstruct its invariants")
        candidate["reconstructed_invariants"] = {
            "h0": float(actual[0]),
            "j_perp": float(actual[1]),
            "p_phi": float(actual[2]),
        }
    return candidate


def _source(row: dict[str, float | int]) -> dict:
    source = dict(row)
    source["p_phi_convention"] = "raw POTATO psi_star=(c/q)P_phi"
    source["flux_units"] = "POTATO native G cm^2"
    return source


def convert_rows(
    rows: list[dict[str, float | int]],
    resolver: Callable,
    *,
    simple_v0_cm_s: float,
    allow_v0_rescale: bool = False,
    max_solutions: int = 64,
    inspector: Callable | None = None,
) -> dict:
    """Resolve all SIMPLE cut states without selecting an orbit topology."""
    if not math.isfinite(simple_v0_cm_s) or simple_v0_cm_s <= 0.0:
        raise ValueError("SIMPLE v0 must be finite and positive")
    converted_rows = []
    for row in rows:
        ratio = row["v0_cm_s"] / simple_v0_cm_s
        if not allow_v0_rescale and not math.isclose(
            ratio, 1.0, rel_tol=NORMALIZATION_TOLERANCE, abs_tol=0.0
        ):
            raise ValueError(
                f"row {row['id']}: reference velocity mismatch; "
                "use --allow-v0-rescale explicitly"
            )
        result = resolver(
            row["H0"],
            row["J_perp"],
            row["psi_star"],
            psi_axis=row["psi_axis"],
            psi_edge=row["psi_edge"],
            v0_ratio=ratio,
            max_solutions=max_solutions,
        )
        if result["status"] != 0 or result["n_solutions"] < 1:
            raise ValueError(
                f"row {row['id']}: SIMPLE inversion failed with status "
                f"{result['status']} and {result['n_solutions']} solutions"
            )
        converted_rows.append(
            {
                "source": _source(row),
                "v0_ratio_potato_over_simple": ratio,
                "simple_invariants": result["simple_invariants"].tolist(),
                "ordering": result["ordering"],
                "candidates": [
                    _candidate(result, index, inspector)
                    for index in range(result["n_solutions"])
                ],
            }
        )
    return {
        "schema": "potato-simple-candidates-v1",
        "time_convention": "tau=v0*t",
        "simple_v0_cm_s": simple_v0_cm_s,
        "topology_policy": "all regular cut roots retained; no outer/inner guess",
        "rows": converted_rows,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("handoff", type=Path, help="POTATO potato_invariants.dat")
    parser.add_argument("equilibrium", type=Path, help="SIMPLE equilibrium input")
    parser.add_argument("--output", type=Path, default=Path("simple_candidates.json"))
    parser.add_argument("--n-e", type=int, default=1, help="particle charge number")
    parser.add_argument("--n-d", type=int, default=2, help="particle mass number")
    parser.add_argument("--facE-al", type=float, default=700.0)
    parser.add_argument("--max-solutions", type=int, default=64)
    parser.add_argument(
        "--allow-v0-rescale",
        action="store_true",
        help="explicitly rescale H0 and J_perp when POTATO and SIMPLE v0 differ",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    import pysimple

    rows = read_handoff(args.handoff)
    pysimple.init(
        str(args.equilibrium),
        deterministic=True,
        ntestpart=1,
        n_e=args.n_e,
        n_d=args.n_d,
        facE_al=args.facE_al,
    )
    output = convert_rows(
        rows,
        pysimple.states_from_potato_invariants,
        simple_v0_cm_s=float(pysimple.params.v0),
        allow_v0_rescale=args.allow_v0_rescale,
        max_solutions=args.max_solutions,
        inspector=pysimple.invariants_from_state,
    )
    output["equilibrium"] = str(args.equilibrium)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"wrote {sum(len(row['candidates']) for row in output['rows'])} candidates to {args.output}")


if __name__ == "__main__":
    main()
