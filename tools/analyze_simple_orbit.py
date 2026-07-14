#!/usr/bin/env python3
"""Extract geometry and trapped-orbit frequencies from SIMPLE ``orbits.nc``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def _length_scale(units: str | None) -> float:
    normalized = (units or "").strip().lower()
    if normalized in {"m", "metre", "meter", "metres", "meters"}:
        return 1.0
    if normalized in {"cm", "centimetre", "centimeter", "centimetres", "centimeters"}:
        return 1.0e-2
    raise ValueError(f"unsupported or missing length unit {units!r}")


def _crossings(time: np.ndarray, values: np.ndarray, *fields: np.ndarray) -> tuple[np.ndarray, ...]:
    """Linearly interpolate zero crossings and accompanying fields."""
    indices = np.flatnonzero(np.signbit(values[:-1]) != np.signbit(values[1:]))
    fractions = -values[indices] / (values[indices + 1] - values[indices])
    result = [time[indices] + fractions * (time[indices + 1] - time[indices])]
    for field in fields:
        result.append(field[indices] + fractions * (field[indices + 1] - field[indices]))
    return tuple(result)


def analyze(path: Path, particle: int = 0) -> dict[str, object]:
    with Dataset(path) as dataset:
        variables = dataset.variables
        required = {"time", "R", "Z", "phi", "v_par"}
        missing = required.difference(variables)
        if missing:
            raise ValueError(f"{path} is missing variables: {', '.join(sorted(missing))}")

        time = np.asarray(variables["time"][:, particle], dtype=float)
        r = np.asarray(variables["R"][:, particle], dtype=float) * _length_scale(
            getattr(variables["R"], "units", None)
        )
        z = np.asarray(variables["Z"][:, particle], dtype=float) * _length_scale(
            getattr(variables["Z"], "units", None)
        )
        phi = np.asarray(variables["phi"][:, particle], dtype=float)
        v_parallel = np.asarray(variables["v_par"][:, particle], dtype=float)

    finite = np.isfinite(time) & np.isfinite(r) & np.isfinite(z) & np.isfinite(phi) & np.isfinite(v_parallel)
    time, r, z, phi, v_parallel = (array[finite] for array in (time, r, z, phi, v_parallel))
    if len(time) < 2 or np.any(np.diff(time) <= 0.0):
        raise ValueError("the finite trajectory must contain at least two strictly increasing times")

    phi = np.unwrap(phi)
    tip_time, tip_r, tip_z, tip_phi = _crossings(time, v_parallel, r, z, phi)
    result: dict[str, object] = {
        "input": str(path),
        "particle": particle,
        "orbit_points": len(time),
        "orbit_type": "trapped" if len(tip_time) else "passing",
        "R_min_m": float(np.min(r)),
        "R_max_m": float(np.max(r)),
        "Z_min_m": float(np.min(z)),
        "Z_max_m": float(np.max(z)),
        "radial_extent_m": float(np.ptp(r)),
        "mirror_crossings": len(tip_time),
        "mirror_times_s": tip_time.tolist(),
        "mirror_R_m": tip_r.tolist(),
        "mirror_Z_m": tip_z.tolist(),
    }

    if len(tip_time) >= 3:
        periods = tip_time[2:] - tip_time[:-2]
        angular_bounce = 2.0 * np.pi / periods
        angular_toroidal = (tip_phi[2:] - tip_phi[:-2]) / periods
        result.update(
            {
                "complete_bounce_periods": len(periods),
                "bounce_period_mean_s": float(np.mean(periods)),
                "bounce_period_std_s": float(np.std(periods)),
                "omega_b_mean_rad_s": float(np.mean(angular_bounce)),
                "omega_b_std_rad_s": float(np.std(angular_bounce)),
                "omega_phi_mean_rad_s": float(np.mean(angular_toroidal)),
                "omega_phi_std_rad_s": float(np.std(angular_toroidal)),
            }
        )
    else:
        result["complete_bounce_periods"] = 0
    return result


def _print_report(result: dict[str, object]) -> None:
    print(f"Orbit: {result['orbit_type']} ({result['orbit_points']} finite points)")
    print(
        "Geometry [m]: "
        f"R=[{result['R_min_m']:.6f}, {result['R_max_m']:.6f}], "
        f"Z=[{result['Z_min_m']:.6f}, {result['Z_max_m']:.6f}], "
        f"radial extent={result['radial_extent_m']:.6f}"
    )
    print(f"Mirror crossings: {result['mirror_crossings']}")
    if result["complete_bounce_periods"]:
        print(f"Complete bounce periods: {result['complete_bounce_periods']}")
        print(
            f"tau_b = {result['bounce_period_mean_s']:.8e} +/- "
            f"{result['bounce_period_std_s']:.2e} s"
        )
        print(
            f"omega_b = {result['omega_b_mean_rad_s']:.8e} +/- "
            f"{result['omega_b_std_rad_s']:.2e} rad/s"
        )
        print(
            f"omega_phi = {result['omega_phi_mean_rad_s']:.8e} +/- "
            f"{result['omega_phi_std_rad_s']:.2e} rad/s"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbits", type=Path, help="SIMPLE orbits.nc file")
    parser.add_argument("--particle", type=int, default=0, help="zero-based particle index")
    parser.add_argument("--json", type=Path, help="also write machine-readable diagnostics")
    args = parser.parse_args()

    result = analyze(args.orbits, args.particle)
    _print_report(result)
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
