#!/usr/bin/env python3
"""Overlay one SIMPLE and one POTATO trajectory in the cylindrical R-Z plane."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


def _to_metres(values: np.ndarray, units: str) -> np.ndarray:
    if units == "m":
        return values
    if units == "cm":
        return values * 1.0e-2
    raise ValueError(f"unsupported length unit {units!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("simple", type=Path, help="SIMPLE orbits.nc")
    parser.add_argument("potato", type=Path, help="POTATO fort.100")
    parser.add_argument("--particle", type=int, default=0)
    parser.add_argument("--potato-units", choices=("cm", "m"), default="cm")
    parser.add_argument("--output", type=Path, default=Path("orbit-overlay.png"))
    args = parser.parse_args()

    with Dataset(args.simple) as dataset:
        r_var, z_var = dataset.variables["R"], dataset.variables["Z"]
        r_simple = _to_metres(np.asarray(r_var[:, args.particle], dtype=float), r_var.units)
        z_simple = _to_metres(np.asarray(z_var[:, args.particle], dtype=float), z_var.units)
    simple_finite = np.isfinite(r_simple) & np.isfinite(z_simple)
    r_simple, z_simple = r_simple[simple_finite], z_simple[simple_finite]

    potato = np.genfromtxt(args.potato)
    r_potato = _to_metres(potato[:, 0], args.potato_units)
    z_potato = _to_metres(potato[:, 2], args.potato_units)
    potato_finite = np.isfinite(r_potato) & np.isfinite(z_potato)
    r_potato, z_potato = r_potato[potato_finite], z_potato[potato_finite]

    if not len(r_simple) or not len(r_potato):
        raise ValueError("both trajectories must contain finite R,Z samples")

    figure, axes = plt.subplots(figsize=(7, 7))
    axes.plot(r_simple, z_simple, label="SIMPLE")
    axes.plot(r_potato, z_potato, "--", label="POTATO")
    axes.scatter(r_simple[0], z_simple[0], marker="o", label="SIMPLE start")
    axes.scatter(r_potato[0], z_potato[0], marker="s", label="POTATO start")
    axes.set(xlabel="R [m]", ylabel="Z [m]", title="Guiding-centre orbit comparison")
    axes.axis("equal")
    axes.grid(True)
    axes.legend()
    figure.tight_layout()
    figure.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
