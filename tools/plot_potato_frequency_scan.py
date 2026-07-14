#!/usr/bin/env python3
"""Plot POTATO's ``freq_scan.dat`` while retaining failed scan points."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scan", type=Path, help="POTATO freq_scan.dat")
    parser.add_argument("--output", type=Path, default=Path("potato-frequency-scan.png"))
    args = parser.parse_args()

    data = np.loadtxt(args.scan)
    data = np.atleast_2d(data)
    rho, omega_b, omega_phi = data[:, 1], data[:, 2], data[:, 3]
    good = data[:, 6].astype(int) == 0

    figure, axes = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    for axis, values, label in zip(
        axes, (omega_b, omega_phi), (r"$\omega_b$ [rad/s]", r"$\omega_\varphi$ [rad/s]")
    ):
        axis.plot(rho[good], values[good], "o-", label="ierr = 0")
        if np.any(~good):
            axis.plot(rho[~good], values[~good], "rx", label="failed")
        axis.set_ylabel(label)
        axis.grid(True)
        axis.legend()
    axes[-1].set_xlabel(r"$\rho_{\rm pol}$")
    figure.tight_layout()
    figure.savefig(args.output, dpi=200)


if __name__ == "__main__":
    main()
