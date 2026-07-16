#!/usr/bin/env python3
"""Convert a tokamak EQDSK g-file into a SIMPLE-ready Boozer chartmap.

This drives the two documented Rung 0 steps in one command:

1. libneo's EQDSK-to-Boozer converter (``libneo.eqdsk_to_boozer_chartmap``,
   merged in itpplasma/libneo#346) produces an axisymmetric chartmap on the
   converter grid.
2. ``rung0/fix_chartmap_for_simple.py`` regrids it onto the axis-complete
   uniform grid SIMPLE expects (uniform rho in [0, 1], nzeta >= 2).

Any axisymmetric tokamak g-file works, not just the circular benchmark case.
The intermediate converter-grid file is kept next to the output for
provenance (suffix ``.converter.nc``).

Requirements: the libneo python package and its compiled ``_efit_to_boozer``
module must be importable. With the sibling layout of this repository::

    PYTHONPATH=../libneo/python:../SIMPLE/build/_deps/libneo-build \
        python tools/eqdsk_to_simple_chartmap.py input.eqdsk output.nc

Geometry in the output is in centimetres and the public radial coordinate is
``rho_tor = sqrt(s_tor)``, as documented in ``rung0/README.md``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("gfile", help="EQDSK g-file input")
    parser.add_argument("output", help="SIMPLE-ready chartmap NetCDF output")
    parser.add_argument("--nrho-converter", type=int, default=20,
                        help="radial surfaces on the converter grid")
    parser.add_argument("--ntheta", type=int, default=48,
                        help="poloidal Boozer points (endpoint-excluded)")
    parser.add_argument("--nrho-out", type=int, default=101,
                        help="uniform radial points in the SIMPLE grid")
    parser.add_argument("--nzeta-out", type=int, default=4,
                        help="toroidal planes in the SIMPLE grid (>= 2)")
    parser.add_argument("--n-s", type=int, default=101,
                        help="uniform s points for the A_phi profile")
    parser.add_argument("--nsurfmax", type=int, default=400,
                        help="starting points for the separatrix search "
                        "(libneo default 10000 is far slower than needed)")
    parser.add_argument("--nlabel", type=int, default=200,
                        help="internal radial grid of the field-line "
                        "integration (libneo default 500)")
    parser.add_argument("--no-psimax", action="store_true",
                        help="do not stop the flux-surface scan at the "
                        "g-file boundary flux (fall back to locating the "
                        "boundary purely by field-line integration)")
    args = parser.parse_args(argv)

    try:
        from libneo.eqdsk_to_boozer_chartmap import convert_eqdsk_to_chartmap
    except ImportError as exc:
        raise SystemExit(
            "libneo python package (with _efit_to_boozer) not importable; "
            "set PYTHONPATH as described in the module docstring"
        ) from exc

    sys.path.insert(0, str(REPO / "rung0"))
    from fix_chartmap_for_simple import fix_chartmap

    output = Path(args.output)
    converter_file = output.with_suffix(".converter.nc")

    psimax = None
    if not args.no_psimax:
        from libneo.eqdsk_base import read_eqdsk
        eq = read_eqdsk(args.gfile)
        psi_axis = float(eq["PsiaxisVs"])
        psi_edge = float(eq["PsiedgeVs"])
        if psi_edge > psi_axis:
            psimax = psi_edge * 1.0e8  # Wb/rad -> G cm^2/rad, absolute
            print(f"flux-surface scan stops at psimax={psimax:.6e} G cm^2/rad")
        else:
            print("g-file psi decreases outward; separatrix located by "
                  "field-line integration only")

    torflux = convert_eqdsk_to_chartmap(
        args.gfile,
        str(converter_file),
        nrho=args.nrho_converter,
        ntheta=args.ntheta,
        nsurfmax=args.nsurfmax,
        nlabel=args.nlabel,
        psimax=psimax,
    )
    print(f"converter chartmap: {converter_file} (torflux={torflux:.6e} G cm^2)")

    fix_chartmap(
        str(converter_file),
        str(output),
        nrho_out=args.nrho_out,
        nzeta_out=args.nzeta_out,
        n_s=args.n_s,
    )
    print(f"SIMPLE chartmap: {output}")


if __name__ == "__main__":
    main()
