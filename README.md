# SIMPLE versus POTATO orbit benchmark

This repository compares collisionless guiding-centre orbits in axisymmetric
tokamak fields:

- [SIMPLE](https://github.com/itpplasma/SIMPLE): symplectic integration in
  canonical flux coordinates;
- [POTATO](https://github.com/itpplasma/NEO-RT): adaptive Runge-Kutta
  integration in cylindrical coordinates with finite orbit width;
- NEO-RT thin-orbit frequencies as a third reference where applicable.

The code implementations live in their own repositories. This repository owns
the common inputs, run procedure, cross-code analysis, plots, and benchmark
documentation.

## Start with Rung 0

Use sibling checkouts:

```text
../benchmark-simple-potato
../SIMPLE
../NEO-RT
```

Build SIMPLE and NEO-RT/POTATO, then follow
[`rung0/README.md`](rung0/README.md). It contains copy-paste run commands,
the input manifest, coordinate conventions, and Rung 0 acceptance checks.

The production equilibrium inputs are already tracked:

- `rung0/circ_chartmap_simple.nc` for SIMPLE (Boozer chart map);
- `rung0/wout_circ.nc` for SIMPLE's native VMEC input path;
- `rung0/potato_run/circ.eqdsk` for POTATO.

They describe the same public synthetic circular reactor-size equilibrium
(`R0 = 6.2 m`, `a = 2.0 m`, `B0 = 5.3 T`) through the codes' native input
formats. `tools/eqdsk_to_simple_chartmap.py` and `tools/eqdsk_to_vmec.py`
convert any axisymmetric tokamak g-file into the same two SIMPLE input
formats. The field-representation error must still be measured before
comparing orbits.

SIMPLE commit
[`08b85a1`](https://github.com/itpplasma/SIMPLE/commit/08b85a11279363693d6b1d5962772ace8df4ea45)
or newer writes direct cylindrical `R` and `Z` variables to `orbits.nc`. The
current chart map and POTATO case use centimetres. No canonical-to-VMEC
conversion is needed for trajectory overlays.

## Documentation

- [`WORKPLAN.md`](WORKPLAN.md): rung ladder, responsibilities, exit criteria,
  and current implementation status;
- [`doc/benchmark.tex`](doc/benchmark.tex): physics definitions, comparison
  metrics, estimators, and acceptance rules;
- [`rung0/README.md`](rung0/README.md): equilibrium provenance and executable
  Rung 0 procedure.

Build the working document from the repository root:

```sh
latexmk -pdf -cd doc/benchmark.tex
```

## Python environment

Create the shared environment after placing SIMPLE and NEO-RT in the sibling
layout:

```sh
./setup-venv.sh
source .venv/bin/activate
```

The setup script installs this repository's Python requirements, SIMPLE's
requirements and editable `pysimple` package, and NEO-RT's Python
requirements. On GNU/Linux it configures the SIMPLE extension for OpenBLAS so
it matches NumPy and SciPy in the environment.

Reusable post-processing scripts live in `tools/`. In particular,
`analyze_simple_orbit.py` extracts interpolated trapped-orbit tips and
frequencies from SIMPLE output, while `compare_orbits.py` overlays matched
SIMPLE and POTATO trajectories in physical cylindrical coordinates.
`potato_invariants_to_simple.py` validates POTATO's versioned invariant
handoff and returns every matching SIMPLE cut state without guessing an
outer/inner orbit class.

If SIMPLE was configured before the environment existed, reactivate `.venv`
and reinstall its Python package:

```sh
source .venv/bin/activate
python -m pip install --no-build-isolation -e ../SIMPLE
```

## Related benchmark

[`benchmark-orbit-proxima`](https://github.com/itpplasma/benchmark-orbit-proxima)
compares SIMPLE guiding-centre and full-orbit trajectories in the W7-X
high-mirror field against ASCOT5 and FIRM3D.
