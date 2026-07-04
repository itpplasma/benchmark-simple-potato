# benchmark-simple-potato

Benchmark of collisionless guiding-centre orbits in axisymmetric tokamak
fields: [SIMPLE](https://github.com/itpplasma/SIMPLE) (symplectic, canonical
flux coordinates) against [POTATO](https://github.com/itpplasma/NEO-RT)
(adaptive RK, cylindrical coordinates, finite orbit width), with NEO-RT
thin-orbit frequencies as third reference.

- `doc/benchmark.tex`: physics, definitions, metrics, literature
  (build: `latexmk -pdf benchmark.tex` in `doc/`).
- `WORKPLAN.md`: execution ladder, current status, task split.
- `rung0/`: common circular-tokamak field inputs for both codes.

## Python environment for the benchmark workflow

This repository can host the recommended local virtual environment for the
benchmark and the related SIMPLE/NEO-RT Python tooling. From the repository
root run:

```sh
./setup-venv.sh
```

By default this script uses the sibling checkouts in `../SIMPLE` and
`../NEO-RT`, installs their Python dependencies into `.venv/`, and installs
SIMPLE's `pysimple` package in editable mode. This avoids the system Python on
PEP 668-managed machines.

Expected checkout layout:

```text
../benchmark-simple-potato
../SIMPLE
../NEO-RT
```

Later, reactivate the same environment with:

```sh
source .venv/bin/activate
```

If SIMPLE previously configured without Python support and printed `Python
f90wrap not found, skipping interface build.`, activate this `.venv` first and
then rerun the SIMPLE Python install inside it:

```sh
source .venv/bin/activate
python -m pip install --no-build-isolation -e ../SIMPLE
```
