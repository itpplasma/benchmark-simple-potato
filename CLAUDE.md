# benchmark-simple-potato

This repository is the **comparison and workflow layer** for orbit benchmarks
between the sibling checkouts:

- `../SIMPLE`
- `../NEO-RT`

It is not the primary home for SIMPLE or NEO-RT implementation work. Use this
repo for:

- benchmark orchestration,
- shared Python environment setup,
- cross-code comparison scripts,
- benchmark inputs and run directories,
- plots, notes, and benchmark documentation.

## Repository purpose

The benchmark compares collisionless guiding-centre orbits in axisymmetric
tokamak fields:

- **SIMPLE**: symplectic, canonical-flux coordinates
- **POTATO / NEO-RT**: adaptive RK, cylindrical coordinates, finite-orbit-width
- **NEO-RT thin-orbit frequencies**: third reference where applicable

## Expected sibling layout

This repository assumes the following checkout structure:

```text
../benchmark-simple-potato
../SIMPLE
../NEO-RT
```

Many scripts and instructions assume those sibling paths exist.

## Shared Python environment

The recommended shared virtual environment for the whole benchmark workspace
lives in this repository:

```text
./.venv
```

Create it from the repository root with:

```bash
./setup-venv.sh
```

That script:

- creates or reuses `./.venv`,
- installs this repo's Python requirements,
- installs SIMPLE's Python requirements,
- installs `../SIMPLE` in editable mode,
- installs NEO-RT's Python requirements.

Reactivate later with:

```bash
source .venv/bin/activate
```

If the SIMPLE Python extension was previously built without Python support,
reactivate `.venv` first and reinstall SIMPLE from there.

## Where the important content lives

- `README.md` - top-level benchmark summary and shared-venv notes
- `WORKPLAN.md` - rung ladder, task split, current benchmark workflow
- `doc/benchmark.tex` - physics definitions, metrics, literature, comparison spec
- `rung0/` - common circular tokamak inputs and example run directories
- `doc/` - benchmark documentation sources

## Where to look next

### If the task is about the benchmark workflow

Start with:

1. `README.md`
2. `WORKPLAN.md`
3. `doc/benchmark.tex`

### If the task turns into a SIMPLE change

Switch to the sibling repo and read:

1. `../SIMPLE/CLAUDE.md`
2. `../SIMPLE/README.md`
3. `../SIMPLE/python/README.md`
4. `../SIMPLE/DOC/coordinates-and-fields.md`

### If the task turns into a NEO-RT change

Switch to the sibling repo and read:

1. `../NEO-RT/CLAUDE.md`
2. `../NEO-RT/README.md`
3. `../NEO-RT/doc/running.md`
4. `../NEO-RT/doc/file_formats.md`
5. `../NEO-RT/doc/library.md`

## Typical benchmark tasks

- set up or refresh the shared `.venv`,
- run or post-process SIMPLE and POTATO benchmark cases,
- compare trajectories, frequencies, or confinement outputs,
- update `WORKPLAN.md`,
- update `doc/benchmark.tex`,
- add benchmark-side tools or plotting scripts.

## Commands commonly relevant here

### Shared Python environment

```bash
./setup-venv.sh
source .venv/bin/activate
```

### Benchmark document

```bash
latexmk -pdf doc/benchmark.tex
```

Run that from the repository root or from `doc/`, depending on your preferred
LaTeX workflow.

## Working rules for this repo

- Treat this repo as the **integration / comparison** layer.
- Keep code changes for SIMPLE or NEO-RT in their own repositories unless the
  task is clearly benchmark-only.
- Keep notes, plots, and scripts here aligned with the current rung structure in
  `WORKPLAN.md`.
- When documenting procedures, prefer sibling-relative paths (`../SIMPLE`,
  `../NEO-RT`) because that is the expected workspace layout.
- Keep the repository and `.venv` group read/writeable in the shared workspace.

## Maintenance

Update this file when any of the following change:

- sibling checkout assumptions,
- shared `.venv` behavior,
- main benchmark doc locations,
- benchmark workflow entrypoints,
- the relationship between this repo and the sibling SIMPLE / NEO-RT repos.
