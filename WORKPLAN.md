# Work plan: SIMPLE vs POTATO orbit benchmark

Companion to `doc/benchmark.tex`, which holds the physics, definitions, and
literature. This file holds the execution ladder, split into two tracks that
run in parallel:

- **Track A (Majid)**: run the codes as a user, produce the comparisons,
  review results and interfaces critically. Every stumbling block is a
  deliverable: file a GitHub issue on the code repo concerned
  ([SIMPLE](https://github.com/itpplasma/SIMPLE),
  [NEO-RT](https://github.com/itpplasma/NEO-RT), or this repo).
- **Track B (Chris)**: implement what Track A needs next; short iteration
  cycles.

Status flags used below: `[works]` runs today as described, `[interim]` works
with a temporary workaround, `[todo:B]` needs Track B implementation first.
No Track A step depends on an unimplemented feature.

Timeframe: July to end of August 2026. Weekly sync on rung status; rungs may
overlap once their inputs exist.

## Rung ladder

### Rung 0: run each code alone on its working case

Goal: both codes built and producing output on cases known to work, before
any cross-comparison.

SIMPLE `[works]`:

```sh
git clone https://github.com/itpplasma/SIMPLE && cd SIMPLE
cmake -S . -B build -G Ninja && cmake --build build -j"$(nproc)"
ctest --test-dir build
```

Then run the circular-tokamak case from this repo: copy `rung0/simple_run/`
somewhere writable and execute `simple.x` in it (`simple.in` and the chart
map are provided). Outputs: `times_lost.dat`, `confined_fraction.dat`,
`orbits.nc`. Plot confined fractions over time as a first sanity check.

POTATO `[works]`:

```sh
git clone https://github.com/itpplasma/NEO-RT && cd NEO-RT/POTATO
cmake -S . -B build -G Ninja && cmake --build build -j"$(nproc)"
ctest --test-dir build
```

Then run the circular case from this repo in a copy of `rung0/potato_run/`:

- `itest_type = 4` in `potato.in`: single orbit; trajectory streams to
  `fort.100` (columns `R phi Z p lambda dlambda/dtau`, NaN-line terminated),
  summary `taub`, `delphi` to stdout.
- `itest_type = 5`: radial frequency scan; writes `freq_scan.dat`
  (`R_start rho_pol omega_b omega_phi taub delphi ierr`).

Plot the single orbit in (R, Z) and `omega_b`, `omega_phi` against `rho_pol`.
Deliverable: short reproduction note (what ran, what the plots show, what was
unclear) plus issues for everything that needed guesswork.

Known broken, do not use: `examples/orbits_and_cuts.py` in SIMPLE fails at
import against the current Python interface `[todo:B]` (see Track B list).
Use `pysimple` (below) instead for Python-driven SIMPLE runs.

### Rung 1: same field in both codes, one orbit overlaid

Goal: verify field identity, then overlay one trapped orbit.

1. Field identity check `[interim]`: sample B and psi from the EQDSK side and
   the chart-map side along matched (R, Z) points; relative agreement at the
   interpolation-error level. Interim chart maps for the circular field are
   in `rung0/` with the scripts that produced them; the clean EQDSK-to-chart-map
   converter in libneo is `[todo:B]` and replaces them without changing this
   rung's procedure.
2. Single banana orbit, both codes, identical (rho_pol, lambda, E), reference
   deuteron at 5 keV. Convert SIMPLE output to (R, Z) via the coordinate
   chain in `pysimple`; overlay footprints; compare turning points and
   midplane width.
3. Invariant drift: energy and p_phi (psi*) drift per bounce for both codes,
   against step size and tolerance.

### Rung 2: frequency comparison at fixed surface

Goal: omega_b and omega_phi from three sources on one plot.

- POTATO `[works]`: `itest_type = 5` scan and the eta-scan diagnostic
  (`canonical_freqs_vs_eta_*.dat`).
- NEO-RT thin orbit `[works, manual]`: `test_freq_scan.f90` builds with the
  test suite and writes `freq_scan_neort.dat` on the same columns; registered
  as an automated test `[todo:B]`.
- SIMPLE `[todo:B]`: Python estimator over trajectory output (tip-crossing
  interpolation, as specified in doc/benchmark.tex section Metrics). Lands in
  `tools/` of this repo; Track A uses and reviews it, checks its
  sampling-resolution convergence.

Pitch scan through the trapped-passing boundary at fixed surface and energy;
points clustered toward the separatrix. Expect and document the logarithmic
tau_b divergence; the resonance-relevant combinations stay finite.

### Rung 3: orbit width and Poincare sections

Banana width against pitch and energy from (R, Z) footprints in both codes;
near-axis potato-class orbits compared explicitly. Tip and toroidal-angle
sections overlaid. SIMPLE-side section extraction goes through the same
trajectory post-processing as rung 2 `[todo:B]`; POTATO sections come from
its native first-return machinery `[works]`.

### Rung 4: sweeps and convergence

Full radial, pitch, and energy sweeps per doc/benchmark.tex; every reported
point backed by a step/tolerance halving scan. Energy sweep checks
convergence of both codes toward the NEO-RT thin-orbit values. Output: the
figure set for the final report.

### Rung 5: realistic equilibrium (stretch)

Repeat rungs 1-4 on an experimental EQDSK equilibrium. Requires non-public
equilibrium data; only after the circular case agrees.

## Track B task list (Chris)

Ordered by what Track A needs next.

1. `[libneo]` EQDSK to Boozer chart-map converter; replaces the interim
   `rung0/` chart maps (tracked as libneo issue 343).
2. `[this repo]` Python comparison package under `tools/`: loaders for
   `freq_scan.dat`, `fort.100`, `orbits.nc`, `freq_scan_neort.dat`; the
   SIMPLE frequency estimator (tip/period crossing interpolation); overlay
   and sweep plotting. Batch-first design: one call processes one run
   directory.
3. `[SIMPLE]` Fix the Python-interface example path: repair or replace
   `examples/orbits_and_cuts.py`; expose a Poincare-cut driver with plain
   scalar/array arguments so f90wrap can wrap it (the current
   `trace_to_cut` signature is not wrappable and silently disappears from
   the generated interface).
4. `[SIMPLE]` Native bounce/precession diagnostic reusing the existing
   tip-detection machinery, once the Python estimator has validated numbers
   (rung 2 exit criterion).
5. `[NEO-RT]` Register `test_freq_scan` as an automated test; align its
   output columns with POTATO's `freq_scan.dat`.
6. `[interfaces]` Keep the adapter layer in this repo until stable, then
   upstream: batch semantics (one call = N particles x full trace), sane
   defaults so a bare call on the circular field works, per-call rather than
   per-timestep language crossings.

## Working agreements

- Comparisons live in this repo; code changes live in the code repos and
  arrive by pull request.
- Every rung ends with a short written result note in `results/<rung>/`
  committed together with the plots and the exact input files.
- Disagreement between codes is reported with the resolution scan that shows
  it is not a discretisation artefact (doc/benchmark.tex, acceptance
  criteria).
