# Work plan: SIMPLE vs POTATO orbit benchmark

`doc/benchmark.tex` defines the physics and metrics. This file gives the
execution order, current status, and exit criterion for each rung.

- **Track A (Majid)**: run the codes as a user, produce the comparisons,
  review results and interfaces critically. Every stumbling block is a
  deliverable: file a GitHub issue on the code repo concerned
  ([SIMPLE](https://github.com/itpplasma/SIMPLE),
  [NEO-RT](https://github.com/itpplasma/NEO-RT), or this repo).
- **Track B (Chris)**: implement what Track A needs next; short iteration
  cycles.

Status flags: `[works]` runs as written, `[manual]` needs manual analysis,
and `[todo:B]` needs Track B implementation. No Track A step depends on a
`[todo:B]` item.

Timeframe: July to end of August 2026. Weekly sync on rung status; rungs may
overlap once their inputs exist.

## Rung ladder

### Rung 0: run each code alone

Goal: both codes built and producing output on cases known to work, before
any cross-comparison.

Build SIMPLE `[works]`:

```sh
cmake -S ../SIMPLE -B ../SIMPLE/build -G Ninja
cmake --build ../SIMPLE/build -j"$(nproc)"
ctest --test-dir ../SIMPLE/build
```

Build POTATO `[works]`:

```sh
cmake -S ../NEO-RT/POTATO -B ../NEO-RT/POTATO/build -G Ninja
cmake --build ../NEO-RT/POTATO/build -j"$(nproc)"
ctest --test-dir ../NEO-RT/POTATO/build
```

Return to this repository and follow the exact staging and run commands in
[`rung0/README.md`](rung0/README.md). Do not copy `rung0/simple_run/` by itself:
its chart map is a relative symlink. The manifest copies the final chart map
as a regular file into a fresh writable run directory.

SIMPLE outputs `times_lost.dat`, `confined_fraction.dat`, and `orbits.nc`.
SIMPLE commit `08b85a1` or newer writes direct cylindrical `R` and `Z` to
`orbits.nc`. For this chart map their units are centimetres. The legacy
variable `s` contains `rho_tor` because the chart-map radial coordinate is
`rho_tor = sqrt(s_tor)`.

POTATO modes used here:

- `itest_type = 4` in `potato.in`: single orbit; trajectory streams to
  `fort.100` (columns `R phi Z p xi dxi/dtau`, NaN-line terminated),
  summary `taub`, `delphi` to stdout.
- `itest_type = 5`: radial frequency scan; writes `freq_scan.dat`
  (`R_start rho_pol omega_b omega_phi taub delphi ierr`).
  `omega_b` and `omega_phi` are physical angular frequencies; the length-like
  `taub` must be divided by POTATO's reference speed `v0` to obtain seconds.

POTATO retains the input key `orbit_lambda`, but its value is the pitch cosine
`xi = v_parallel/v`. Reserve `Lambda = mu B0/E` for the analytic trapping
parameter in `doc/benchmark.tex`.

Exit: both programs run from fresh directories; SIMPLE has finite `R,Z` and
the trapped `xi = 0.3` orbit changes sign in `v_par`; POTATO writes a finite
single orbit. Record code commits and exact input files in the Rung 0 result
note. Also record output units and anything that required guesswork.

For Python-driven SIMPLE runs, use the `pysimple` interface from a local virtual
environment (`./setup-venv.sh` in this repo installs it from `../SIMPLE`).
`examples/orbits_and_cuts.py` is now a `pysimple.trace_orbit()` plotting example;
for the shortest supported entry point, start from `examples/simple_api.py`.

### Rung 1: measure field agreement and overlay one orbit

Goal: quantify the representation error between the EQDSK and chart map, then
overlay one trapped orbit.

1. Field representation check `[manual]`: sample `B` and flux from the EQDSK
   and chart-map sides at matched `(R,Z)` points. Report maximum and RMS
   differences. The chart maps in `rung0/` were produced by libneo's merged
   [EQDSK-to-Boozer converter](https://github.com/itpplasma/libneo/pull/346),
   then resampled for SIMPLE as documented in `rung0/README.md`.
2. Single banana orbit `[manual]`: use a deuteron with `E = 5 keV` and the same
   pitch cosine `xi` in both codes. Run SIMPLE first. Read the first finite
   `R,Z` sample from `orbits.nc` and use it as POTATO's `orbit_Rstart` and
   `orbit_Zstart`. Do not set `rho_tor` and `rho_pol` numerically equal; they
   label different fluxes. Overlay SIMPLE's direct `R,Z` variables with
   columns 1 and 3 of POTATO's `fort.100`.
3. Invariant drift: energy and p_phi (psi*) drift per bounce for both codes,
   against step size and tolerance.

Exit: the field mismatch is finite and documented; both trajectories start at
the same physical point; the overlay, turning-point differences, midplane
width difference, and invariant drifts are reported with a resolution scan.

### Rung 2: frequency comparison at fixed surface

Goal: omega_b and omega_phi from three sources on one plot.

- POTATO `[works]`: `itest_type = 5` scan and the eta-scan diagnostic
  (`canonical_freqs_vs_eta_*.dat`).
- NEO-RT thin orbit `[works, manual]`: `test_freq_scan.f90` builds with the
  test suite and writes `freq_scan_neort.dat` on the same columns; registered
  as an automated test `[todo:B]`.
- SIMPLE `[works, manual]`: `tools/analyze_simple_orbit.py`
  performs tip-crossing interpolation as specified in `doc/benchmark.tex`
  section Metrics. Track A checks its sampling-resolution convergence before
  frequency differences are accepted.

Scan the pitch cosine `xi` through the trapped-passing boundary at fixed
surface and energy, with points clustered toward the separatrix. Document the
logarithmic `tau_b` divergence; the resonance-relevant combinations stay
finite.

Exit: all three frequency sources use the same physical surface, particle, and
pitch; the SIMPLE estimator is converged in trajectory sampling; discrepancies
are compared with the field-representation and numerical-error floors.

### Rung 3: orbit width and Poincaré sections

Banana width against pitch and energy from (R, Z) footprints in both codes;
near-axis potato-class orbits compared explicitly. Tip and toroidal-angle
sections overlaid. SIMPLE-side section extraction goes through the same
trajectory post-processing as rung 2 `[todo:B]`; POTATO sections come from
its native first-return machinery `[works]`.

Exit: widths and section locations are reported with matched initial
conditions and a resolution scan; near-axis potato-class cases are labelled
separately.

### Rung 4: sweeps and convergence

Full radial, pitch, and energy sweeps per doc/benchmark.tex; every reported
point backed by a step/tolerance halving scan. Energy sweep checks
convergence of both codes toward the NEO-RT thin-orbit values. Output: the
figure set for the final report.

Exit: every plotted point carries its input files and convergence evidence;
failed or lost trajectories remain visible rather than being dropped.

### Rung 5: realistic equilibrium (stretch)

Repeat rungs 1-4 on an experimental EQDSK equilibrium. Requires non-public
equilibrium data; only after the circular case agrees.

Exit: the equilibrium has an agreed sharing policy and provenance, and the
circular-case acceptance checks still pass with the same analysis code.

## Track B task list (Chris)

Ordered by what Track A needs next.

1. `[done:SIMPLE]` Direct per-step `R,Z` output in `orbits.nc` through the
   active reference geometry (SIMPLE commit `08b85a1`).
2. `[done:libneo]` EQDSK-to-Boozer chart-map converter (libneo PR 346,
   closing issue 343).
3. `[this repo]` Python comparison package under `tools/`: loaders for
   `freq_scan.dat`, `fort.100`, `orbits.nc`, `freq_scan_neort.dat`; the
   SIMPLE frequency estimator (tip/period crossing interpolation); overlay
   and sweep plotting. Batch-first design: one call processes one run
   directory.
4. `[SIMPLE]` Keep the repaired `examples/orbits_and_cuts.py` on the supported
   `pysimple.trace_orbit()` path, and still expose a dedicated Poincaré-cut
   driver with plain scalar/array arguments so f90wrap can wrap it (the current
   `trace_to_cut` signature is not wrappable and silently disappears from the
   generated interface).
5. `[SIMPLE]` Native bounce/precession diagnostic reusing the existing
   tip-detection machinery, once the Python estimator has validated numbers
   (rung 2 exit criterion).
6. `[NEO-RT]` Register `test_freq_scan` as an automated test; align its
   output columns with POTATO's `freq_scan.dat`.
7. `[interfaces]` Keep the adapter layer in this repo until stable, then
   upstream: batch semantics (one call = N particles x full trace), defaults
   that run the circular field without extra configuration, and per-call
   rather than per-timestep language crossings.

## Working agreements

- Comparisons live in this repo; code changes live in the code repos and
  arrive by pull request.
- Every rung ends with a short written result note in `results/<rung>/`
  committed together with the plots and the exact input files.
- Disagreement between codes is reported with the resolution scan that shows
  it is not a discretisation artefact (doc/benchmark.tex, acceptance
  criteria).
