# Rung 0 inputs and runs

Rung 0 runs SIMPLE and POTATO separately on the public synthetic circular
tokamak equilibrium. Run both cases before changing particle parameters or
comparing trajectories.

## Equilibrium files

| Path | Purpose |
|---|---|
| `potato_run/circ.eqdsk` | POTATO equilibrium: synthetic 129 x 129 EQDSK with `R0 = 6.20 m`, `a = 2.00 m`, and `B0 = 5.3 T` (reactor size) |
| `circ_chartmap_simple.converter.nc` | Direct EQDSK-to-Boozer-chart-map output on the converter grid; provenance input, not the SIMPLE run input |
| `circ_chartmap_simple.nc` | SIMPLE input: axis-complete uniform grid with 101 radial points, 48 poloidal points, and 4 toroidal planes |
| `wout_circ.nc` | The same field written in VMEC wout format (direct geometry transcription, no equilibrium solve); exercises SIMPLE's native VMEC input path |
| `simple_run/circ_chartmap.nc` | Relative symlink to `circ_chartmap_simple.nc` under the filename used by `simple.in` |

The equilibrium is reactor-size so that both the 5 keV deuterium benchmark
particle and SIMPLE's default 3.5 MeV alpha are confined at mid-radius; on the
earlier small machine (`R0 = 1.6 m`) a 3.5 MeV alpha left through the edge
before one bounce.

The generator is [`gen_circular_eqdsk.py`](gen_circular_eqdsk.py) in this
directory. It derives from NEO-RT's public POTATO gate generator but writes
reactor-scale parameters and a poloidal flux that realizes the flux-surface
safety factor `q(r) = 1.5 + 2.5 (r/a)^2` exactly on the concentric-circle
model (the NEO-RT original uses a cylindrical midplane estimate), and it
continues `psi` smoothly beyond the LCFS instead of clamping it at the
boundary value (the clamp puts a gradient kink at the LCFS that spline
readers smear over one or two grid cells, corrupting flux-surface
quantities in the outermost few percent of the plasma). The file is
therefore no longer byte-identical to NEO-RT's golden-record input.
Regenerating the EQDSK needs libneo with the fixed-column g-file header
writer (itpplasma/libneo#396); regenerating the chart map needs libneo with
the corrected covariant `B_theta` and the `psimax` boundary stop in the
EQDSK-to-chartmap converter (itpplasma/libneo#399).

## Converting tokamak equilibria (any g-file)

Two repository scripts turn an axisymmetric EQDSK g-file into SIMPLE inputs;
both work for shaped tokamaks, not only this circular case:

```sh
# EQDSK -> Boozer chart map on the SIMPLE grid (libneo converter + regrid)
PYTHONPATH=../libneo/python:../SIMPLE/build/_deps/libneo-build \
    python tools/eqdsk_to_simple_chartmap.py input.eqdsk output_chartmap.nc

# EQDSK -> VMEC wout NetCDF (direct transcription of the g-file field)
PYTHONPATH=../libneo/python:../SIMPLE/build/_deps/libneo-build \
    python tools/eqdsk_to_vmec.py input.eqdsk wout_output.nc
```

`eqdsk_to_simple_chartmap.py` reproduces `circ_chartmap_simple.nc` (and keeps
the converter-grid intermediate next to it). It stops the converter's
flux-surface scan at the g-file boundary flux (`psimax`). This matters for
synthetic equilibria without an X-point, such as this circular case: their
`psi` keeps rising beyond the LCFS, so the converter's default field-line
"separatrix search" runs to the computational box edge and would place `s = 1`
on the wrong surface (besides taking hours instead of seconds).

`eqdsk_to_vmec.py` reproduces `wout_circ.nc`. Its default mode transcribes
the g-file's own flux-surface geometry into wout Fourier variables using the
straight-field-line angle (`lambda = 0`), with `iota` from the flux-surface
line integral `q = f/(2 pi) * loop-integral dl/(R |grad psi|)`. No
equilibrium is solved, so SIMPLE's VMEC input path traces the SAME field as
the g-file and POTATO. With `--resolve`, VMEC++ instead re-solves the
fixed-boundary equilibrium from the extracted boundary, pressure, and iota
profile; a model g-file that is not a Grad-Shafranov solution then acquires
its Shafranov shift (~12 cm axis shift here) and orbit quantities move at
the percent level (trapped precession by ~13 %). Use `--resolve` only when
a force-balanced variant is wanted explicitly.

Committed equilibrium checksums:

```text
1f1a0e9876ee79cdeb26ef9b96b8944f9811222f9b38793dfade33ab77fb1b25  potato_run/circ.eqdsk
5b7f96a04bb0b42e0585fe8beabb3b3925115354ad5bb9aba6402880cdfe7189  circ_chartmap_simple.converter.nc
0af60ff1b301de513f4681658e7b03c3a786621beb20585e50caab4a49d5dd0e  circ_chartmap_simple.nc
13fdc38cf1847d3e56df550ab3b9f2e55a6910e49719e75262e86badbac9d630  wout_circ.nc
```

Three-way validation of the committed files (5 keV deuteron at
`s_tor = 0.3`; SIMPLE canonical frequency API on both inputs vs POTATO
`itest_type = 4` on the same EQDSK, all in rad/s):

| quantity | SIMPLE chart map | SIMPLE wout | POTATO |
|---|---|---|---|
| trapped `xi = 0.3` `omega_b` | 1.2809e4 | 1.2877e4 | 1.288e4 |
| trapped `xi = 0.3` `omega_phi` | -273.3 | -268.4 | -272.6 |
| passing `xi = 0.8` `omega_b` | 3.7358e4 | 3.7371e4 | 3.738e4 |
| passing `xi = 0.8` `omega_phi` | 8.4556e4 | 8.5058e4 | 8.506e4 |

All pairwise differences are below 2 % (most below 0.5 %); the remaining
spread reflects interpolation and the different radial/angle
discretizations, not physics differences. The alpha seeds (trapped and
passing) also return `status = SUCCESS` on both SIMPLE inputs.

All chart-map geometry uses centimetres. The radial coordinate is
`rho_tor = sqrt(s_tor)`. The VMEC file uses VMEC conventions (SI, `s_tor`).

No experimental equilibrium is stored here. Rung 5 remains blocked until a
shareable equilibrium and its provenance are agreed.

## Benchmark particle

The Rung 0 particle is a 5 keV deuteron: `n_e = 1`, `n_d = 2`,
`facE_al = 700` in `simple.in`, `E_alpha = 5d3` eV with `A = 2`, `Z = 1` in
`potato.in`. The same overrides apply to the Python API:

```python
pysimple.init("rung0/circ_chartmap_simple.nc", deterministic=True,
              n_e=1, n_d=2, facE_al=700.0, ...)
```

SIMPLE's default 3.5 MeV alpha is also confined on this equilibrium, so quick
API smoke tests without overrides now work; quantitative comparisons against
POTATO still require the deuteron settings above.

Mind the radial coordinate of the seed: the chart map's public coordinate is
`rho_tor = sqrt(s_tor)` (benchmark surface: `sqrt(0.3) = 0.5477...`), while
the VMEC file uses `s_tor` directly (same surface: `0.3`).

## Prepare writable run directories

Run these commands from the repository root after building the sibling
checkouts described in `WORKPLAN.md`. Keep the same shell for the subsequent
SIMPLE and POTATO commands because they use `run_root`.

```sh
run_root=$(mktemp -d "${TMPDIR:-/tmp}/simple-potato-rung0.XXXXXX")

mkdir -p "$run_root/simple" "$run_root/potato"
cp rung0/simple_run/simple.in rung0/simple_run/start.dat "$run_root/simple/"
cp rung0/circ_chartmap_simple.nc "$run_root/simple/circ_chartmap.nc"
cp rung0/potato_run/circ.eqdsk \
   rung0/potato_run/convexwall.dat \
   rung0/potato_run/field_divB0.inp \
   rung0/potato_run/potato.in \
   rung0/potato_run/profile_poly.in \
   "$run_root/potato/"

printf 'Run directory: %s\n' "$run_root"
```

Copying the files explicitly avoids breaking the relative chart-map symlink
in `simple_run/`. `eqmagprofs.dat` is generated by POTATO and is not an input.

## Run SIMPLE

```sh
SIMPLE_EXE=$(realpath ../SIMPLE/build/simple.x)
(cd "$run_root/simple" && "$SIMPLE_EXE")
```

The supplied trapped case has

```text
s_tor = 0.3
rho_tor = sqrt(0.3) = 0.5477225575051661
xi = v_parallel / v = 0.3
```

For chart-map `startmode = 2`, column 1 of `start.dat` is `rho_tor`, not
`s_tor`; column 5 is `xi`. A passing check with `xi = 0.8` also requires
`contr_pp = -1d10` in `simple.in`, because the default `contr_pp = -1` skips
deeply passing particles and leaves their trajectory as NaN.

Successful output includes `orbits.nc`, `times_lost.dat`, and
`confined_fraction.dat`. With SIMPLE commit `08b85a1` or newer,
`orbits.nc` contains direct cylindrical `R` and `Z` variables. Check their
`units` attributes before conversion; they are `cm` for this chart map. The
file also records

```text
coordinate_type = "chartmap"
radial_coordinate = "rho"
```

The legacy variable named `s` therefore contains `rho_tor` in this run.
Untraced samples remain NaN.

Check the output contract:

```sh
python - "$run_root/simple/orbits.nc" <<'PY'
import sys

import netCDF4
import numpy as np

with netCDF4.Dataset(sys.argv[1]) as data:
    assert data.coordinate_type == "chartmap"
    assert data.radial_coordinate == "rho"
    assert data.variables["R"].units == "cm"
    assert data.variables["Z"].units == "cm"
    radius = np.asarray(data.variables["R"][:])
    height = np.asarray(data.variables["Z"][:])

finite = np.isfinite(radius) & np.isfinite(height)
assert finite.any()
print(f"finite cylindrical samples: {finite.sum()}")
PY
```

## Run POTATO

```sh
POTATO_EXE=$(realpath ../NEO-RT/POTATO/build/potato.x)
(cd "$run_root/potato" && "$POTATO_EXE")
```

With `itest_type = 4`, the trajectory is written to `fort.100`. Its columns
are `R phi Z p xi dxi/dtau`, although the POTATO input and output interfaces
retain the historical name `lambda` for `xi`. `R` and `Z` are in centimetres.

The last row of `profile_poly.in` is deliberately zero: Rung 0 has no radial
electric potential, matching the SIMPLE Hamiltonian. A nonzero last row can
change the orbit class and makes POTATO's normalized momentum vary.
`gen_circular_eqdsk.py` now writes the zero row directly.

## Post-process the trajectories

The benchmark-side analysis tools grew out of Majid Khan's first quantitative
comparison. Run them from the repository root after staging both outputs:

```sh
python tools/analyze_simple_orbit.py "$run_root/simple/orbits.nc" \
  --json "$run_root/simple/orbit-diagnostics.json"
python tools/compare_orbits.py \
  "$run_root/simple/orbits.nc" "$run_root/potato/fort.100" \
  --output "$run_root/orbit-overlay.png"
```

For a POTATO `itest_type = 5` run:

```sh
python tools/plot_potato_frequency_scan.py \
  "$run_root/potato/freq_scan.dat" \
  --output "$run_root/potato-frequency-scan.png"
```

The SIMPLE analyzer follows the benchmark definition in
`doc/benchmark.tex`: it interpolates the `v_par = 0` tips and measures a full
bounce from tip `k` to tip `k+2`. POTATO's `find_bounce` uses a refined return
section, but both estimators represent one physical bounce. A discrepancy is
therefore a comparison diagnostic, not a reason to silently change the SIMPLE
event definition.

In `freq_scan.dat`, `omega_b` and `omega_phi` are already physical angular
frequencies. The `taub` column is POTATO's length-like normalized bounce
quantity; physical time is `taub/v0`, not `taub` itself. Match the first finite
SIMPLE `R,Z` point to POTATO before comparing rows—matching a radial minimum or
assigning equal numerical values to `rho_tor` and `rho_pol` is not equivalent.

## Rung 0 acceptance

Rung 0 is complete when:

1. both programs exit successfully in fresh run directories;
2. the SIMPLE trapped orbit has finite `R` and `Z`, and `v_par` changes sign;
3. POTATO writes a finite single-orbit trajectory terminated by its NaN line;
4. the input paths, code commits, and output units are recorded in the result
   note.

The next step is Rung 1 in `WORKPLAN.md`: measure the EQDSK-to-chart-map field
representation error, match one physical initial point, and overlay the two
orbits in `(R,Z)`.
