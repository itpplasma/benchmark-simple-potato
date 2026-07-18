# ==========================================================
# Test the new SIMPLE canonical frequency API
# ==========================================================

from pathlib import Path

import numpy as np
import pysimple

print("Using pysimple from:")
print(pysimple.__file__)
print()

print()
print("="*60)
print("      SIMPLE CANONICAL FREQUENCY TEST")
print("="*60)
print()
# ----------------------------------------------------------------------
# Benchmark particle (5 keV deuteron)
# ----------------------------------------------------------------------
particle_kwargs = dict(n_e=1, n_d=2, facE_al=700.0)
xi = 0.8        # passing
# xi = 0.3      # trapped


# ==========================================================
#  equilibrium: "chartmap" (circular benchmark) or "vmec"
# ==========================================================
CASE ="vmec"     # "chartmap" or "vmec"

REPO = Path(__file__).resolve().parent

if CASE =="chartmap":
    equilibrium = str(REPO / "rung0" / "circ_chartmap_simple.nc")
    # Benchmark particle:
    # 5 keV deuteron (matches Rung 0 benchmark).
    # Without these overrides SIMPLE uses its default
    # 3.5 MeV alpha particle.
    # Benchmark seed: rho_tor = sqrt(0.3), theta = 0, phi = 0, v/v0 = 1.
    particle = np.array([np.sqrt(0.3), 0.0, 0.0, 1.0, xi])

else:
    equilibrium = str(REPO / "rung0" / "wout_circ.nc")
    print(f"Loading equilibrium: {equilibrium}")
    # particle_kwargs = {}   # default 3.5 MeV alpha
    particle = np.array([0.3, 0.0, 0.0, 1.0, xi])  

# ==========================================================
# Initialize SIMPLE
# ==========================================================

print("Initializing SIMPLE...")
print()
npoiper2 = 1024   # 2048
trace_time = 1e-3
print("Before init")
print(f"Equilibrium file : {equilibrium}")
pysimple.init(
    equilibrium,
    deterministic=True,
    ntestpart=1,
    npoiper2=npoiper2,
    trace_time=trace_time,
    **particle_kwargs,
)

print("Initialization completed.")
print()

# ==========================================================
# Compute canonical frequencies
# ==========================================================

print("Computing canonical frequencies...")
print()

result = pysimple.compute_canonical_frequencies(
    particle,
    integrator="midpoint",
    n_periods=1,
)

print("Calculation completed.")
print()

print("="*72)
print("Initial Particle")
print("="*72)

print(f"s / rho_tor                : {particle[0]:.6f}")
print(f"theta (rad)                : {particle[1]:.6f}")
print(f"phi (rad)                  : {particle[2]:.6f}")
print(f"v/v0                       : {particle[3]:.6f}")
print(f"xi                         : {particle[4]:.6f}")
print()

# ==========================================================
# Results
# ==========================================================

print("="*72)
print("          CANONICAL ORBIT FREQUENCY RESULTS")
print("="*72)

status_string = "SUCCESS" if result["status"] == 0 else f"ERROR ({result['status']})"

equilibrium_name = "VMEC" if CASE == "vmec" else "Chart-map"

print(f"{'Status':30}: {status_string}")
print(f"{'Orbit class':30}: {result['orbit_class']}")
print(f"{'Equilibrium':30}: {equilibrium_name}")
print()

print(f"Bounce period (s)            : {result['period']:.10e}")
print(f"Bounce period std (s)        : {result['period_std']:.10e}")
print()

print(f"Bounce frequency (rad/s)     : {result['omega_b']:.10e}")
print(f"Toroidal frequency (rad/s)   : {result['omega_phi']:.10e}")
print()

print(f"Toroidal displacement (rad)  : {result['delta_phi']:.10e}")
print(f"Displacement std (rad)       : {result['delta_phi_std']:.10e}")
print()

print(f"Parallel direction           : {result['parallel_direction']}")
print(f"Number of periods            : {result['n_periods']}")
print(f"Integration steps            : {result['n_steps']}")

print("="*72)
print("Run Settings")
print("="*72)

print(f"Integrator                 : midpoint")
print(f"n_periods                  : {result['n_periods']}")
print(f"npoiper2                   : {npoiper2}")
print(f"Trace time (s)             : {trace_time}")
print()

# ==========================================================
# FREQUENCY Consistency check (only meaningful on success)
# ==========================================================

if result["status"] == 0:
    omega_check = 2.0 * np.pi / result["period"]
    difference = omega_check - result["omega_b"]

    print()
    print("="*72)
    print("               CONSISTENCY CHECK")
    print("="*72)

    print(f"2*pi/period (rad/s)          : {omega_check:.10e}")
    print(f"Returned omega_b (rad/s)     : {result['omega_b']:.10e}")
    print(f"Difference                   : {difference:.3e}")

    if abs(difference) < 1.0e-8:
        print("\nPASS: omega_b agrees with  2*pi/period")
    else:
        print("\nWARNING: omega_b differs from 2*pi/period")
else:
    print()
    print(f"Non-success status {result['status']}: consistency check skipped.")
