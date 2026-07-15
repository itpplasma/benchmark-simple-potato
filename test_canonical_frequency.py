# ==========================================================
# Test the new SIMPLE canonical frequency API
# ==========================================================

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

# ==========================================================
#  equilibrium  #VMEC    or   Circ Chartmap
# ==========================================================

#equilibrium = "/proj/plasma/CODE/khanm/SIMPLE-canonical-frequencies/test/test_data/wout.nc"
equilibrium = "/proj/plasma/CODE/khanm/benchmark-simple-potato/rung0/circ_chartmap_simple.nc"

# ==========================================================
# Initialize SIMPLE
# ==========================================================

print("Initializing SIMPLE...")
print()
npoiper2=1024   #2048
print("Before init")
print(f"Equilibrium file : {equilibrium}") # Require .nc file  
pysimple.init(
    equilibrium,
    deterministic=True,
    ntestpart=1,
    npoiper2=npoiper2, #2048,
    trace_time=1e-3,
)

print("Initialization completed.")
print()
# ==========================================================
# Trapped particle (VMEC test case)
# ==========================================================

particle = np.array([0.4, 0.7, 0.1, 1.0, 0.1])   # Trapped
#particle = np.array([0.4, 0.7, 0.1, 1.0, 0.9])  # Passing
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

print(f"s                          : {particle[0]:.6f}")
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

print(f"Status                       : {status_string}")
print(f"Orbit class                  : {result['orbit_class']}")
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
print(f"Trace time (s)             : 1.0e-3")
print()
# ==========================================================
# FREQUENCY Consistency check
# ==========================================================

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
