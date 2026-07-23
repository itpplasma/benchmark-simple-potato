# ==========================================================
# Resonance Map Generator
#
# This program scans invariant space (P_phi, J_perp) at
# fixed energy H0=1, reconstructs particle states using
# states_from_invariants(), computes the canonical orbit
# frequencies, and produces resonance maps.
# ==========================================================
from pathlib import Path

import numpy as np
import pysimple

# ==========================================================
# Workflow
#
# 1. Initialize SIMPLE
# 2. Define equilibrium and particle energy
# 3. Specify an invariant point (P_phi, J_perp)
# 4. Reconstruct particle state
# 5. Compute canonical frequencies
# 6. Extend to a grid in (P_phi, J_perp)
# 7. Produce resonance maps
# ==========================================================

print("Using pysimple from:")
print(pysimple.__file__)
print()

print()
print("="*72)
print("                 RESONANCE MAP GENERATOR")
print("="*72)
print()
print()
# ----------------------------------------------------------------------
# Reference particle
#----------------------------------------------------------------------
# Used only to define the reference energy H0.
# Future versions will work directly in (P_phi, J_perp) space.
# ----------------------------------------------------------------------

particle_kwargs = dict(n_e=1, n_d=2, facE_al=700.0)
#xi = 0.8        # passing
xi = 0.3      # trapped


# ==========================================================
#  equilibrium: "chartmap" (circular benchmark) or "vmec"
# ==========================================================
CASE ="chartmap"     # "chartmap" or "vmec"

REPO = Path(__file__).resolve().parent

if CASE =="chartmap":
    equilibrium = str(REPO / "rung0" / "aug_30835_chartmap.nc")
    # Benchmark particle:
    # 5 keV deuteron (matches Rung 0 benchmark).
    # Without these overrides SIMPLE uses its default
    # 3.5 MeV alpha particle.
    # Benchmark seed: rho_tor = sqrt(0.3), theta = 0, phi = 0, v/v0 = 1.
    particle = np.array([0.5520740743749941,1.182620756878805e-08,0.0,1.0,xi])
    #particle = np.array([np.sqrt(0.3), 0.0, 0.0, 1.0, xi])

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
# Compute invariants of the initial particle
# ==========================================================

print("="*72)
print("INITIAL PARTICLE INVARIANTS")
print("="*72)

inv = pysimple.invariants_from_state(particle)
# ==========================================================
# Reference invariants
# ==========================================================

H0 = inv["h0"]
Jperp = inv["j_perp"]
Pphi = inv["p_phi"]
print(inv)
print()

# ---------------------------------------------------------------
# Sensitivity tests:
# ---------------------------------------------------------------
Pphi_test_factors = [
    0.99,
    0.995,
    1.000,
    1.005,
    1.010,
]

Jperp_test_factors = [
    0.90,
    0.95,
    1.00,
    1.05,
    1.10,
]

print("=" * 72)
print("REFERENCE INVARIANTS")
print("=" * 72)

print(f"H0     = {H0:.12e}")
print(f"Pphi   = {Pphi:.12e}")
print(f"Jperp  = {Jperp:.12e}")

# ==========================================================
# Invariant point to analyse
# ==========================================================

#H0_scan = H0
#Jperp_scan = Jperp
#Pphi_scan = Pphi

# ==========================================================
# Scan offsets in invariant space
# ==========================================================
#Pphi_offsets = np.linspace(-5e-3, 5e-3, 21)
Jperp_offsets = np.linspace(-2e-5, 2e-5, 11)
Pphi_offsets = np.linspace(-0.01 * Pphi, 0.01 * Pphi, 21)

print()
print("=" * 72)
print("RELATIVE SCAN WIDTH")
print("=" * 72)

print(f"Pphi range  = {Pphi_offsets.min():+.3e} ... {Pphi_offsets.max():+.3e}")
print(f"Jperp range = {Jperp_offsets.min():+.3e} ... {Jperp_offsets.max():+.3e}")

print()
print(f"ΔPphi/Pphi   = {Pphi_offsets.max()/abs(Pphi):.3e}")
print(f"ΔJperp/Jperp = {Jperp_offsets.max()/abs(Jperp):.3e}")
# ==========================================================
# Storage for scan results
# ==========================================================
print()
print("=" * 72)
print("Pphi SENSITIVITY TEST")
print("=" * 72)

for factor in Pphi_test_factors:

    Pphi_test = Pphi * factor

    states = pysimple.states_from_invariants(
        H0,
        Jperp,
        Pphi_test,
    )

    if len(states["states"]) == 0:
        print(f"{factor:6.3f}   No physical orbit")
        continue

    particle_test = states["states"][0]

    result = pysimple.compute_canonical_frequencies(
        particle_test,
        integrator="midpoint",
        n_periods=1,
    )

    print(
        f"{factor:6.3f}   "
        f"Pphi={Pphi_test:.6e}   "
        f"{result['orbit_class']:8s}   "
        f"ωb={result['omega_b']:.6e}   "
        f"ωφ={result['omega_phi']:.6e}"
    )

print("=" * 72)
print()

print()
print("=" * 72)
print("Jperp SENSITIVITY TEST")
print("=" * 72)

for factor in Jperp_test_factors:

    Jperp_test = Jperp * factor

    states = pysimple.states_from_invariants(
        H0,
        Jperp_test,
        Pphi,
    )

    if len(states["states"]) == 0:
        print(f"{factor:6.3f}   No physical orbit")
        continue

    particle_test = states["states"][0]

    result = pysimple.compute_canonical_frequencies(
        particle_test,
        integrator="midpoint",
        n_periods=1,
    )

    print(
        f"{factor:6.3f}   "
        f"Jperp={Jperp_test:.6e}   "
        f"{result['orbit_class']:8s}   "
        f"ωb={result['omega_b']:.6e}   "
        f"ωφ={result['omega_phi']:.6e}"
    )

print("=" * 72)
print()


scan_data = []
for dJperp in Jperp_offsets:
    success_count = 0
    for dPphi in Pphi_offsets:

    # Always start from the reference invariants
        H0_scan = H0
        Jperp_scan = Jperp + dJperp
        Pphi_scan = Pphi + dPphi

        print("="*72)
        print(f"dJperp = {dJperp:.2e}")
        print(f"dPphi = {dPphi:.2e}")

        states = pysimple.states_from_invariants(
            H0_scan,
            Jperp_scan,
            Pphi_scan,
        )

        if len(states["states"]) == 0:
            print("No physical orbit found.")
            print()
            continue

        reconstructed_particle = states["states"][0]


        print("Computing canonical frequencies...")
        print()

        result = pysimple.compute_canonical_frequencies(
            reconstructed_particle,
            integrator="midpoint",
            n_periods=1,
        )
        success_count += 1
        scan_data.append(
            {
                "Pphi": Pphi_scan,
                "Jperp": Jperp_scan,

                "dPphi": dPphi,
                "dJperp": dJperp,

                "omega_b": result["omega_b"],
                "omega_phi": result["omega_phi"],
                "orbit_class": result["orbit_class"],
                "orbit_flag": 1 if result["orbit_class"] == "passing" else 0,
            }
        )

        print("Calculation completed.")
        print()

        print(
            f"dJperp={dJperp:+.2e}   "
            f"dPphi={dPphi:+.2e}   "
            f"class={result['orbit_class']:8s}   "
            f"omega_b={result['omega_b']:.6e}   "
            f"omega_phi={result['omega_phi']:.6e}"
        )

        print()
        print("-" * 72)
        print(f"dJperp = {dJperp:+.2e} : {success_count}/{len(Pphi_offsets)} successful")
        print("-" * 72)
        print()
print()
print("="*72)
print("STORED SCAN RESULTS")
print("="*72)

print(f"Number of successful orbit calculations: {len(scan_data)}")
print()

print("First five scan points:")

for point in scan_data[:5]:
    print(point)

import pandas as pd

df = pd.DataFrame(scan_data)

print()
print("=" * 72)
print("DATAFRAME")
print("=" * 72)

with pd.option_context(
    "display.precision", 12,
    "display.max_columns", None,
):
    print(df.head())

df.to_csv("resonance_scan.csv", index=False)

print()
print("Data saved to resonance_scan.csv")


#import matplotlib.pyplot as plt
#plt.figure(figsize=(7,5))

#dPphi_values = [point["dPphi"] for point in scan_data]
#omega_b_values = [point["omega_b"] for point in scan_data]
#omega_phi_values = [point["omega_phi"] for point in scan_data]

#plt.plot(dPphi_values, omega_b_values, "o-", label=r"$\omega_b$")
#plt.plot(dPphi_values, omega_phi_values, "s-", label=r"$\omega_\phi$")
#plt.xlabel(r"$\Delta P_\phi$")
#plt.ylabel("Frequency (rad/s)")
#plt.legend()

#plt.grid(True)

#plt.tight_layout()
#plt.show()
