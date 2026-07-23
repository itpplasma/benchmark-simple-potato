import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plot_resonance(ax, x, y, omega_b, omega_phi, m, n,
                   color="red", label=None):

    R = m * omega_b - n * omega_phi

    print()
    print(f"{m}:{n} resonance")
    print(f"Minimum = {R.min():.6e}")
    print(f"Maximum = {R.max():.6e}")

    # Only draw if the resonance actually crosses zero
    if R.min() <= 0.0 <= R.max():

        cs = ax.contour(
            x,
            y,
            R,
            levels=[0.0],
            colors=color,
            linewidths=2.5,
        )

        if label is None:
            label = rf"$({n},{m})$"

        ax.clabel(
            cs,
            fmt={0.0: label},
            fontsize=9,
        )

    else:
        print("No resonance in scanned region.")

# ---------------------------------------------------------
# Read resonance scan
# ---------------------------------------------------------

df = pd.read_csv("resonance_scan.csv")
df["omega_ratio"] = df["omega_phi"] / df["omega_b"]

print()
print("=" * 72)
print("ORBIT CLASS COUNTS")
print("=" * 72)

print(df["orbit_class"].value_counts())

df = df.sort_values(["dJperp", "dPphi"]).reset_index(drop=True)
print()
print("=" * 72)
print("DATA")
print("=" * 72)

print(df.head())
print()
print("=" * 72)
print("GRID INFORMATION")
print("=" * 72)

print(f"Total rows : {len(df)}")

Pphi_values = sorted(df["Pphi"].unique())
Jperp_values = sorted(df["Jperp"].unique())

print(f"Unique Pphi values  : {len(Pphi_values)}")
print(f"Unique Jperp values : {len(Jperp_values)}")

print()
print("First five Pphi values:")
print(Pphi_values[:5])

print()
print("All Jperp values:")
for value in Jperp_values:
    print(f"{value:.12e}")

print()
print("=" * 72)
print("BUILDING 2D ARRAYS")
print("=" * 72)

omega_b_grid = (
    df
    .pivot(index="Jperp", columns="Pphi", values="omega_b")
    .sort_index()
)

omega_phi_grid = (
    df
    .pivot(index="Jperp", columns="Pphi", values="omega_phi")
    .sort_index()
)

orbit_flag_grid = (
    df
    .pivot(index="Jperp", columns="Pphi", values="orbit_flag")
    .sort_index()
)

print(f"omega_b grid shape   : {omega_b_grid.shape}")
print(f"omega_phi grid shape : {omega_phi_grid.shape}")

Pphi = omega_b_grid.columns.to_numpy()
# Reference value (middle of the scan)
Pphi0 = np.mean(Pphi)

# Plotting coordinate
dPphi = Pphi - Pphi0

Jperp = omega_b_grid.index.to_numpy()

print()
print(f"Pphi array length   : {len(Pphi)}")
print(f"Jperp array length  : {len(Jperp)}")

omega_b = omega_b_grid.to_numpy()
omega_phi = omega_phi_grid.to_numpy()
orbit_flag = orbit_flag_grid.to_numpy()


# ==========================================================
# Grid spacing
# ==========================================================

dPphi_step = dPphi[1] - dPphi[0]
Jperp_step = Jperp[1] - Jperp[0]

print()
print("=" * 72)
print("GRID SPACING")
print("=" * 72)
print(f"ΔPphi step  = {dPphi_step:.6e}")
print(f"ΔJperp step = {Jperp_step:.6e}")
print()

# ==========================================================
# Frequency gradients
# ==========================================================

domega_b_dJperp, domega_b_dPphi = np.gradient(
    omega_b,
    Jperp_step,
    dPphi_step,
)

domega_phi_dJperp, domega_phi_dPphi = np.gradient(
    omega_phi,
    Jperp_step,
    dPphi_step,
)

# ==========================================================
# Dimensionless sensitivities
# ==========================================================

Swb_Pphi = (dPphi[np.newaxis, :] / omega_b) * domega_b_dPphi
Swb_Jperp = (Jperp[:, np.newaxis] / omega_b) * domega_b_dJperp

Swphi_Pphi = (dPphi[np.newaxis, :] / omega_phi) * domega_phi_dPphi
Swphi_Jperp = (Jperp[:, np.newaxis] / omega_phi) * domega_phi_dJperp

print("=" * 72)
print("FREQUENCY GRADIENTS")
print("=" * 72)

print(f"|dωb/dPphi|   max = {np.max(np.abs(domega_b_dPphi)):.6e}")
print(f"|dωb/dJperp|  max = {np.max(np.abs(domega_b_dJperp)):.6e}")
print()

print(f"|dωφ/dPphi|   max = {np.max(np.abs(domega_phi_dPphi)):.6e}")
print(f"|dωφ/dJperp|  max = {np.max(np.abs(domega_phi_dJperp)):.6e}")
print("=" * 72)
print()

print("=" * 72)
print("DIMENSIONLESS SENSITIVITIES")
print("=" * 72)

print(f"|Sωb(Pphi)|   max = {np.max(np.abs(Swb_Pphi)):.6e}")
print(f"|Sωb(Jperp)|  max = {np.max(np.abs(Swb_Jperp)):.6e}")
print()

print(f"|Sωφ(Pphi)|   max = {np.max(np.abs(Swphi_Pphi)):.6e}")
print(f"|Sωφ(Jperp)|  max = {np.max(np.abs(Swphi_Jperp)):.6e}")
print("=" * 72)
print()

# ==========================================================
# Final resonance map
# ==========================================================

fig, ax = plt.subplots(figsize=(8,6))

filled = ax.contourf(
    dPphi,
    Jperp,
    omega_b,
    levels=20,
)

plt.colorbar(
    filled,
    ax=ax,
    label=r"$\omega_b$ (rad/s)"
)
ax.contour(
    dPphi,
    Jperp,
    orbit_flag,
    levels=[0.5],
    colors="black",
    linestyles="--",
    linewidths=2.5,
)

# ==========================================================
# Resonance family for ASDEX Upgrade
# Fixed toroidal mode number n = 2
# ==========================================================

colors = [
    "red",
    "blue",
    "green",
    "magenta",
    "orange",
    "cyan",
]

n = 2

for m, color in zip(range(6), colors):

    plot_resonance(
        ax,
        dPphi,
        Jperp,
        omega_b,
        omega_phi,
        m,
        n,
        color=color,
    )


ax.set_xlabel(r"$\Delta P_\phi$")
ax.set_ylabel(r"$J_\perp$")
ax.set_title(r"Resonance map in invariant space")

legend_handles = [

    Line2D(
        [0], [0],
        color="black",
        linestyle="--",
        linewidth=2.5,
        label="Trapped–passing boundary",
    ),

    Line2D(
        [0], [0],
        color="red",
        linewidth=2.5,
        label="Resonances (n,m)",
    ),

]

ax.legend(
    handles=legend_handles,
    loc="upper right",
)
ax.grid(True)

plt.tight_layout()
plt.show()
