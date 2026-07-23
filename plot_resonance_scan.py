import pandas as pd

# ---------------------------------------------------------
# Read resonance scan
# ---------------------------------------------------------

df = pd.read_csv("resonance_scan.csv")
df["omega_ratio"] = df["omega_phi"] / df["omega_b"]

df = df.sort_values(["dJperp", "dPphi"]).reset_index(drop=True)
print()
print("=" * 72)
print("DATA")
print("=" * 72)

print(df.head())

print()
print(df.info())
print()
print("=" * 72)
print("ORBIT CLASS COUNTS")
print("=" * 72)

print(df["orbit_class"].value_counts())
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Select one Jperp slice
# ---------------------------------------------------------

dJ = -2.0e-5

slice_df = df[df["dJperp"] == dJ]

print()
print("=" * 72)
print(f"SLICE: dJperp = {dJ:.2e}")
print("=" * 72)

print(f"Number of points: {len(slice_df)}")

print()
print(slice_df.head())
# ==========================================================
# Figure 1
# Single Jperp slice
# ==========================================================
plt.figure(figsize=(8,5))

plt.plot(
    slice_df["dPphi"],
    slice_df["omega_b"],
    "o-",
    linewidth=2,
    markersize=5,
    label=r"$\omega_b$"
)

plt.plot(
    slice_df["dPphi"],
    slice_df["omega_phi"],
    "s-",
    linewidth=2,
    markersize=5,
    label=r"$\omega_\phi$"
)

plt.xlabel(r"$\Delta P_\phi$")
plt.ylabel("Frequency (rad/s)")
plt.title(f"dJperp = {dJ:.2e}")

plt.grid(True)
plt.legend()

plt.tight_layout()
#plt.show()

# ==========================================================
# Figure 2
# All Jperp slices
# ==========================================================
plt.figure(figsize=(8,5))

for dJ in sorted(df["dJperp"].unique()):

    slice_df = df[df["dJperp"] == dJ]

    plt.plot(
        slice_df["dPphi"],
        slice_df["omega_b"],
        marker="o",
        linewidth=1.5,
        markersize=3,
        label=f"{dJ:.1e}"
    )

plt.xlabel(r"$\Delta P_\phi$")
plt.ylabel(r"$\omega_b$ (rad/s)")
plt.title(r"Bounce frequency for all $J_\perp$ slices")

plt.grid(True)

plt.legend(
    title=r"$\Delta J_\perp$",
    fontsize=8,
    ncol=2
)


# ---------------------------------------------------------
# Check scan grid
# ---------------------------------------------------------

print()
print("=" * 72)
print("SCAN GRID")
print("=" * 72)

print("Unique dJperp values:")
print(sorted(df["dJperp"].unique()))

print()

print("Unique dPphi values:")
print(sorted(df["dPphi"].unique()))

# ==========================================================
# Figure 3
# Successful orbit locations
# ==========================================================

plt.figure(figsize=(7,6))

plt.scatter(
    df["dPphi"],
    df["dJperp"],
    c="blue",
    s=40
)

plt.xlabel(r"$\Delta P_\phi$")
plt.ylabel(r"$\Delta J_\perp$")

plt.title("Successful orbit calculations")

plt.grid(True)

plt.tight_layout()
# ==========================================================
# Figure 4
# Bounce frequency map
# ==========================================================
plt.figure(figsize=(8,6))

sc = plt.scatter(
    df["dPphi"],
    df["dJperp"],
    c=df["omega_b"],
    cmap="viridis",
    s=70
)

plt.xlabel(r"$\Delta P_\phi$")
plt.ylabel(r"$\Delta J_\perp$")

plt.title(r"Bounce frequency $\omega_b$")

plt.grid(True)

cbar = plt.colorbar(sc)
cbar.set_label(r"$\omega_b$ (rad/s)")

plt.tight_layout()
# ==========================================================
# Figure 5
# Orbit class map
# ==========================================================

plt.figure(figsize=(8,6))

passing = df[df["orbit_class"] == "passing"]
trapped = df[df["orbit_class"] == "trapped"]

plt.scatter(
    passing["dPphi"],
    passing["dJperp"],
    marker="o",
    s=70,
    label="Passing"
)

plt.scatter(
    trapped["dPphi"],
    trapped["dJperp"],
    marker="s",
    s=70,
    label="Trapped"
)

plt.xlabel(r"$\Delta P_\phi$")
plt.ylabel(r"$\Delta J_\perp$")

plt.title("Orbit classification")

plt.grid(True)
plt.legend()

# ==========================================================
# Figure 6
# Frequency ratio map
# ==========================================================
plt.figure(figsize=(8,6))

sc = plt.scatter(
    df["dPphi"],
    df["dJperp"],
    c=df["omega_ratio"],
    cmap="plasma",
    s=70
)

plt.xlabel(r"$\Delta P_\phi$")
plt.ylabel(r"$\Delta J_\perp$")

plt.title(r"$\omega_\phi/\omega_b$")

plt.grid(True)

cbar = plt.colorbar(sc)
cbar.set_label(r"$\omega_\phi/\omega_b$")

plt.tight_layout()
plt.show()
