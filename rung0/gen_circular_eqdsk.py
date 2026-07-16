#!/usr/bin/env python3
"""Generate the reactor-size circular tokamak benchmark inputs.

Derived from NEO-RT's public POTATO gate generator
(POTATO/test/golden_record_resonance/gen_circular_eqdsk.py) with two
benchmark-specific changes:

1. Reactor-scale parameters (ITER-like R0 = 6.2 m, a = 2.0 m, B0 = 5.3 T),
   so that both the 5 keV deuterium benchmark particle and SIMPLE's default
   3.5 MeV alpha are confined at mid-radius.
2. A corrected flux-surface safety factor. For concentric circular surfaces
   with psi per radian and B_pol = psi'(r)/R, the field-line q is

       q(r) = f * r / (psi'(r) * sqrt(R0^2 - r^2)),   f = R0*B0,

   so the generator inverts this exactly: psi'(r) = f*r/(q(r)*sqrt(R0^2-r^2)).
   (The NEO-RT original uses the cylindrical midplane estimate, which is fine
   for its resonance-count gate but makes the written q column inconsistent
   with the psi map.)

The equilibrium remains a model field: concentric circular flux surfaces are
not a Grad-Shafranov solution, so the VMEC re-solve of the same boundary and
profiles (tools/eqdsk_to_vmec.py) acquires a Shafranov shift. This is
documented in rung0/README.md.

The benchmark case has no radial electric field: the potential row of
profile_poly.in is written as zeros (see the 2026-07-10 orbit-class defect
documented in the project notes; a nonzero potential makes POTATO's xi=0.3
seed passing).

Outputs (written next to this script):
  potato_run/circ.eqdsk     g-file consumed by POTATO via field_divB0.inp
  potato_run/convexwall.dat circular stretch-coords wall (cm)
  potato_run/profile_poly.in monotonic density/temperature, zero potential
"""
import os
import sys

import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
try:
    from libneo.eqdsk_base import read_eqdsk, write_eqdsk
except ImportError:  # fall back to sibling checkouts
    for candidate in (
        os.path.join(THIS, "..", "..", "libneo", "python"),
        os.path.join(THIS, "..", "..", "SIMPLE", "build", "_deps",
                     "libneo-src", "python"),
    ):
        if os.path.isdir(candidate):
            sys.path.insert(0, candidate)
            break
    from libneo.eqdsk_base import read_eqdsk, write_eqdsk  # noqa: E402

# --- circular equilibrium parameters (public-safe synthetic) -----------------
R0 = 6.20        # magnetic axis major radius [m]
A = 2.00         # minor radius of the boundary flux surface [m]
B0 = 5.3         # toroidal field on axis [T]
P0 = 5.0e3       # on-axis pressure [Pa]; low beta
# Monotonically rising (sheared) safety factor q(r) = Q0 + (QA-Q0)*(r/a)^2.
Q0 = 1.5         # safety factor on axis
QA = 4.0         # safety factor at the boundary r=a
NR = 129
NZ = 129


def q_of_r(r):
    """Sheared flux-surface safety factor as a function of minor radius."""
    return Q0 + (QA - Q0) * (r / A) ** 2


def psi_of_r(r_max=A):
    """Poloidal flux psi(r) (per radian) that realizes q_of_r exactly.

    For concentric circular surfaces the flux-surface q is
    q = f*r/(psi'(r)*sqrt(R0^2-r^2)) with f = R0*B0, so
    psi'(r) = R0*B0*r/(q(r)*sqrt(R0^2-r^2)). Integrate outward, psi(0)=0.

    r_max may exceed the boundary radius A: the same smooth formula then
    continues psi outside the LCFS. This matters for the psi MAP: clamping
    psi at its boundary value outside the LCFS puts a gradient kink at
    r = A that the bicubic spline of any reader smears over one or two
    grid cells, corrupting |grad psi| (and with it every flux-surface
    quantity) in the outermost few percent of the plasma.
    """
    n = 4096
    r = np.linspace(0.0, r_max, n)
    dpsidr = np.zeros(n)
    dpsidr[1:] = (R0 * B0 * r[1:]
                  / (q_of_r(r[1:]) * np.sqrt(R0 ** 2 - r[1:] ** 2)))
    psi = np.concatenate([[0.0], np.cumsum(0.5 * (dpsidr[1:] + dpsidr[:-1])
                                          * np.diff(r))])
    return r, psi


def build_eqdsk():
    # Box generously around the boundary so orbits stay inside the grid.
    rboxleft = R0 - 1.6 * A
    rboxlength = 3.2 * A
    zboxmid = 0.0
    zboxlength = 3.2 * A

    R = rboxleft + np.linspace(0.0, rboxlength, NR)
    Z = zboxmid - 0.5 * zboxlength + np.linspace(0.0, zboxlength, NZ)
    RR, ZZ = np.meshgrid(R, Z)  # shape (NZ, NR), matching read/write layout

    rho = np.sqrt((RR - R0) ** 2 + ZZ ** 2)
    r_tab, psi_tab = psi_of_r(r_max=1.02 * float(rho.max()))
    psi_edge = float(np.interp(A, r_tab, psi_tab))
    PsiVs = np.interp(rho, r_tab, psi_tab)

    # Toroidal field: constant f = R0*B0 (vacuum, low beta), so B_tor = f/R.
    fprof = np.full(NR, R0 * B0)
    fdfdpsiprof = np.zeros(NR)  # f constant -> f f' = 0

    # Small parabolic pressure in normalized flux s=psi/psi_edge, zero at edge.
    s = np.linspace(0.0, 1.0, NR)
    ptotprof = P0 * (1.0 - s)
    dpressdpsiprof = np.full(NR, -P0 / psi_edge)

    # q profile on the uniform s grid (s = psi/psi_edge -> r via the inverse).
    r_on_s = np.interp(s * psi_edge, psi_tab, r_tab)
    qprof = q_of_r(r_on_s)

    # Plasma current from Ampere on the boundary flux surface:
    # closed-loop integral of B_pol dl = psi'(a)*a*2*pi/sqrt(R0^2-a^2) = mu0*Ip.
    mu0 = 4.0e-7 * np.pi
    dpsida = R0 * B0 * A / (QA * np.sqrt(R0 ** 2 - A ** 2))
    Ip = 2.0 * np.pi * A * dpsida / (mu0 * np.sqrt(R0 ** 2 - A ** 2))

    # Circular boundary and a slightly larger circular limiter.
    theta = np.linspace(0.0, 2.0 * np.pi, 129)
    lcfs = np.column_stack([R0 + A * np.cos(theta), A * np.sin(theta)])
    lim = np.column_stack(
        [R0 + 1.1 * A * np.cos(theta), 1.1 * A * np.sin(theta)]
    )

    eqdata = {
        "header": "benchmark circular reactor      CIRC 00000        ",
        "nrgr": NR,
        "nzgr": NZ,
        "rboxlength": rboxlength,
        "zboxlength": zboxlength,
        "R0": R0,
        "rboxleft": rboxleft,
        "zboxmid": zboxmid,
        "Rpsi0": R0,
        "Zpsi0": 0.0,
        "PsiaxisVs": 0.0,
        "PsiedgeVs": psi_edge,
        "Btor_at_R0": B0,
        "Ip": Ip,
        "fprof": fprof,
        "ptotprof": ptotprof,
        "fdfdpsiprof": fdfdpsiprof,
        "dpressdpsiprof": dpressdpsiprof,
        "PsiVs": PsiVs,
        "qprof": qprof,
        "npbound": lcfs.shape[0],
        "nplimiter": lim.shape[0],
        "Lcfs": lcfs,
        "Limiter": lim,
    }
    return eqdata, R, Z


def write_convexwall(path):
    """Circular stretch-coords wall (cm) enclosing the boundary with margin."""
    rw = 1.12 * A * 100.0
    r0_cm = R0 * 100.0
    theta = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False)
    with open(path, "w") as f:
        for t in theta:
            f.write(f"{r0_cm + rw*np.cos(t):24.16e}{rw*np.sin(t):24.16e}\n")


def write_profile_poly(path):
    """Monotonic species profiles; the potential row is zero by definition.

    POTATO reads ten coefficients per array (descending powers of s_pol).
    Rung 0 has no radial electric field, matching the SIMPLE Hamiltonian.
    """
    def lin(slope, intercept):
        return [0.0] * 8 + [slope, intercept]

    rows = [
        lin(-2.5e13, 5.0e13),   # density [cm^-3]
        [0.0] * 10,             # dummy (reader discards)
        lin(-1.4e3, 2.0e3),     # temperature [eV]
        [0.0] * 10,             # potential: zero radial electric field
    ]
    with open(path, "w") as f:
        f.write("% Public-safe synthetic circular-case profiles "
                "(NOT AUG #30835)\n")
        f.write("% ten coefficients per line, descending powers of s_pol; "
                "order: density, dummy, temperature, potential\n")
        for row in rows:
            f.write(" ".join(f"{c:.16e}" for c in row) + "\n")


def main():
    eqdata, R, Z = build_eqdsk()
    out_dir = os.path.join(THIS, "potato_run")
    eqdsk_path = os.path.join(out_dir, "circ.eqdsk")
    write_eqdsk(eqdsk_path, eqdata)

    psi_edge = eqdata["PsiedgeVs"]

    # Round-trip check: read back and confirm the structural fields survive.
    back = read_eqdsk(eqdsk_path)
    assert back["nrgr"] == NR and back["nzgr"] == NZ
    assert abs(back["R0"] - R0) < 1e-6
    assert abs(back["Btor_at_R0"] - B0) < 1e-6
    assert abs(back["PsiedgeVs"] - psi_edge) < 1e-6
    assert np.allclose(back["fprof"], R0 * B0, atol=1e-6)
    assert np.allclose(back["PsiVs"], eqdata["PsiVs"], atol=1e-6 * psi_edge)

    write_convexwall(os.path.join(out_dir, "convexwall.dat"))
    write_profile_poly(os.path.join(out_dir, "profile_poly.in"))

    print("wrote circ.eqdsk, convexwall.dat, profile_poly.in")
    print(f"  R0={R0} a={A} B0={B0} psi_edge={psi_edge:.4f} "
          f"q_axis={eqdata['qprof'][0]:.3f} q_edge={eqdata['qprof'][-1]:.3f} "
          f"Ip={eqdata['Ip']:.3e} A")


if __name__ == "__main__":
    main()
