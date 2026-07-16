#!/usr/bin/env python3
"""Convert a tokamak EQDSK g-file into a VMEC wout NetCDF with VMEC++.

The converter extracts from the g-file:

- the last closed flux surface, refit as a Fourier boundary,
- the pressure profile,
- the safety-factor profile, used as a fixed iota profile (ncurr = 0),
- the toroidal flux, from phi_tor(psi) = 2*pi * integral q(psi) dpsi.

VMEC++ then re-solves the fixed-boundary equilibrium on that data. The
result is a proper force-balanced equilibrium with the same boundary,
pressure, and rotational transform; it is not bit-identical to the g-file
psi map (a model g-file that is not a Grad-Shafranov solution acquires,
for example, its Shafranov shift here).

Profiles are entered as cubic splines on normalized toroidal flux, so any
monotone axisymmetric g-file works, including shaped and up-down
asymmetric ones (detected automatically; sets lasym).

Requirements: the libneo python package (reader) and a working vmecpp.
With the sibling layout of this repository::

    PYTHONPATH=../libneo/python python tools/eqdsk_to_vmec.py \
        rung0/potato_run/circ.eqdsk rung0/wout_circ.nc
"""

from __future__ import annotations

import argparse

import numpy as np


def fourier_fit_boundary(lcfs: np.ndarray, mpol: int):
    """Fit R(theta), Z(theta) Fourier coefficients from closed LCFS points.

    Returns (rbc, rbs, zbc, zbs) arrays of length mpol and the maximum
    relative fit residual. Theta is the geometric poloidal angle about the
    boundary centroid, oriented counterclockwise in the (R, Z) plane to
    match VMEC's angle convention for tokamaks.
    """
    points = np.asarray(lcfs, dtype=float)
    if np.allclose(points[0], points[-1]):
        points = points[:-1]
    r_center = 0.5*(points[:, 0].max() + points[:, 0].min())
    z_center = 0.5*(points[:, 1].max() + points[:, 1].min())
    theta = np.unwrap(np.arctan2(points[:, 1] - z_center,
                                 points[:, 0] - r_center))
    if theta[-1] < theta[0]:
        points = points[::-1]
        theta = theta[::-1]

    theta_uniform = np.linspace(0.0, 2.0*np.pi, 512, endpoint=False)
    order = np.argsort(np.mod(theta, 2.0*np.pi))
    theta_sorted = np.mod(theta, 2.0*np.pi)[order]
    r_sorted = points[order, 0]
    z_sorted = points[order, 1]
    wrap = ([theta_sorted[-1] - 2.0*np.pi], [theta_sorted[0] + 2.0*np.pi])
    theta_ext = np.concatenate([wrap[0], theta_sorted, wrap[1]])
    r_ext = np.concatenate([[r_sorted[-1]], r_sorted, [r_sorted[0]]])
    z_ext = np.concatenate([[z_sorted[-1]], z_sorted, [z_sorted[0]]])
    r_u = np.interp(theta_uniform, theta_ext, r_ext)
    z_u = np.interp(theta_uniform, theta_ext, z_ext)

    rbc = np.zeros(mpol)
    rbs = np.zeros(mpol)
    zbc = np.zeros(mpol)
    zbs = np.zeros(mpol)
    n = theta_uniform.size
    rbc[0] = r_u.mean()
    zbc[0] = z_u.mean()
    for m in range(1, mpol):
        cos_m = np.cos(m*theta_uniform)
        sin_m = np.sin(m*theta_uniform)
        rbc[m] = 2.0*np.dot(r_u, cos_m)/n
        rbs[m] = 2.0*np.dot(r_u, sin_m)/n
        zbc[m] = 2.0*np.dot(z_u, cos_m)/n
        zbs[m] = 2.0*np.dot(z_u, sin_m)/n

    m_grid = np.arange(mpol)[:, None]
    r_fit = (rbc[:, None]*np.cos(m_grid*theta_uniform)
             + rbs[:, None]*np.sin(m_grid*theta_uniform)).sum(axis=0)
    z_fit = (zbc[:, None]*np.cos(m_grid*theta_uniform)
             + zbs[:, None]*np.sin(m_grid*theta_uniform)).sum(axis=0)
    scale = max(r_u.max() - r_u.min(), z_u.max() - z_u.min())
    residual = max(np.abs(r_fit - r_u).max(), np.abs(z_fit - z_u).max())/scale
    return rbc, rbs, zbc, zbs, residual


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("gfile", help="EQDSK g-file input")
    parser.add_argument("wout", help="VMEC wout NetCDF output")
    parser.add_argument("--mpol", type=int, default=8,
                        help="poloidal Fourier resolution")
    parser.add_argument("--ntor", type=int, default=1,
                        help="toroidal Fourier resolution; the field stays "
                             "axisymmetric (all n /= 0 amplitudes are zero), "
                             "but SIMPLE's VMEC reader needs ntor >= 1 to "
                             "spline over phi")
    parser.add_argument("--ns", type=int, default=101,
                        help="radial surfaces on the finest multigrid stage")
    parser.add_argument("--asym-tol", type=float, default=1e-8,
                        help="relative threshold for enabling lasym")
    args = parser.parse_args(argv)

    from libneo.eqdsk_base import read_eqdsk
    import vmecpp

    eq = read_eqdsk(args.gfile)

    psi = np.linspace(eq["PsiaxisVs"], eq["PsiedgeVs"], eq["nrgr"])
    pressure = np.asarray(eq["ptotprof"], dtype=float)

    # Toroidal flux from the psi map itself: phi_tor(psi_k) is the area
    # integral of B_tor = f(psi)/R over the region inside the psi_k level and
    # inside the LCFS. This stays correct even when a g-file's q column is
    # inconsistent with its own psi map (as in model equilibria); q and iota
    # then follow from dphi/dpsi.
    from matplotlib.path import Path as MplPath
    from scipy.interpolate import RectBivariateSpline

    refine = 32
    r_fine = np.linspace(eq["R"][0], eq["R"][-1], refine*eq["nrgr"])
    z_fine = np.linspace(eq["Z"][0], eq["Z"][-1], refine*eq["nzgr"])
    psi_spline = RectBivariateSpline(eq["Z"], eq["R"], eq["PsiVs"])
    psi_fine = psi_spline(z_fine, r_fine)
    rr, zz = np.meshgrid(r_fine, z_fine)

    lcfs_path = MplPath(np.asarray(eq["Lcfs"]))
    inside_lcfs = lcfs_path.contains_points(
        np.column_stack([rr.ravel(), zz.ravel()]), radius=1e-9
    ).reshape(rr.shape)

    f_of_psi = np.asarray(eq["fprof"], dtype=float)
    f_fine = np.interp(psi_fine, psi, f_of_psi,
                       left=f_of_psi[0], right=f_of_psi[-1])
    btor_over_cell = f_fine/rr*(r_fine[1] - r_fine[0])*(z_fine[1] - z_fine[0])

    # Cumulative toroidal flux over cells sorted by psi gives phi_tor at every
    # level in one pass.
    sign = np.sign(eq["PsiedgeVs"] - eq["PsiaxisVs"])
    psi_cells = (sign*psi_fine[inside_lcfs]).ravel()
    flux_cells = btor_over_cell[inside_lcfs].ravel()
    order = np.argsort(psi_cells)
    psi_sorted = psi_cells[order]
    flux_cumulative = np.cumsum(flux_cells[order])
    phi_tor = np.interp(sign*psi, psi_sorted, flux_cumulative,
                        left=0.0, right=flux_cumulative[-1])
    phi_tor -= phi_tor[0]
    phiedge = phi_tor[-1]
    s_tor = phi_tor/phiedge

    # q = dphi_tor / dPsi_pol_total with Psi_pol_total = 2*pi*psi.
    q = np.gradient(phi_tor, psi)/(2.0*np.pi)
    # The innermost flux levels enclose only a handful of refined grid cells,
    # so replace them by a smooth extrapolation from the well-resolved region.
    fit = np.polynomial.Polynomial.fit(psi[2:8], q[2:8], 2)
    q[:2] = fit(psi[:2])
    q_file = np.asarray(eq["qprof"], dtype=float)
    print(f"q from psi map: axis={q[0]:.3f} edge={q[-1]:.3f} "
          f"(g-file column: {q_file[0]:.3f}..{q_file[-1]:.3f})")

    rbc, rbs, zbc, zbs, residual = fourier_fit_boundary(eq["Lcfs"], args.mpol)
    scale = max(abs(rbc).max(), abs(zbs).max())
    lasym = bool(max(abs(rbs[1:]).max(), abs(zbc[1:]).max(),
                     abs(zbc[0])) > args.asym_tol*scale)
    print(f"boundary fit: mpol={args.mpol} residual={residual:.2e} "
          f"lasym={lasym}")

    ntor = args.ntor
    ncol = 2*ntor + 1

    def sparse(coeffs: np.ndarray) -> np.ndarray:
        full = np.zeros((args.mpol, ncol))
        full[:, ntor] = coeffs
        return full

    def axis(value: float) -> np.ndarray:
        full = np.zeros(ntor + 1)
        full[0] = value
        return full

    inp = vmecpp.VmecInput(
        lasym=lasym,
        nfp=1,
        mpol=args.mpol,
        ntor=ntor,
        ns_array=np.array([16, 51, args.ns]),
        ftol_array=np.array([1e-8, 1e-10, 1e-12]),
        niter_array=np.array([2000, 4000, 8000]),
        phiedge=float(phiedge),
        ncurr=0,
        raxis_c=axis(eq["Rpsi0"]),
        zaxis_s=axis(0.0),
        raxis_s=axis(0.0) if lasym else None,
        zaxis_c=axis(eq["Zpsi0"]) if lasym else None,
        rbc=sparse(rbc),
        zbs=sparse(zbs),
        rbs=sparse(rbs) if lasym else None,
        zbc=sparse(zbc) if lasym else None,
    )

    iota = 1.0/q
    inp = vmecpp.set_profile(inp, "pressure",
                             lambda s: np.interp(s, s_tor, pressure))
    inp = vmecpp.set_profile(inp, "iota",
                             lambda s: np.interp(s, s_tor, iota))

    output = vmecpp.run(inp)
    output.wout.save(args.wout)
    wout = output.wout
    print(f"wrote {args.wout}: volume={wout.volume_p:.4e} m^3, "
          f"b0={wout.b0:.4f} T, aspect={wout.aspect:.3f}")


if __name__ == "__main__":
    main()
