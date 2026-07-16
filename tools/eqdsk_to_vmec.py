#!/usr/bin/env python3
"""Convert a tokamak EQDSK g-file into a VMEC wout NetCDF.

Default mode is a DIRECT conversion without any equilibrium solve: the
flux surfaces of the g-file's own psi map are traced (libneo
efit_to_boozer), expressed in the straight-field-line angle (so
lambda = 0 exactly), and written as wout Fourier geometry together with
iota(s) and the toroidal flux. SIMPLE reconstructs B purely from that
geometry and the fluxes, so the traced field is the SAME field as the
g-file - orbit results agree with a cylindrical-coordinate reference
(POTATO) at the interpolation level, with no Shafranov systematic.

With --resolve, VMEC++ instead re-solves the fixed-boundary equilibrium
from the extracted boundary, pressure, and iota profile. That yields a
proper force-balanced equilibrium; a model g-file that is not a
Grad-Shafranov solution then acquires, for example, its Shafranov shift,
so orbits differ from the g-file field at the percent level.

Requirements: the libneo python package with its compiled
_efit_to_boozer module (direct mode), or a working vmecpp (--resolve).
With the sibling layout of this repository::

    PYTHONPATH=../libneo/python:../SIMPLE/build/_deps/libneo-build \
        python tools/eqdsk_to_vmec.py \
        rung0/potato_run/circ.eqdsk rung0/wout_circ.nc
"""

from __future__ import annotations

import argparse
import math
import os
import tempfile

import numpy as np

TWOPI = 2.0*math.pi


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


def toroidal_flux_from_psi_map(eq):
    """Toroidal flux and safety factor profiles from the g-file psi map.

    Returns (psi, phi_tor, q, s_tor) on the g-file's radial psi levels.
    phi_tor(psi_k) is the area integral of B_tor = f(psi)/R over the region
    inside the psi_k level and inside the LCFS. This stays correct even when
    a g-file's q column is inconsistent with its own psi map (as in model
    equilibria); q and iota then follow from dphi/dpsi.
    """
    from matplotlib.path import Path as MplPath
    from scipy.interpolate import RectBivariateSpline

    psi = np.linspace(eq["PsiaxisVs"], eq["PsiedgeVs"], eq["nrgr"])

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

    # q = dphi_tor / dPsi_pol_total with Psi_pol_total = 2*pi*psi.
    q = np.gradient(phi_tor, psi)/(2.0*np.pi)
    # The innermost flux levels enclose only a handful of refined grid cells,
    # so replace them by a smooth extrapolation from the well-resolved region.
    fit = np.polynomial.Polynomial.fit(psi[2:8], q[2:8], 2)
    q[:2] = fit(psi[:2])
    # The outermost levels are cut raggedly by the LCFS polygon, which makes
    # dphi/dpsi noisy there; extrapolate them from the clean region as well.
    # (For diverted equilibria whose edge q steepens toward the separatrix
    # this smooths the last few percent of s; acceptable at the resolution
    # the profile is used here.)
    ntail = max(3, eq["nrgr"]//20)
    window = slice(-(ntail + 8), -ntail)
    fit = np.polynomial.Polynomial.fit(psi[window], q[window], 2)
    q[-ntail:] = fit(psi[-ntail:])

    s_tor = phi_tor/phi_tor[-1]
    return psi, phi_tor, q, s_tor


def convert_direct(args) -> None:
    """Write the g-file field itself as wout geometry (no equilibrium solve).

    Surfaces come from libneo's efit_to_boozer field-line integration,
    sampled in the symmetry-flux (straight-field-line) angle, so the wout
    lambda is identically zero and B follows from geometry, iota, and the
    toroidal flux alone.
    """
    import _efit_to_boozer as _etb
    from libneo.eqdsk_base import read_eqdsk
    from libneo.eqdsk_to_boozer_chartmap import (
        _write_convex_wall_from_lcfs,
        _write_field_divB0_inp,
        _write_inp,
    )
    import netCDF4

    gfile = os.path.abspath(args.gfile)
    eq = read_eqdsk(gfile)

    inp_kwargs = dict(nsurfmax=args.nsurfmax, nlabel=args.nlabel)
    if not args.no_psimax and eq["PsiedgeVs"] > eq["PsiaxisVs"]:
        inp_kwargs["psimax"] = float(eq["PsiedgeVs"])*1.0e8
        print(f"flux-surface scan stops at psimax={inp_kwargs['psimax']:.6e} "
              "G cm^2/rad")

    from scipy.interpolate import RectBivariateSpline

    ns = args.ns
    mpol = args.mpol
    ntor = args.ntor
    ntheta = 256
    npsi = 4*ns
    s_grid = np.linspace(0.0, 1.0, ns)
    theta = np.linspace(0.0, TWOPI, ntheta, endpoint=False)

    # Fine poloidal-flux grid (SI Wb/rad relative to the axis) on which the
    # surface geometry is sampled and q is computed.
    psi_lcfs = float(eq["PsiedgeVs"] - eq["PsiaxisVs"])
    psi_levels = np.linspace(0.0, psi_lcfs, npsi + 1)[1:]

    psi_spl = RectBivariateSpline(np.asarray(eq["Z"]), np.asarray(eq["R"]),
                                  np.asarray(eq["PsiVs"]))
    psi_axis_prof = np.linspace(eq["PsiaxisVs"], eq["PsiedgeVs"], eq["nrgr"])
    f_of_psi = np.asarray(eq["fprof"], dtype=float)

    R_pt = np.zeros((npsi, ntheta))
    Z_pt = np.zeros((npsi, ntheta))

    orig_dir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        try:
            _write_inp("efit_to_boozer.inp", gfile, **inp_kwargs)
            lcfs = eq["Lcfs"]
            _write_convex_wall_from_lcfs("convexwall.dat",
                                         lcfs[:, 0], lcfs[:, 1])
            _write_field_divB0_inp("field_divB0.inp", gfile,
                                   convexfile="convexwall.dat")
            _etb.efit_to_boozer.init()
            for k in range(npsi):
                psi_cgs = psi_levels[k]*1.0e8
                for j, th in enumerate(theta):
                    s_f = np.float64(0.0)
                    psi_f = np.float64(psi_cgs)
                    th_f = np.float64(th)
                    out = _etb.efit_to_boozer.magdata(2, s_f, psi_f, th_f)
                    R_pt[k, j] = out[7]/100.0  # cm -> m
                    Z_pt[k, j] = out[10]/100.0
            _etb.efit_to_boozer.deinit()
        finally:
            os.chdir(orig_dir)

    # Safety factor per surface from the flux-surface line integral
    # q = f/(2 pi) * loop integral of dl/(R |grad psi|), robust against the
    # noise of finite-difference or level-set-area constructions.
    q_psi = np.zeros(npsi)
    for k in range(npsi):
        R_s, Z_s = R_pt[k], Z_pt[k]
        dR = (np.roll(R_s, -1) - np.roll(R_s, 1))/2.0
        dZ = (np.roll(Z_s, -1) - np.roll(Z_s, 1))/2.0
        dl = np.hypot(dR, dZ)
        grad_psi = np.hypot(psi_spl(Z_s, R_s, dx=1, grid=False),
                            psi_spl(Z_s, R_s, dy=1, grid=False))
        f_here = np.interp(psi_levels[k] + eq["PsiaxisVs"],
                           psi_axis_prof, f_of_psi)
        q_psi[k] = abs(f_here)/TWOPI*np.sum(dl/(R_s*grad_psi))

    # Toroidal flux from the same q: dPhi = 2 pi q dpsi. This makes the s
    # labeling exactly consistent with iota.
    from scipy.integrate import cumulative_trapezoid
    phi_of_psi = TWOPI*cumulative_trapezoid(q_psi, psi_levels, initial=0.0)
    # add the small axis segment (psi from 0 to the first level)
    phi_of_psi += TWOPI*q_psi[0]*psi_levels[0]
    phiedge = float(phi_of_psi[-1])
    s_of_psi = phi_of_psi/phiedge

    # Fourier decomposition per fine surface (stellarator-symmetric in the
    # m-decomposition; asymmetric shaping would need rmns/zmnc too).
    m_arr = np.arange(mpol)
    cos_mt = np.cos(np.outer(m_arr, theta))
    sin_mt = np.sin(np.outer(m_arr, theta))
    rmnc_p = R_pt @ cos_mt.T*(2.0/ntheta)
    rmnc_p[:, 0] /= 2.0
    zmns_p = Z_pt @ sin_mt.T*(2.0/ntheta)
    zmns_p[:, 0] = 0.0

    # Orient theta counterclockwise (Z ~ +sin(theta)) like VMEC tokamak
    # convention; the field-line integration may deliver either direction.
    if zmns_p[-1, 1] < 0.0:
        zmns_p = -zmns_p
        rmns_check = R_pt @ sin_mt.T*(2.0/ntheta)
        if np.abs(rmns_check[-1, 1:]).max() > 1e-6*np.abs(rmnc_p[-1, 1]):
            raise SystemExit("theta flip requires up-down symmetry; "
                             "asymmetric case not supported in direct mode")

    # Interpolate harmonics and iota from the fine psi grid onto the uniform
    # s grid of the output file.
    rmnc_m = np.zeros((ns, mpol))
    zmns_m = np.zeros((ns, mpol))
    for m in range(mpol):
        rmnc_m[:, m] = np.interp(s_grid, s_of_psi, rmnc_p[:, m])
        zmns_m[:, m] = np.interp(s_grid, s_of_psi, zmns_p[:, m])
    iotaf = 1.0/np.interp(s_grid, s_of_psi, q_psi)

    # Axis row: m=0 from the traced innermost surfaces, m>0 vanish.
    from numpy.polynomial import Polynomial
    axis_fit = Polynomial.fit(s_of_psi[:6], rmnc_p[:6, 0], 2)
    rmnc_m[0, :] = 0.0
    rmnc_m[0, 0] = axis_fit(0.0)
    zmns_m[0, :] = 0.0
    iota_fit = Polynomial.fit(s_of_psi[:6], 1.0/q_psi[:6], 2)
    iotaf[0] = iota_fit(0.0)

    phi = s_grid*phiedge

    # Assemble the (m, n) tables in the standard wout layout; all n /= 0
    # amplitudes are zero (axisymmetric), but SIMPLE's reader needs
    # ntor >= 1 to spline over phi.
    xm, xn = [], []
    for m in range(mpol):
        for n in range(-ntor, ntor + 1):
            if m == 0 and n < 0:
                continue
            xm.append(m)
            xn.append(n)
    xm = np.array(xm, dtype=float)
    xn = np.array(xn, dtype=float)
    mnmax = len(xm)

    rmnc = np.zeros((ns, mnmax))
    zmns = np.zeros((ns, mnmax))
    for i, (m, n) in enumerate(zip(xm.astype(int), xn.astype(int))):
        if n == 0:
            rmnc[:, i] = rmnc_m[:, m]
            zmns[:, i] = zmns_m[:, m]

    with netCDF4.Dataset(args.wout, "w") as d:
        d.createDimension("radius", ns)
        d.createDimension("mn_mode", mnmax)
        d.createDimension("n_tor", ntor + 1)

        def scalar(name, value, dtype="i4"):
            v = d.createVariable(name, dtype)
            v[...] = value

        scalar("lasym__logical__", 0)
        scalar("nfp", 1)
        scalar("ns", ns)
        scalar("mpol", mpol)
        scalar("ntor", ntor)
        scalar("mnmax", mnmax)
        scalar("signgs", -1)
        scalar("Rmajor_p", rmnc_m[0, 0], "f8")
        for name, data, dims in (
            ("phi", phi, ("radius",)),
            ("iotaf", iotaf, ("radius",)),
            ("iotas", iotaf, ("radius",)),
            ("xm", xm, ("mn_mode",)),
            ("xn", xn, ("mn_mode",)),
            ("rmnc", rmnc, ("radius", "mn_mode")),
            ("zmns", zmns, ("radius", "mn_mode")),
            ("lmns", np.zeros((ns, mnmax)), ("radius", "mn_mode")),
            ("raxis_cc", np.array([rmnc_m[0, 0]] + [0.0]*ntor), ("n_tor",)),
            ("zaxis_cs", np.zeros(ntor + 1), ("n_tor",)),
        ):
            v = d.createVariable(name, "f8", dims)
            v[:] = data
        d.eqdsk2wout_source = gfile
        d.conversion = "direct geometry transcription, no equilibrium solve"

    a_edge = rmnc_m[-1, 1]
    print(f"wrote {args.wout}: direct transcription, R_axis={rmnc_m[0,0]:.4f} m, "
          f"a={a_edge:.4f} m, phiedge={phiedge:.4f} Wb, "
          f"iota axis/edge={iotaf[0]:.4f}/{iotaf[-1]:.4f}")


def convert_resolve(args) -> None:
    """Re-solve the fixed-boundary equilibrium with VMEC++ (force balance)."""
    from libneo.eqdsk_base import read_eqdsk
    import vmecpp

    eq = read_eqdsk(args.gfile)

    pressure = np.asarray(eq["ptotprof"], dtype=float)

    psi, phi_tor, q, s_tor = toroidal_flux_from_psi_map(eq)
    phiedge = phi_tor[-1]
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


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("gfile", help="EQDSK g-file input")
    parser.add_argument("wout", help="VMEC wout NetCDF output")
    parser.add_argument("--resolve", action="store_true",
                        help="re-solve the fixed-boundary equilibrium with "
                             "VMEC++ instead of transcribing the g-file "
                             "field directly (adds force balance and its "
                             "Shafranov shift)")
    parser.add_argument("--mpol", type=int, default=8,
                        help="poloidal Fourier resolution")
    parser.add_argument("--ntor", type=int, default=1,
                        help="toroidal Fourier resolution; the field stays "
                             "axisymmetric (all n /= 0 amplitudes are zero), "
                             "but SIMPLE's VMEC reader needs ntor >= 1 to "
                             "spline over phi")
    parser.add_argument("--ns", type=int, default=101,
                        help="radial surfaces (finest multigrid stage with "
                             "--resolve)")
    parser.add_argument("--asym-tol", type=float, default=1e-8,
                        help="relative threshold for enabling lasym "
                             "(--resolve only)")
    parser.add_argument("--nsurfmax", type=int, default=400,
                        help="direct mode: starting points for the "
                             "separatrix search in efit_to_boozer")
    parser.add_argument("--nlabel", type=int, default=200,
                        help="direct mode: internal radial grid of the "
                             "field-line integration")
    parser.add_argument("--no-psimax", action="store_true",
                        help="direct mode: do not stop the flux-surface "
                             "scan at the g-file boundary flux")
    args = parser.parse_args(argv)

    if args.resolve:
        convert_resolve(args)
    else:
        convert_direct(args)


if __name__ == "__main__":
    main()
