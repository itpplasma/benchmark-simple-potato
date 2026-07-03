"""Post-process an axisymmetric Boozer chartmap from the EQDSK converter
for SIMPLE compatibility.

SIMPLE/libneo requirements:
  - rho: uniform in [0, 1] with rho(1)=0 satisfying rho(1) in [0, 0.01]
  - s (for A_phi): uniform in [rho_min^2, rho_max^2] = [0, 1]
  - theta, zeta: uniform endpoint-excluded grids
  - nzeta >= 2 (Boozer convention)
  - A_phi(s=0) = 0 by definition

Strategy:
  - Use a uniform rho grid with nrho_out points in [0, 1], step = 1/nrho_out.
    With nrho_out=101: step=0.01, rho(1)=0 (satisfies [0,0.01]).
  - Use a separate uniform s grid for A_phi with n_s points in [0, 1].
  - Interpolate R,Z,Bmod,B_theta,B_phi from the 20 converter surfaces.
  - Expand nzeta=1 -> nzeta_out by rotation (axisymmetric tokamak).

Usage: python fix_chartmap_for_simple.py input.nc output.nc [nrho_out [nzeta_out [n_s]]]
"""
import sys
import numpy as np
import netCDF4
from scipy.interpolate import CubicSpline


def fix_chartmap(input_path, output_path, nrho_out=101, nzeta_out=4, n_s=101):
    with netCDF4.Dataset(input_path) as ds:
        rho_in = np.asarray(ds.variables['rho'][:])
        s_in = np.asarray(ds.variables['s'][:])
        theta = np.asarray(ds.variables['theta'][:])
        # NetCDF (zeta, theta, rho) -> Python (rho, theta, zeta)
        x_in = np.transpose(np.asarray(ds.variables['x'][:]), (2, 1, 0))
        y_in = np.transpose(np.asarray(ds.variables['y'][:]), (2, 1, 0))
        z_in = np.transpose(np.asarray(ds.variables['z'][:]), (2, 1, 0))
        Bmod_in = np.transpose(np.asarray(ds.variables['Bmod'][:]), (2, 1, 0))
        A_phi_in = np.asarray(ds.variables['A_phi'][:])
        B_theta_in = np.asarray(ds.variables['B_theta'][:])
        B_phi_in = np.asarray(ds.variables['B_phi'][:])
        attrs = {a: getattr(ds, a) for a in ds.ncattrs()}

    nrho_in, ntheta, _ = x_in.shape

    # Cylindrical coords at zeta=0
    R_in = np.sqrt(x_in[:, :, 0]**2 + y_in[:, :, 0]**2)  # (nrho_in, ntheta)
    Z_in = z_in[:, :, 0]

    # Magnetic axis from linear extrapolation
    drho0 = rho_in[1] - rho_in[0]
    R_axis_arr = R_in[0, :] - rho_in[0] * (R_in[1, :] - R_in[0, :]) / drho0
    Z_axis_arr = Z_in[0, :] - rho_in[0] * (Z_in[1, :] - Z_in[0, :]) / drho0
    R_axis_val = float(np.mean(R_axis_arr))
    Z_axis_val = float(np.mean(Z_axis_arr))
    Bmod_axis_arr = Bmod_in[0, :, 0]

    print(f"Axis: R={R_axis_val:.4f} cm, Z={Z_axis_val:.4f} cm")

    # Prepend axis for spline
    rho_ext = np.concatenate([[0.0], rho_in])
    R_ext = np.vstack([R_axis_val * np.ones(ntheta), R_in])
    Z_ext = np.vstack([Z_axis_val * np.ones(ntheta), Z_in])
    Bmod_ext = np.vstack([Bmod_axis_arr, Bmod_in[:, :, 0]])
    B_theta_ext = np.concatenate([[B_theta_in[0]], B_theta_in])
    B_phi_ext = np.concatenate([[B_phi_in[0]], B_phi_in])

    # Uniform output rho grid [0, 1]
    rho_out = np.linspace(0.0, 1.0, nrho_out)

    # Interpolate geometry and field
    R_out = np.zeros((nrho_out, ntheta))
    Z_out = np.zeros((nrho_out, ntheta))
    Bmod_out = np.zeros((nrho_out, ntheta))

    for j in range(ntheta):
        R_out[:, j] = CubicSpline(rho_ext, R_ext[:, j], extrapolate=True)(rho_out)
        Z_out[:, j] = CubicSpline(rho_ext, Z_ext[:, j], extrapolate=True)(rho_out)
        Bmod_out[:, j] = CubicSpline(rho_ext, Bmod_ext[:, j], extrapolate=True)(rho_out)

    B_theta_out = CubicSpline(rho_ext, B_theta_ext, extrapolate=True)(rho_out)
    B_phi_out = CubicSpline(rho_ext, B_phi_ext, extrapolate=True)(rho_out)

    # Uniform s grid for A_phi: [0, 1]
    s_out = np.linspace(0.0, 1.0, n_s)
    s_ext = np.concatenate([[0.0], s_in])
    A_phi_ext = np.concatenate([[0.0], A_phi_in])
    A_phi_out = CubicSpline(s_ext, A_phi_ext, extrapolate=True)(s_out)
    A_phi_out[0] = 0.0  # exact

    # Expand to nzeta by rotation
    zeta_out = np.linspace(0.0, 2*np.pi, nzeta_out, endpoint=False)
    x_out = np.zeros((nrho_out, ntheta, nzeta_out))
    y_out = np.zeros((nrho_out, ntheta, nzeta_out))
    z_out = np.zeros((nrho_out, ntheta, nzeta_out))
    Bmod_3d = np.zeros((nrho_out, ntheta, nzeta_out))

    for iz, ze in enumerate(zeta_out):
        x_out[:, :, iz] = R_out * np.cos(ze)
        y_out[:, :, iz] = R_out * np.sin(ze)
        z_out[:, :, iz] = Z_out
        Bmod_3d[:, :, iz] = Bmod_out

    print(f"  rho: {rho_out[0]:.6f} to {rho_out[-1]:.6f}, step={rho_out[1]-rho_out[0]:.6f}")
    print(f"  s:   {s_out[0]:.6f} to {s_out[-1]:.6f}, step={s_out[1]-s_out[0]:.6f}")
    print(f"  nrho={nrho_out}, n_s={n_s}, ntheta={ntheta}, nzeta={nzeta_out}")
    print(f"  R: {R_out.min():.1f} to {R_out.max():.1f} cm")
    print(f"  A_phi: {A_phi_out[0]:.4e} to {A_phi_out[-1]:.4e}")

    with netCDF4.Dataset(output_path, 'w', format='NETCDF4') as ds:
        ds.createDimension('rho', nrho_out)
        ds.createDimension('s', n_s)
        ds.createDimension('theta', ntheta)
        ds.createDimension('zeta', nzeta_out)

        v = ds.createVariable('rho', 'f8', ('rho',))
        v[:] = rho_out
        v = ds.createVariable('s', 'f8', ('s',))
        v[:] = s_out
        v = ds.createVariable('theta', 'f8', ('theta',))
        v[:] = theta
        v = ds.createVariable('zeta', 'f8', ('zeta',))
        v[:] = zeta_out

        for name, data in [('x', x_out), ('y', y_out), ('z', z_out)]:
            v = ds.createVariable(name, 'f8', ('zeta', 'theta', 'rho'))
            v.units = 'cm'
            v[:] = np.transpose(data, (2, 1, 0))

        v = ds.createVariable('Bmod', 'f8', ('zeta', 'theta', 'rho'))
        v[:] = np.transpose(Bmod_3d, (2, 1, 0))

        v = ds.createVariable('A_phi', 'f8', ('s',))
        v.radial_abscissa = 's'
        v[:] = A_phi_out

        v = ds.createVariable('B_theta', 'f8', ('rho',))
        v[:] = B_theta_out
        v = ds.createVariable('B_phi', 'f8', ('rho',))
        v[:] = B_phi_out

        v = ds.createVariable('num_field_periods', 'i4')
        v.assignValue(np.int32(1))

        for k, val in attrs.items():
            setattr(ds, k, val)

    print(f"Written: {output_path}")


if __name__ == '__main__':
    inp = sys.argv[1]
    out = sys.argv[2]
    nrho = int(sys.argv[3]) if len(sys.argv) > 3 else 101
    nz = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    ns = int(sys.argv[5]) if len(sys.argv) > 5 else 101
    fix_chartmap(inp, out, nrho, nz, ns)
