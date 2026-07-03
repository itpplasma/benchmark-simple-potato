"""Expand an axisymmetric Boozer chartmap from nzeta=1 to nzeta>=4.

For a tokamak, the x,y,z geometry at zeta_i is obtained by rotating
the zeta=0 geometry by angle zeta_i around the z-axis.
Bmod is independent of zeta.

Usage: python expand_chartmap_nzeta.py input.nc output.nc [nzeta]
"""
import sys
import numpy as np
import netCDF4

def expand(input_path, output_path, nzeta_out=4):
    with netCDF4.Dataset(input_path) as ds:
        rho = np.asarray(ds.variables['rho'][:])
        s = np.asarray(ds.variables['s'][:])
        theta = np.asarray(ds.variables['theta'][:])
        zeta0 = np.asarray(ds.variables['zeta'][:])  # shape (1,)
        # NetCDF shape is (zeta, theta, rho) -> transpose to (rho, theta, zeta)
        x0 = np.transpose(np.asarray(ds.variables['x'][:]), (2, 1, 0))  # (nrho, ntheta, 1)
        y0 = np.transpose(np.asarray(ds.variables['y'][:]), (2, 1, 0))
        z0 = np.transpose(np.asarray(ds.variables['z'][:]), (2, 1, 0))
        Bmod0 = np.transpose(np.asarray(ds.variables['Bmod'][:]), (2, 1, 0))
        A_phi = np.asarray(ds.variables['A_phi'][:])
        B_theta = np.asarray(ds.variables['B_theta'][:])
        B_phi = np.asarray(ds.variables['B_phi'][:])

        attrs = {a: getattr(ds, a) for a in ds.ncattrs()}
        aphi_attrs = dict(ds.variables['A_phi'].__dict__)

    nrho = len(rho)
    ntheta = len(theta)

    # x,y at zeta=0 give R,Z
    R_map = np.sqrt(x0[:, :, 0]**2 + y0[:, :, 0]**2)  # (nrho, ntheta)
    Z_map = z0[:, :, 0]                                  # (nrho, ntheta)

    # New zeta grid over one full period [0, 2*pi)
    zeta_new = np.linspace(0.0, 2*np.pi, nzeta_out, endpoint=False)

    x_new = np.zeros((nrho, ntheta, nzeta_out))
    y_new = np.zeros((nrho, ntheta, nzeta_out))
    z_new = np.zeros((nrho, ntheta, nzeta_out))
    Bmod_new = np.zeros((nrho, ntheta, nzeta_out))

    for iz, ze in enumerate(zeta_new):
        x_new[:, :, iz] = R_map * np.cos(ze)
        y_new[:, :, iz] = R_map * np.sin(ze)
        z_new[:, :, iz] = Z_map
        Bmod_new[:, :, iz] = Bmod0[:, :, 0]

    with netCDF4.Dataset(output_path, 'w', format='NETCDF4') as ds:
        ds.createDimension('rho', nrho)
        ds.createDimension('s', len(s))
        ds.createDimension('theta', ntheta)
        ds.createDimension('zeta', nzeta_out)

        v = ds.createVariable('rho', 'f8', ('rho',))
        v[:] = rho
        v = ds.createVariable('s', 'f8', ('s',))
        v[:] = s
        v = ds.createVariable('theta', 'f8', ('theta',))
        v[:] = theta
        v = ds.createVariable('zeta', 'f8', ('zeta',))
        v[:] = zeta_new

        for name, data in [('x', x_new), ('y', y_new), ('z', z_new)]:
            v = ds.createVariable(name, 'f8', ('zeta', 'theta', 'rho'))
            v.units = 'cm'
            v[:] = np.transpose(data, (2, 1, 0))

        v = ds.createVariable('Bmod', 'f8', ('zeta', 'theta', 'rho'))
        v[:] = np.transpose(Bmod_new, (2, 1, 0))

        v = ds.createVariable('A_phi', 'f8', ('s',))
        v.radial_abscissa = 's'
        v[:] = A_phi

        v = ds.createVariable('B_theta', 'f8', ('rho',))
        v[:] = B_theta
        v = ds.createVariable('B_phi', 'f8', ('rho',))
        v[:] = B_phi

        v = ds.createVariable('num_field_periods', 'i4')
        v.assignValue(np.int32(1))

        for k, val in attrs.items():
            setattr(ds, k, val)

    print(f"Expanded {input_path} (nzeta=1) -> {output_path} (nzeta={nzeta_out})")


if __name__ == '__main__':
    inp = sys.argv[1]
    out = sys.argv[2]
    nz = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    expand(inp, out, nz)
