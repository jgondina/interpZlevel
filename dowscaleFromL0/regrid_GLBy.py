import xarray as xr
import xesmf
import numpy as np

def regrid_GLBy(src_grd, dst_grd, var, method='nearest_s2d', fillValue = 1e31):

    srcCoords = {'lat': src_grd.hgrid.lat_rho, 'lon': src_grd.hgrid.lon_rho}
    dstCoords = {'lat': dst_grd.hgrid.lat_rho, 'lon': dst_grd.hgrid.lon_rho}

    regrid = xesmf.Regridder(
        srcCoords, dstCoords,
        method=method,
        periodic=False,
        filename='regrid_t.nc',
        reuse_weights=False
    )

    # Fills nans and other invalid values.
    var2 = np.zeros(var.shape)
    var2[~np.isnan(var)] = var[~np.isnan(var)]
    var2[np.abs(var - fillValue) > 1e-10] = 0.0

    print(var)
    tdest = regrid(var2)

    return tdest
