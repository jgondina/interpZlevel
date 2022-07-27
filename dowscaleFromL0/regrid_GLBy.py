import xarray as xr
import xesmf
import numpy as np
import matplotlib.pyplot as plt

def regrid_GLBy(src_grd, dst_grd, var, method='nearest_s2d', fillValue = 1e31):

    # srcMask = (np.abs(var - fillValue) > 1e-10) | np.isnan(var)
    #
    plt.imshow(src_grd.hgrid.mask_rho)
    plt.figure()
    plt.imshow(dst_grd.hgrid.mask_rho)
    plt.show()

    # print(np.sum(var[0,:,:][src_grd.hgrid.mask_rho.shape]))

    print(src_grd.hgrid.lat_rho.shape, src_grd.hgrid.mask_rho.shape)
    print(dst_grd.hgrid.lat_rho.shape, dst_grd.hgrid.mask_rho.shape)
    srcCoords = {'lat': src_grd.hgrid.lat_rho, 'lon': src_grd.hgrid.lon_rho, 'mask': src_grd.hgrid.mask_rho.astype(np.int)}
    dstCoords = {'lat': dst_grd.hgrid.lat_rho, 'lon': dst_grd.hgrid.lon_rho, 'mask': dst_grd.hgrid.mask_rho.astype(np.int)}

    regrid = xesmf.Regridder(
        srcCoords, dstCoords,
        method=method,
        # extrap_method = 'nearest_s2d',
        periodic=False,
        filename='regrid_t.nc',
        reuse_weights=False
    )

    # Fills nans and other invalid values.

    var = var[:]
    var2 = np.zeros(var.shape)
    var2[~np.isnan(var)] = var[~np.isnan(var)]
    var2[np.abs(var - fillValue) > 1e-10] = 0.0

    print(var2[:,:].shape)
    tdest = regrid(var)

    return tdest
