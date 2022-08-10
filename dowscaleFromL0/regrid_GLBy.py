import sys
import xesmf
import numpy as np
import matplotlib.pyplot as plt


def regrid_GLBy(src_grd, dst_grd, var, method='nearest_s2d', fillValue = 1e31, varType = 'rho'):

    # Computes the masks
    # srcMask = (np.abs(var - fillValue) < 1e-10*(np.abs(fillValue))) | np.isnan(var) | (src_grd.hgrid.mask_rho == 1.0)
    # dstMask = (np.abs(var - fillValue) < 1e-10*(np.abs(fillValue))) | np.isnan(var) | (src_grd.hgrid.mask_rho == 1.0)

    if varType == 'rho':
        srcCoords = {'lat': src_grd.hgrid.lat_rho, 'lon': src_grd.hgrid.lon_rho, 'mask': src_grd.hgrid.mask_rho.astype(np.int)}
        dstCoords = {'lat': dst_grd.hgrid.lat_rho, 'lon': dst_grd.hgrid.lon_rho, 'mask': 1+0*dst_grd.hgrid.mask_rho.astype(np.int)}
    elif varType == 'u':
        srcCoords = {'lat': src_grd.hgrid.lat_u  , 'lon': src_grd.hgrid.lon_u  , 'mask': src_grd.hgrid.mask_u  .astype(np.int)}
        dstCoords = {'lat': dst_grd.hgrid.lat_rho, 'lon': dst_grd.hgrid.lon_rho, 'mask': dst_grd.hgrid.mask_rho.astype(np.int)}
    elif varType == 'v':
        srcCoords = {'lat': src_grd.hgrid.lat_v  , 'lon': src_grd.hgrid.lon_v  , 'mask': src_grd.hgrid.mask_v  .astype(np.int)}
        dstCoords = {'lat': dst_grd.hgrid.lat_rho, 'lon': dst_grd.hgrid.lon_rho, 'mask': dst_grd.hgrid.mask_rho.astype(np.int)}
    else:
        print('ERROR: Invalid varType. Should be one of rho, u, v')
        sys.exit(1)

    # Computes the regridder. If the weights file does not exist or has any problem, recreates it.
    try:
        regrid = xesmf.Regridder(
            srcCoords, dstCoords,
            method=method,
            extrap_method = 'nearest_s2d',
            periodic=False,
            filename='./regrid_%s.nc' % varType,
            reuse_weights=True)
    except:
        regrid = xesmf.Regridder(
            srcCoords, dstCoords,
            method=method,
            extrap_method = 'nearest_s2d',
            periodic=False,
            filename='./regrid_%s.nc' % varType,
            reuse_weights=False)


    # Converts a possible masked array into a regular one.
    try:
        var = var[:].data
    except:
        pass


    tdest = regrid(var)

    print('KKKK', tdest.shape, tdest[0,:,:].shape)
    print('KKKK', tdest[0, :, :].shape)
    plt.imshow(tdest[0,:,:])
    plt.show()

    # try:
    #     plt.imshow(tdest[:, :])
    # except:
    #     plt.ion()
    #     plt.imshow(tdest[0, :, :],vmin=-1,vmax=1)
    #     plt.figure()
    #     plt.imshow(tdest[10, :, :],vmin=-1,vmax=1)
    #     plt.figure()
    #     plt.ioff()

    # plt.show()


    return tdest
