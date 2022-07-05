import numpy as np
from netCDF4 import Dataset
from setDepth import *
import oceanIn
import scipy.interpolate as spinterp


oceanInFilename = './ocean.in'
gridFileName    = './forecast_snb_his_2022062706.nc'
zout = (-4, -3, -2, -1, -0.2)


oceanInput = oceanIn.OceanIn(oceanInFilename)
gridData = Dataset(gridFileName)


h    = gridData.variables['h'   ][:]
zeta = gridData.variables['zeta'][0,:,:]

Igridtypes = (IGrid.density, IGrid.streamFunc, IGrid.uVel, IGrid.vVel, IGrid.wVel)

z = {}
for igrid in Igridtypes:

    z[igrid] = setDepth(oceanInput.Vtransform, oceanInput.Vstretching, oceanInput.theta_s, oceanInput.theta_b,
                        oceanInput.tcline, oceanInput.N, igrid, h, zeta = zeta, report = False)

    print('Generated stretched z levels for grid = %s, size %s' % (IGrid.strIGridPointType[igrid-1], repr(z[igrid].shape)))

print('Finished generating stretched z levels')



print('Starting vertical interpolation...')
for igrid in Igridtypes:
    z_igrid = z[igrid]
    for i in range(z_igrid.shape[0]):
        print(f"Processing \r{(100.0*i)/z_igrid.shape[0]:.1f}%", end="")

        for j in range(z_igrid.shape[1]):
            zin = z_igrid[i,j,:]
            # print(zin)
            f = spinterp.interp1d(zin, zin**2, bounds_error = False, fill_value = np.nan)
            yout = f(zout)

pass