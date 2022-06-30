import numpy as np
from netCDF4 import Dataset
from setDepth import *
import oceanIn

oceanInFilename = './ocean.in'
gridFileName    = './Sandy_roms_grid.nc'


oceanInput = oceanIn.OceanIn(oceanInFilename)
gridData = Dataset(gridFileName)

h = gridData.variables['h'][:]

Igridtypes = (IGrid.density, IGrid.streamFunc, IGrid.uVel, IGrid.vVel, IGrid.wVel)

z = {}
for igrid in Igridtypes:

    z[igrid] = setDepth(oceanInput.Vtransform, oceanInput.Vstretching, oceanInput.theta_s, oceanInput.theta_b,
                        oceanInput.tcline, oceanInput.N, igrid, h, report = False)

    print('Generated stretched z levels for grid = %s, size %s' % (IGrid.strIGridPointType[igrid-1], repr(z[igrid].shape)))

print('Finished generating stretched z levels')
