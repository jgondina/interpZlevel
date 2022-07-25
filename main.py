import numpy as np
from netCDF4 import Dataset
from setDepth import *
import oceanIn
import scipy.interpolate as spinterp
import ray
import time

Nproc = 15
oceanInFilename = './ocean.in'
gridFileName    = './forecast_snb_his_2022062706.nc'
dataFileName    = './forecast_snb_his_2022062706.nc'
zout = (-4, -3, -2, -1, -0.2)

ray.init() # Initializes the parallel package ray.

t = time.time()

oceanInput = oceanIn.OceanIn(oceanInFilename)
gridData = Dataset(gridFileName)
varData  = Dataset(dataFileName)

vars = {'temp':           IGrid.density,
        'salt':           IGrid.density,
        'ubar_eastward':  IGrid.density,
        'vbar_northward': IGrid.density}

IGridTypes = (IGrid.density, IGrid.streamFunc, IGrid.uVel, IGrid.vVel, IGrid.wVel)


h    = gridData.variables['h'   ][:]
zeta = gridData.variables['zeta'][0,:,:]

# Igridtypes = (IGrid.density, IGrid.streamFunc, IGrid.uVel, IGrid.vVel, IGrid.wVel)

z = {}
v = {}
for igrid in IGridTypes:
    z[igrid] = setDepth(oceanInput.Vtransform, oceanInput.Vstretching, oceanInput.theta_s, oceanInput.theta_b,
                        oceanInput.tcline, oceanInput.N, igrid, h, zeta = zeta, report = False)

    print('Generated stretched z levels for grid = %s, size %s' % (IGrid.strIGridPointType[igrid-1], repr(z[igrid].shape)))

print('Finished generating stretched z levels')


@ray.remote
def procInterp(var_proc, z_igrid_proc, zout):
    res = np.zeros((len(zout), z_igrid_proc.shape[1]))
    for j in range(z_igrid_proc.shape[1]):
        zin = z_igrid_proc[:,j]
        varin = var_proc[:,j]
        f = spinterp.interp1d(zin, varin, bounds_error = False, fill_value = np.nan)
        res[:,j] = f(zout)
    return res

print('Starting vertical interpolation...')
for varName in vars:
    igrid = vars[varName]
    z_igrid = z[igrid]
    var = varData.variables[varName]
    # for i in range(z_igrid.shape[0]):
    #     print(f"Processing \r{(100.0*i)/z_igrid.shape[0]:.1f}%", end="")
    #     for j in range(z_igrid.shape[1]):
    #         zin = z_igrid[i,j,:]
    #         # print(zin)
    #         f = spinterp.interp1d(zin, zin**2, bounds_error = False, fill_value = np.nan)
    #         yout = f(zout)

    for time_idx in range(var.shape[0]):

        v_igrid = np.zeros((len(zout),) + z_igrid.shape[1:])
        v[igrid] = v_igrid
        for i in range(0,z_igrid.shape[1],Nproc):
            print(f"Processing \r{(100.0 * i) / z_igrid.shape[1]:.1f}% using {Nproc} processes", end="")

            N = int(z_igrid.shape[1]/Nproc)
            future = ()
            for i2 in range(i, min(i+Nproc, z_igrid.shape[1])):
                z_igrid_proc = z_igrid[:,i2,:]
                var_proc = var[time_idx,:,i2,:]
                future += (procInterp.remote(var_proc, z_igrid_proc, zout),)
            # z_igrid_proc = z_igrid[i, (Nproc-1)*N:, :]
            # future += (procInterp.remote(z_igrid_proc, zout),)
            j = 0
            for i2 in range(i, min(i+Nproc, z_igrid.shape[1])):
                v_igrid[:,i2,:] = ray.get(future[j])
                j += 1
            # v_igrid[i, (Nproc-1)*N:, :] = ray.get(future[Nproc-1])

print('')
print ('time = %.1f s' % (time.time() - t))


pass