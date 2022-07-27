import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
import time
from datetime import datetime
from matplotlib.dates import date2num, num2date

import pyroms
import pyroms_toolbox
from pyroms import _remapping
import xarray as xr
import xesmf
import pdb
from regrid_GLBy import regrid_GLBy

class nctime(object):
    pass

def remap_clm(src_file, src_varname, src_grd, dst_grd, dxy=20, cdepth=0, kk=0, dst_dir='./'):

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'

    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname][0]
    time = cdf.variables['ocean_time'][0]


    # create IC file
    dst_file = src_file.rsplit('/')[-1]
    dst_file = dst_dir + dst_file[:-3] + '_' + src_varname + '_clim_' + dst_grd.name + '.nc'
    print('\nCreating file', dst_file)
    if os.path.exists(dst_file) is True:
        os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open IC file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    #load var
    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    #get missing value
    spval = src_var._FillValue
    src_var = src_var[0]

    # determine variable dimension
    ndim = len(src_var.shape)

    pos = 't'
    Cpos = 'rho'
    z = src_grd.vgrid.z_r
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    if src_varname == 'zeta':
        dst_varname = 'zeta'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'free-surface'
        units = 'meter'
        field = 'free-surface, scalar, series'
        vartime = 'ocean_time'
    elif src_varname == 'temp':
        dst_varname = 'temp'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'potential temperature'
        units = 'Celsius'
        field = 'temperature, scalar, series'
        vartime = 'ocean_time'
    elif src_varname == 'salt':
        dst_varname = 'salt'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'salinity'
        units = 'PSU'
        field = 'salinity, scalar, series'
        vartime = 'ocean_time'
    else:
        raise ValueError('Undefined src_varname')


    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1,0,0]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # create variable in file
    print('Creating variable', dst_varname)
    nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    nc.variables[dst_varname].long_name = long_name
    nc.variables[dst_varname].units = units
    nc.variables[dst_varname].field = field
    nc.variables[dst_varname].time = vartime


    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, 'to', dst_grd.name)
    print('time =', time)

    if ndim == 3:
        # flood the grid
        print('flood the grid, spval = ', spval)
        src_varz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_var, src_grd, pos=pos, spval=spval, \
                                dxy=dxy, cdepth=cdepth, kk=kk)
#       flooded = xr.DataArray(src_varz)
#       flooded.to_netcdf("flooded.nc")
    else:
        src_varz = src_var

    return src_varz

    # horizontal interpolation using xesmf
    print('horizontal interpolation using xesmf')

    dst_varz = regrid_GLBy(src_varz, method='bilinear')

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        print('before vertical remap 1', dst_varz[:,928,324])
        print('before vertical remap 2', dst_varz[:,460,350])
        print('vertical interpolation from standard z level to sigma')
        dst_var = pyroms.remapping.z2roms(dst_varz[::-1,:,:], dst_grdz, \
                          dst_grd, Cpos=Cpos, spval=spval, flood=False)
        print('after vertical remap 1', dst_var[:,928,324])
        print('after vertical remap 2', dst_var[:,460,350])

    # land mask
        idxu = np.where(dst_grd.hgrid.mask_rho == 0)
        for n in range(dst_grd.vgrid.N):
            dst_var[n,idxu[0], idxu[1]] = spval
    else:
        dst_var = dst_varz

    # write data in destination file
    print('write data in destination file')
    nc.variables['ocean_time'][0] = time
    nc.variables[dst_varname][0] = dst_var

    # close destination file
    nc.close()

    if src_varname == 'ssh':
        return dst_varz


def remap_clm_uv(src_file, src_grd, dst_grd, dxy=20, cdepth=0, kk=0, dst_dir='./'):

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    cdf = netCDF.Dataset(src_file)
    src_varu = cdf.variables['u'][0]
    src_varv = cdf.variables['v'][0]
    time = cdf.variables['ocean_time'][0]

    # get dimensions
    Mp, Lp = dst_grd.hgrid.mask_rho.shape

    # create destination file
    dst_file = src_file.rsplit('/')[-1]
    dst_fileu = dst_dir + dst_file[:-3] + '_u_clim_' + dst_grd.name + '.nc'
    print('\nCreating destination file', dst_fileu)
    if os.path.exists(dst_fileu) is True:
        os.remove(dst_fileu)
    pyroms_toolbox.nc_create_roms_file(dst_fileu, dst_grd, nctime)
    dst_filev = dst_dir + dst_file[:-3] + '_v_clim_' + dst_grd.name + '.nc'
    print('Creating destination file', dst_filev)
    if os.path.exists(dst_filev) is True:
        os.remove(dst_filev)
    pyroms_toolbox.nc_create_roms_file(dst_filev, dst_grd, nctime)

    # open destination file
    ncu = netCDF.Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
    ncv = netCDF.Dataset(dst_filev, 'a', format='NETCDF3_64BIT')

    #load var
    cdf = netCDF.Dataset(src_file)
    src_varu = cdf.variables['u']
    src_varv = cdf.variables['v']

    #get missing value
    spval = src_varu._FillValue
    src_varu = src_varu[0]
    src_varv = src_varv[0]

    # build intermediate zgrid
    zlevel = -src_grd.z_t[::-1,0,0]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

    # create variable in destination file
    print('Creating variable u')
    ncu.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['u'].long_name = '3D u-momentum component'
    ncu.variables['u'].units = 'meter second-1'
    ncu.variables['u'].field = 'u-velocity, scalar, series'
    ncu.variables['u'].time = 'ocean_time'
    # create variable in destination file
    print('Creating variable ubar')
    ncu.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['ubar'].long_name = '2D u-momentum component'
    ncu.variables['ubar'].units = 'meter second-1'
    ncu.variables['ubar'].field = 'ubar-velocity,, scalar, series'
    ncu.variables['ubar'].time = 'ocean_time'

    print('Creating variable v')
    ncv.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['v'].long_name = '3D v-momentum component'
    ncv.variables['v'].units = 'meter second-1'
    ncv.variables['v'].field = 'v-velocity, scalar, series'
    ncv.variables['v'].time = 'ocean_time'
    print('Creating variable vbar')
    ncv.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['vbar'].long_name = '2D v-momentum component'
    ncv.variables['vbar'].units = 'meter second-1'
    ncv.variables['vbar'].field = 'vbar-velocity,, scalar, series'
    ncv.variables['vbar'].time = 'ocean_time'


    # remaping
    print('remapping and rotating u and v from', src_grd.name, \
                      'to', dst_grd.name)
    print('time =', time)


    # flood the grid
    print('flood the grid')
    src_uz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varu, src_grd, pos='t', \
                spval=spval, dxy=dxy, cdepth=cdepth, kk=kk)
    src_vz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varv, src_grd, pos='t', \
                spval=spval, dxy=dxy, cdepth=cdepth, kk=kk)

    # horizontal interpolation using xesmf
    print('horizontal interpolation using xesmf')
    dst_uz = regrid_GLBy(src_uz, method='bilinear')
    dst_vz = regrid_GLBy(src_vz, method='bilinear')

    # vertical interpolation from standard z level to sigma
    print('vertical interpolation from standard z level to sigma')
    print('before vertical remap u', dst_uz[:,928,324])
    print('before vertical remap v', dst_vz[:,928,324])
    dst_u = pyroms.remapping.z2roms(dst_uz[::-1,:,:], dst_grdz, \
                        dst_grd, Cpos='rho', spval=spval, flood=False)
    dst_v = pyroms.remapping.z2roms(dst_vz[::-1,:,:], dst_grdz,  \
                        dst_grd, Cpos='rho', spval=spval, flood=False)
    print('after vertical remap u', dst_uz[:,928,324])
    print('after vertical remap v', dst_vz[:,928,324])


    # rotate u,v fields
    src_angle = regrid_GLBy(src_grd.angle, method='bilinear')

    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
    U = dst_u + dst_v*1j
    eitheta = np.exp(-1j*angle[:,:,:])
    U = U * eitheta
    dst_u = np.real(U)
    dst_v = np.imag(U)


    # move back to u,v points
    dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
    dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])

    # spval
    idxu = np.where(dst_grd.hgrid.mask_u == 0)
    idxv = np.where(dst_grd.hgrid.mask_v == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u[n,idxu[0], idxu[1]] = spval
        dst_v[n,idxv[0], idxv[1]] = spval


    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
    z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])

    dst_ubar = np.zeros((dst_u.shape[1], dst_u.shape[2]))
    dst_vbar = np.zeros((dst_v.shape[1], dst_v.shape[2]))

    for i in range(dst_ubar.shape[1]):
        for j in range(dst_ubar.shape[0]):
            dst_ubar[j,i] = (dst_u[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]

    for i in range(dst_vbar.shape[1]):
        for j in range(dst_vbar.shape[0]):
            dst_vbar[j,i] = (dst_v[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]

    # spval
    dst_ubar[idxu[0], idxu[1]] = spval
    dst_vbar[idxv[0], idxv[1]] = spval

    # write data in destination file
    print('write data in destination file')
    ncu.variables['ocean_time'][0] = time
    ncu.variables['u'][0] = dst_u
    ncu.variables['ubar'][0] = dst_ubar

    ncv.variables['ocean_time'][0] = time
    ncv.variables['v'][0] = dst_v
    ncv.variables['vbar'][0] = dst_vbar

    # close destination file
    ncu.close()
    ncv.close()