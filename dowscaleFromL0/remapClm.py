import sys

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
from remove import z22roms



ncAttribsList = {
    'zeta': {
        'dst_varname': 'zeta',
        'dimensions':  ('ocean_time', 'eta_rho', 'xi_rho'),
        'long_name':   'free-surface',
        'units':      'meter',
        'field':      'free-surface, scalar, series',
        'vartime':     'ocean_time'},
    'temp': {
        'dst_varname': 'temp',
        'dimensions':  ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
        'long_name':   'potential temperature',
        'units':       'Celsius',
        'field':       'temperature, scalar, series',
        'vartime':     'ocean_time'},
    'salt': {
        'dst_varname': 'salt',
        'dimensions':  ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
        'long_name':   'salinity',
        'units':       'PSU',
        'field':       'salinity, scalar, series',
        'vartime':     'ocean_time'},
    'u': {
        'dst_varname': 'u',
        'dimensions': ('ocean_time', 's_rho', 'eta_u', 'xi_u'),
        'long_name':  '3D u-momentum component',
        'units':      'meter second-1',
        'field':      'u-velocity, scalar, series',
        'vartime':       'ocean_time'},
    'ubar': {
        'dst_varname': 'ubar',
        'dimensions':  ('ocean_time', 'eta_u', 'xi_u'),
        'long_name':   '2D u-momentum component',
        'units':       'meter second-1',
        'field':       'ubar-velocity,, scalar, series',
        'vartime':        'ocean_time'},
    'v': {
        'dst_varname': 'v',
        'dimensions': ('ocean_time', 's_rho', 'eta_v', 'xi_v'),
        'long_name':  '3D v-momentum component',
        'units':      'meter second-1',
        'field':      'u-velocity, scalar, series',
        'vartime':       'ocean_time'},
    'vbar': {
        'dst_varname': 'vbar',
        'dimensions':  ('ocean_time', 'eta_v', 'xi_v'),
        'long_name':   '2D v-momentum component',
        'units':       'meter second-1',
        'field':       'ubar-velocity,, scalar, series',
        'vartime':        'ocean_time'}

}

# Dictionary with all the netCDF files that have been already created.
createdFiles = {}


class nctime(object):
    pass

def remapClimate2D(src_file, src_varname, src_grd, dst_grd, oceanTimes,dst_dir='./', idxTime = None):

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'

    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]

    if idxTime is None:
        idxTime = 0
    procTime = cdf.variables['ocean_time'][idxTime]
    print('2D rho-var interpolation of %s at time = %f' % (src_varname, procTime))
    
    # create IC file
    dst_file = src_file.rsplit('/')[-1]
    dst_file = dst_dir + dst_file[:-3] + '_' + src_varname + '_clim_' + dst_grd.name + '.nc'
    if dst_file not in createdFiles:
        print('Creating file', dst_file)
        if os.path.exists(dst_file) is True:
            os.remove(dst_file)
        pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)
        createdFiles[dst_file] = 'done'

    # open IC file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    #get missing value
    spval = src_var._FillValue
    src_var = src_var[0]

    # Check variable dimension
    assert len(src_var.shape) == 2

    pos = 't'
    Cpos = 'rho'
    z = src_grd.vgrid.z_r
    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    try:
        ncAttribs = ncAttribsList[src_varname]
        dst_varname = ncAttribs['dst_varname']
    except:
        print('ERROR: INVALID SOURCE VARIABLE: %s' % src_varname)
        sys.exit(1)

    # create variable in file (if needed)
    if (dst_varname not in nc.variables.keys()):
        print('Creating variable', dst_varname)
        nc.createVariable(dst_varname, 'f8',  ncAttribs['dimensions'])
        nc.variables[dst_varname].long_name = ncAttribs['long_name']
        nc.variables[dst_varname].units     = ncAttribs['units']
        nc.variables[dst_varname].field     = ncAttribs['field']
        nc.variables[dst_varname].time      = ncAttribs['vartime']
        nc.variables['ocean_time'] = oceanTimes

    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, 'to', dst_grd.name)
    print('time =', procTime)

    # horizontal interpolation using xesmf
    print('horizontal interpolation using xesmf')
    dst_var = regrid_GLBy(src_grd, dst_grd, src_var, method='bilinear', fillValue=spval)

    # write data in destination file
    print('write data in destination file')
    nc.variables['ocean_time'][idxTime] = procTime
    nc.variables[dst_varname][idxTime] = dst_var
    
    # close destination file
    nc.close()

    return dst_var


def remapClimate3D(src_file, src_varname, src_grd, dst_grd, oceanTimes, dst_dir='./', idxTime = None, z = None):
    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'

    cdf = netCDF.Dataset(src_file)
    src_var = cdf.variables[src_varname]
    spval = src_var._FillValue

    if idxTime is None:
        idxTime = 0
    procTime = cdf.variables['ocean_time'][idxTime]
    src_var = src_var[idxTime]

    print('3D rho-var interpolation of %s at time = %f' % (src_varname, procTime))


    # create IC file
    dst_file = src_file.rsplit('/')[-1]
    dst_file = dst_dir + dst_file[:-3] + '_' + src_varname + '_clim_' + dst_grd.name + '.nc'
    if dst_file not in createdFiles:
        print('Creating file', dst_file)
        if os.path.exists(dst_file) is True:
            os.remove(dst_file)
        pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)
        createdFiles[dst_file] = 'done'

    # open IC file
    nc = netCDF.Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    # Check variable dimension
    assert len(src_var.shape) == 3

    pos = 't'
    Cpos = 'rho'
    if z is None:
        z = src_grd.vgrid.z_r[:]

    Mp, Lp = dst_grd.hgrid.mask_rho.shape
    try:
        ncAttribs = ncAttribsList[src_varname]
        dst_varname = ncAttribs['dst_varname']
    except:
        print('ERROR: INVALID SOURCE VARIABLE: %s' % src_varname)
        sys.exit(1)

    # create variable in file (if needed)
    if (dst_varname not in nc.variables.keys()):
        print('Creating variable', dst_varname)
        nc.createVariable(dst_varname, 'f8',  ncAttribs['dimensions'])
        nc.variables[dst_varname].long_name = ncAttribs['long_name']
        nc.variables[dst_varname].units     = ncAttribs['units']
        nc.variables[dst_varname].field     = ncAttribs['field']
        nc.variables[dst_varname].time      = ncAttribs['vartime']
        nc.variables['ocean_time'] = oceanTimes


    # build intermediate zgrid, with the same horizontal nodes as L1, but the depths of L0
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, z, z.shape[0])

    # print('>>>>>>>', dst_zcoord.__dict__)

    # theta_s = 8.0
    # theta_b = 0.4
    # Tcline = 20.0
    # L1_N = 15
    # Vstretching = 4
    # dst_zcoord = pyroms.vgrid.s_coordinate_2(dst_grd.vgrid.h, theta_b, theta_s, Tcline, L1_N, zeta=z)

    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)


    # remapping
    print('remapping', dst_varname, 'from', src_grd.name, 'to', dst_grd.name)
    print('time =', procTime)

    # # flood the grid
    # print('flood the grid, spval = ', spval)
    # src_varz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_var, src_grd, pos=pos, spval=spval, \
    #                         dxy=dxy, cdepth=cdepth, kk=kk)


    # horizontal interpolation using xesmf
    print('horizontal interpolation using xesmf')

    dst_varz = regrid_GLBy(src_grd, dst_grd, src_var, method='bilinear', fillValue=spval)
    print('KKKK2', dst_varz.shape)

    # vertical interpolation from standard z level to sigma
    print('vertical interpolation from standard z level to sigma')


    dst_var = z22roms(dst_varz[::-1, :, :], dst_grdz,
                      dst_grd, Cpos=Cpos, spval=spval, flood=False)

    # dst_var = pyroms.remapping.z2roms(dst_varz[::-1, :, :], dst_grdz,
    #                       dst_grd, Cpos=Cpos, spval=spval, flood=False)

    # land mask
    # idx = np.where(dst_grd.hgrid.mask_rho == 0)
    # for n in range(dst_grd.vgrid.N):
        # dst_var[n, idx[0], idx[1]] = np.nan


    # write data in destination file
    print('write data in destination file at time idx = %i (%f)' % (idxTime, procTime))
    nc.variables['ocean_time'][idxTime] = procTime
    nc.variables[dst_varname][idxTime,:,:,:] = dst_var

    # close destination file
    nc.close()


def remapClimateUV(src_file, src_grd, dst_grd, oceanTimes, dst_dir='./', idxTime = None, z = None):
    print('3D velocity interpolation')

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'
    cdf = netCDF.Dataset(src_file)

    if idxTime is None:
        idxTime = 0
    procTime = cdf.variables['ocean_time'][idxTime]


    # create destination file
    dst_file = src_file.rsplit('/')[-1]
    dst_fileu = dst_dir + dst_file[:-3] + '_u_clim_' + dst_grd.name + '.nc'
    if dst_fileu not in createdFiles:
        print('Creating file', dst_fileu)
        if os.path.exists(dst_fileu) is True:
            os.remove(dst_fileu)
        pyroms_toolbox.nc_create_roms_file(dst_fileu, dst_grd, nctime)
        createdFiles[dst_fileu] = 'done'

    dst_filev = dst_dir + dst_file[:-3] + '_v_clim_' + dst_grd.name + '.nc'
    if dst_filev not in createdFiles:
        print('Creating file', dst_filev)
        if os.path.exists(dst_filev) is True:
            os.remove(dst_filev)
        pyroms_toolbox.nc_create_roms_file(dst_filev, dst_grd, nctime)
        createdFiles[dst_filev] = 'done'

    # open destination file
    ncu = netCDF.Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
    ncv = netCDF.Dataset(dst_filev, 'a', format='NETCDF3_64BIT')

    #load var
    cdf = netCDF.Dataset(src_file)
    src_varu = cdf.variables['u']
    src_varv = cdf.variables['v']

    #get missing value
    fillValue = src_varu._FillValue
    src_varu = src_varu[idxTime]
    src_varv = src_varv[idxTime]


    # build intermediate zgrid
    zlevel = -src_grd.vgrid.z_r[:]
    zlevel = zlevel[::-1,0,0]

    # nzlevel = len(zlevel)
    # dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    # dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

    # build intermediate zgrid, with the same horizontal nodes as L1, but the depths of L0
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, z, z.shape[0])
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name + '_Z', dst_grd.hgrid, dst_zcoord)

    # create variables in destination file
    print('Creating variables u, v, ubar, vbar')
    varList = {'u': ncu, 'v': ncv, 'ubar': ncu, 'vbar': ncv}
    for src_varname in varList:
        nc = varList[src_varname]
        try:
            ncAttribs = ncAttribsList[src_varname]
            dst_varname = ncAttribs['dst_varname']
        except:
            print('ERROR: INVALID SOURCE VARIABLE: %s' % src_varname)
            sys.exit(1)

        # create variable in file (if needed)
        if (dst_varname not in nc.variables.keys()):
            print('Creating variable', dst_varname)
            nc.createVariable(dst_varname, 'f8',  ncAttribs['dimensions'])
            nc.variables[dst_varname].long_name = ncAttribs['long_name']
            nc.variables[dst_varname].units =     ncAttribs['units']
            nc.variables[dst_varname].field =     ncAttribs['field']
            nc.variables[dst_varname].time =      ncAttribs['vartime']
            nc.variables['ocean_time'] = oceanTimes


    # remaping
    print('remapping and rotating u and v from', src_grd.name, 'to', dst_grd.name)
    print('time =', procTime)
    print('horizontal interpolation using xesmf')
    dst_uz = regrid_GLBy(src_grd, dst_grd, src_varu, method='bilinear', varType='u', fillValue=fillValue)
    dst_vz = regrid_GLBy(src_grd, dst_grd, src_varv, method='bilinear', varType='v', fillValue=fillValue)

    # plt.imshow(dst_uz[0,:,:])
    # plt.show()

    print('Vertical interpolation from standard z level to sigma')
    dst_u = z22roms(dst_uz[::-1, :, :], dst_grdz, dst_grd, Cpos='rho', spval=fillValue, flood=False)
    dst_v = z22roms(dst_uz[::-1, :, :], dst_grdz, dst_grd, Cpos='rho', spval=fillValue, flood=False)
    # dst_u = pyroms.remapping.z2roms(dst_uz[::-1,:,:], dst_grdz, dst_grd, Cpos='rho', spval=fillValue, flood=False)
    # dst_v = pyroms.remapping.z2roms(dst_vz[::-1,:,:], dst_grdz, dst_grd, Cpos='rho', spval=fillValue, flood=False)


    # plt.imshow(dst_u[0, :, :])
    # plt.show()


    print('Rotating u, v fields')
    print('  Interpolating the angle')
    src_angle = regrid_GLBy(src_grd, dst_grd, src_grd.hgrid.angle_rho, method='bilinear')

    print('  Computing the difference of angles')
    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))

    print('  Actual rotation')
    U = dst_u + dst_v*1j
    eitheta = np.exp(-1j*angle[:,:,:])
    U = U * eitheta
    dst_u = np.real(U)
    dst_v = np.imag(U)

    print('  Move back to U, V points')
    dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
    dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])

    print('Putting FillValue in the masked nodes')
    idxu = (dst_grd.hgrid.mask_u == 0)
    idxv = (dst_grd.hgrid.mask_v == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u[n,idxu] = fillValue
        dst_v[n,idxv] = fillValue


    # compute depth average velocity ubar and vbar
    # get z at the U, V positions
    z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
    z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])

    print('Computes barotropic velocities')
    diffZ_u = np.diff(z_u[:,:,:], 1, 0)
    diffZ_v = np.diff(z_v[:,:,:], 1, 0)
    print(diffZ_u.shape, dst_u.shape)
    print(dst_u[:,1,1].shape, np.diff(z_u[:,1,1]).shape)
    dst_ubar = np.sum(dst_u*diffZ_u, 0) / -z_u[0,:,:]
    dst_vbar = np.sum(dst_v*diffZ_v, 0) / -z_v[0,:,:]

    # for i in range(dst_ubar.shape[1]):
    #     # print(i)
    #     for j in range(dst_ubar.shape[0]):
    #         dst_ubar[j,i] = (dst_u[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]
    #
    # for i in range(dst_vbar.shape[1]):
    #     for j in range(dst_vbar.shape[0]):
    #         dst_vbar[j,i] = (dst_v[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]


    # fillValue
    dst_ubar[idxu] = fillValue
    dst_vbar[idxv] = fillValue

    print('Write data in destination file')
    ncu.variables['ocean_time'][idxTime] = procTime
    ncu.variables['u'][idxTime] = dst_u
    ncu.variables['ubar'][idxTime] = dst_ubar

    ncv.variables['ocean_time'][idxTime] = procTime
    ncv.variables['v'][idxTime] = dst_v
    ncv.variables['vbar'][idxTime] = dst_vbar

    # close destination file
    ncu.close()
    ncv.close()


# def remap_clm_uv(src_file, src_grd, dst_grd, dxy=20, cdepth=0, kk=0, dst_dir='./'):
#
#     # get time
#     nctime.long_name = 'time'
#     nctime.units = 'days since 1900-01-01 00:00:00'
#     time = cdf.variables['ocean_time'][0]
#
#     cdf = netCDF.Dataset(src_file)
#     src_varu = cdf.variables['u']
#     src_varv = cdf.variables['v']
#     spval = src_varu._FillValue
#     src_varu = src_varu[0]
#     src_varv = src_varv[0]
#
#     # create destination file
#     dst_file = src_file.rsplit('/')[-1]
#     dst_fileu = dst_dir + dst_file[:-3] + '_u_clim_' + dst_grd.name + '.nc'
#     print('\nCreating destination file', dst_fileu)
#     if os.path.exists(dst_fileu) is True:
#         os.remove(dst_fileu)
#     pyroms_toolbox.nc_create_roms_file(dst_fileu, dst_grd, nctime)
#
#     dst_filev = dst_dir + dst_file[:-3] + '_v_clim_' + dst_grd.name + '.nc'
#     print('Creating destination file', dst_filev)
#     if os.path.exists(dst_filev) is True:
#         os.remove(dst_filev)
#     pyroms_toolbox.nc_create_roms_file(dst_filev, dst_grd, nctime)
#
#     # open destination file
#     ncu = netCDF.Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
#     ncv = netCDF.Dataset(dst_filev, 'a', format='NETCDF3_64BIT')
#
#     # build intermediate zgrid
#     zlevel = -src_grd.z_t[::-1,0,0]
#     nzlevel = len(zlevel)
#     dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
#     dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)
#
#     # create variable in destination file
#     print('Creating variable u')
#     ncu.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
#     ncu.variables['u'].long_name = '3D u-momentum component'
#     ncu.variables['u'].units = 'meter second-1'
#     ncu.variables['u'].field = 'u-velocity, scalar, series'
#     ncu.variables['u'].time = 'ocean_time'
#     # create variable in destination file
#     print('Creating variable ubar')
#     ncu.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
#     ncu.variables['ubar'].long_name = '2D u-momentum component'
#     ncu.variables['ubar'].units = 'meter second-1'
#     ncu.variables['ubar'].field = 'ubar-velocity,, scalar, series'
#     ncu.variables['ubar'].time = 'ocean_time'
#
#     print('Creating variable v')
#     ncv.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
#     ncv.variables['v'].long_name = '3D v-momentum component'
#     ncv.variables['v'].units = 'meter second-1'
#     ncv.variables['v'].field = 'v-velocity, scalar, series'
#     ncv.variables['v'].time = 'ocean_time'
#     print('Creating variable vbar')
#     ncv.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
#     ncv.variables['vbar'].long_name = '2D v-momentum component'
#     ncv.variables['vbar'].units = 'meter second-1'
#     ncv.variables['vbar'].field = 'vbar-velocity,, scalar, series'
#     ncv.variables['vbar'].time = 'ocean_time'
#
#
#     # remaping
#     print('remapping and rotating u and v from', src_grd.name, 'to', dst_grd.name)
#     print('time =', time)
#
#
#     # flood the grid
#     print('flood the grid')
#     src_uz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varu, src_grd, pos='t', spval=spval, dxy=dxy, cdepth=cdepth, kk=kk)
#     src_vz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varv, src_grd, pos='t', spval=spval, dxy=dxy, cdepth=cdepth, kk=kk)
#
#     # horizontal interpolation using xesmf
#     print('horizontal interpolation using xesmf')
#     dst_uz = regrid_GLBy(src_uz, method='bilinear')
#     dst_vz = regrid_GLBy(src_vz, method='bilinear')
#
#     # vertical interpolation from standard z level to sigma
#     print('vertical interpolation from standard z level to sigma')
#     dst_u = pyroms.remapping.z2roms(dst_uz[::-1,:,:], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=False)
#     dst_v = pyroms.remapping.z2roms(dst_vz[::-1,:,:], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=False)
#
#
#
#     # rotate u,v fields
#     src_angle = regrid_GLBy(src_grd.angle, method='bilinear')
#
#     dst_angle = dst_grd.hgrid.angle_rho
#     angle = dst_angle - src_angle
#     angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
#     U = dst_u + dst_v*1j
#     eitheta = np.exp(-1j*angle[:,:,:])
#     U = U * eitheta
#     dst_u = np.real(U)
#     dst_v = np.imag(U)
#
#
#     # move back to u,v points
#     dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
#     dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])
#
#     # spval
#     idxu = np.where(dst_grd.hgrid.mask_u == 0)
#     idxv = np.where(dst_grd.hgrid.mask_v == 0)
#     for n in range(dst_grd.vgrid.N):
#         dst_u[n,idxu[0], idxu[1]] = spval
#         dst_v[n,idxv[0], idxv[1]] = spval
#
#
#     # compute depth average velocity ubar and vbar
#     # get z at the right position
#     z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
#     z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])
#
#     dst_ubar = np.zeros((dst_u.shape[1], dst_u.shape[2]))
#     dst_vbar = np.zeros((dst_v.shape[1], dst_v.shape[2]))
#
#     for i in range(dst_ubar.shape[1]):
#         for j in range(dst_ubar.shape[0]):
#             dst_ubar[j,i] = (dst_u[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]
#
#     for i in range(dst_vbar.shape[1]):
#         for j in range(dst_vbar.shape[0]):
#             dst_vbar[j,i] = (dst_v[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]
#
#     # spval
#     dst_ubar[idxu[0], idxu[1]] = spval
#     dst_vbar[idxv[0], idxv[1]] = spval
#
#     # write data in destination file
#     print('write data in destination file')
#     ncu.variables['ocean_time'][0] = time
#     ncu.variables['u'][0] = dst_u
#     ncu.variables['ubar'][0] = dst_ubar
#
#     ncv.variables['ocean_time'][0] = time
#     ncv.variables['v'][0] = dst_v
#     ncv.variables['vbar'][0] = dst_vbar
#
#     # close destination file
#     ncu.close()
#     ncv.close()