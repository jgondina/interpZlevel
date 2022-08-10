# This routine creates boundary and initial condition files for ROMS:
# coawst_bdy.nc  coawst_ini.nc on a user-defined grid for a user-defined date.
#
#  This is currently set up to use L0 grid outputs to force L1 grids for the
#  NOPP-NHCI project, based on efforts by:
#  written by Maitane Olabarrieta 09/26/2021

import sys
sys.path.append(r'/home/jo.gonzalez/projects/interpZlevel')

import datetime
import numpy as np
import pyroms
import pyroms_toolbox
import xarray as xr
import xesmf
from remapClm import remapClimate2D, remapClimate3D, remapClimateUV
from setDepth import setDepth




# Select analysis period
timeIni = datetime.datetime(2022, 7, 2, 0, 0, 0)
timeEnd = datetime.datetime(2022, 7, 7, 0, 0, 0)
timeInterval = datetime.timedelta(days=1)

# Select L0 parent grid and output file
L0_grid = r'/orange/olabarrieta/share/BRYINI_example/useast_grd5_2_cnapsv2.nc'
L0_out  = r'/orange/olabarrieta/share/BRYINI_example/noppL0_ocean_his.nc'
L0_Vtransform  = 2
L0_Vstretching = 4
L0_theta_s = 8.0
L0_theta_b = 4.0


# Input ROMS grid to which interpolate operational L0 grid results (L1 grid?)
L1_grid = r'/orange/olabarrieta/share/BRYINI_example/GOMSAB_1km_ext.nc'
L1_out  = r'/blue/olabarrieta/molabarrieta/PROJECTS/NOPP/L1_GOMSAB_1km/ocean_his_00001.nc'


# Enter grid vertical coordinate parameters --These need to be consistent with the ROMS setup.
theta_s = 8.0
theta_b = 0.4
Tcline = 20.0
L1_N = 15
Vtransform  = 2       # Vertical transformation equation
Vstretching = 4       # Vertical stretching function

# Input paths and names of the ini, clm and bry files to generate
init_file = 'L0_L1_GOMSAB_1km_ini.nc'
clm_file  = 'L0_L1_GOMSAB_1km_clm.nc'
bry_file  = 'L0_L1_GOMSAB_1km_bry.nc'

# Extract the characteristics of the child grids
#
# print('Getting roms grid dimensions ...')
#
# modelgridData = Dataset(L1_grid)
#
# h = modelgridData['h']
# hmin = h.min()
# if Vtransform == 1:
#
#     hc = min(max(hmin, 0), Tcline)
# elif Vtransform == 2:
#
#     hc = Tcline
# else:
#     print('ERROR: invalic Vtransform = ', Vtransform)
#     sys.exit(1)


# Time ranges
time   = np.arange(timeIni, timeEnd, datetime.timedelta(hours=3))
time4d = np.arange(timeIni, timeEnd, datetime.timedelta(days=1))


## pyroms.grid.get_ROMS_grid('L0', grid_file = '/orange/olabarrieta/share/BRYINI_example/GOMSAB_1km_ext.nc', hist_file = '/orange/olabarrieta/share/BRYINI_example/noppL0_ocean_his.nc')

# READ L0 GRID INFORMATION AND EXTRACT VARIABLES
gridL1 = pyroms.grid.get_ROMS_grid('L1', grid_file = L1_grid, hist_file = L1_out)
L1_lat_rho = gridL1.hgrid.lat_rho
L1_lat_u   = gridL1.hgrid.lat_u
L1_lat_v   = gridL1.hgrid.lat_v
L1_lon_rho = gridL1.hgrid.lon_rho
L1_lon_u   = gridL1.hgrid.lon_u
L1_lon_v   = gridL1.hgrid.lon_v

nxr, nyr = L1_lat_rho.shape
nxu, nyu = L1_lat_u.  shape
nxv, nyv = L1_lat_v.  shape

L1_u       = np.zeros((nxu, nyu, L1_N, len(time4d)))
L1_v       = np.zeros((nxv, nyv, L1_N, len(time4d)))
L1_temp    = np.zeros((nxr, nyr, L1_N, len(time4d)))
L1_salt    = np.zeros((nxr, nyr, L1_N, len(time4d)))
L1_ubar_4d = np.zeros((nxu, nyu,       len(time4d)))
L1_vbar_4d = np.zeros((nxv, nyv,       len(time4d)))
L1_zeta_4d = np.zeros((nxr, nyr,       len(time4d)))
L1_ubar    = np.zeros((nxu, nyu,       len(time)))
L1_vbar    = np.zeros((nxv, nyv,       len(time)))
L1_zeta    = np.zeros((nxr, nyr,       len(time)))

# Computes the domain ranges (with 0.1 deg of margin).
L1_lon_max = np.max(L1_lon_rho[:]) + 0.1
L1_lon_min = np.min(L1_lon_rho[:]) - 0.1
L1_lat_max = np.max(L1_lat_rho[:]) + 0.1
L1_lat_min = np.min(L1_lat_rho[:]) - 0.1

# READ L0 GRID INFORMATION AND EXTRACT VARIABLES
gridL0 = pyroms.grid.get_ROMS_grid('L0', grid_file =L0_grid, hist_file = L0_out)
L0_lat_rho = gridL0.hgrid.lat_rho
L0_lon_rho = gridL0.hgrid.lon_rho

# Gets the indices of the minimum set of L0 nodes that cover mesh L1 (this is to reduce the number of calculations of the interpolation).
IIr, JJr = np.where((L0_lon_rho<=L1_lon_max) & (L0_lon_rho>=L1_lon_min) &
                    (L0_lat_rho<=L1_lat_max) & (L0_lat_rho>=L1_lat_min))

# Indices limits
L0_xinir = IIr.min()
L0_xendr = IIr.max()
L0_yinir = JJr.min()
L0_yendr = JJr.max()

# Limits for the u and v nodes
L0_xiniu = L0_xinir
L0_xendu = L0_xendr - 1
L0_yiniu = L0_yinir
L0_yendu = L0_yendr
L0_xiniv = L0_xinir
L0_xendv = L0_xendr
L0_yiniv = L0_yinir
L0_yendv = L0_yendr - 1

L0_lonr = L0_lon_rho[L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_latr = L0_lat_rho[L0_xinir:L0_xendr, L0_yinir:L0_yendr]

L0_latr = gridL0.hgrid.lat_rho[L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_lonr = gridL0.hgrid.lon_rho[L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_latu = gridL0.hgrid.lat_u  [L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_lonu = gridL0.hgrid.lon_u  [L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_latv = gridL0.hgrid.lat_v  [L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_lonv = gridL0.hgrid.lon_v  [L0_xinir:L0_xendr, L0_yinir:L0_yendr]

L0_nxr, L0_nyr = L0_latr.shape
L0_nxu, L0_nyu = L0_latu.shape
L0_nxv, L0_nyv = L0_latv.shape




L0_maskr  = gridL0.hgrid.mask_rho [L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_angler = gridL0.hgrid.angle_rho[L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_masku  = gridL0.hgrid.mask_u   [L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_maskv  = gridL0.hgrid.mask_v   [L0_xinir:L0_xendr, L0_yinir:L0_yendr]

L0_h      = gridL0.vgrid.h        [L0_xinir:L0_xendr, L0_yinir:L0_yendr]
L0_hc     = gridL0.vgrid.hc
L0_N      = gridL0.vgrid.s_rho[:].size

L1_h      = gridL1.vgrid.h
L1_hc     = gridL1.vgrid.hc






# for zz in 1:L0_N:
#     L0_lonr_z(:,:,zz)=L0_lonr
#     L0_latr_z(:,:,zz)=L0_latr
#     L0_lonu_z(:,:,zz)=L0_lonu
#     L0_latu_z(:,:,zz)=L0_latu
#     L0_lonv_z(:,:,zz)=L0_lonv
#     L0_latv_z(:,:,zz)=L0_latv
#
#
# ocean_time=double(ncread(L0_out,'ocean_time'))./(3600*24)+datenum(1858,11,17,00,00,00)
# dum=find(ocean_time<=time_end & ocean_time>=time_ini)
# L0_nt=length(dum)
# L0_tini=dum(1)
# L0_tend=dum(end)
# L0_ocean_time=ocean_time(dum)
#
# for zz in 1:L1_N:
#     lon_rho_romsz(:,:,zz)=gn.lon_rho
#     lat_rho_romsz(:,:,zz)=gn.lat_rho
#
#
# tic
#
#
#
#
#
# # -------------------------------------------------------------------------------------------------
# # -------------------------------------------------------------------------------------------------
# # -------------------------------------------------------------------------------------------------
# # -------------------------------------------------------------------------------------------------

oceanTimes = [5000000.0,5000001.0,5000002.0]



for idxTime, time in enumerate(time):
    print('processing time: %s' % time)

    print('Interpolating 2D + time variables')
    L0_zeta_to_L1 = remapClimate2D(L0_out, 'zeta', gridL0, gridL1, oceanTimes, dst_dir='./', idxTime = idxTime)


    # L0_zr = setDepth(L0_Vtransform, L0_Vstretching, L0_theta_s, L0_theta_b, L0_hc, L0_N, 1, L0_h, zeta = L0_zeta)
    # print('ddddddddddd', L0_zr.shape)
    # print('ddddddddddd', L0_zr[:,100,100])
    # print('dddddddd', L0_zr[:, :,:])

    L1_zr = setDepth(Vtransform, Vstretching, theta_s, theta_b, L0_hc, L0_N, 1, L1_h, zeta=L0_zeta_to_L1)


    # L0_UV = remapClimateUV(L0_out, gridL0, gridL1, oceanTimes, dst_dir='./', idxTime = idxTime, z = L1_zr)

    # s_coordinate_2(h, theta_b, theta_s, Tcline, N, hraw=None, zeta=None):


    # L0_zr = setDepth(Vtransform, Vstretching, theta_s, theta_b, L0_hc, L0_N, 1, L0_h, zeta = L0_zeta[L0_xinir:L0_xendr, L0_yinir:L0_yendr])

    # print(':::::::::', L0_zr.shape)

    # L0_zr = setDepth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, zeta=None, report=False):
    print('Interpolating 3D + time variables')
    L0_temp = remapClimate3D(L0_out, 'temp', gridL0, gridL1, oceanTimes, dst_dir='./', idxTime = idxTime, z = L1_zr)
    L0_salt = remapClimate3D(L0_out, 'salt', gridL0, gridL1, oceanTimes, dst_dir='./', idxTime = idxTime, z = L1_zr)


#
#     aa(:,:)=squeeze(L0_zeta(:,:,1))
#     aa(L0_maskr==0)=nan
#     aa(:,:)=maplev(aa)
#     L0_zr(:,:,:)=set_depth(L0_Vtransform,L0_Vstretching,L0_theta_s,L0_theta_b,L0_hc,L0_N, ...
#         1,L0_h,aa)
#     clear aa
#
#
#     for i = 1,
#
#
#     L1_temp(:,:,:,tt)=griddata(L0_lonr_z,L0_latr_z,L0_zr,temp2,lon_rho_romsz,lat_rho_romsz,gn.z_r)
#     clear temp2 %F
#
#     datax = smooth3(L1_temp(:,:,:,tt),'box',5)
#     L1_temp(:,:,:,tt) = datax
#     clear datax
#
#     L1_salt(:,:,:,tt)=griddata(L0_lonr_z,L0_latr_z,L0_zr,salt2,lon_rho_romsz,lat_rho_romsz,gn.z_r)
#     clear salt2 %F
#
#     datax = smooth3(L1_salt(:,:,:,tt),'box',5)
#     L1_salt(:,:,:,tt) = datax
#     clear datax
#
#     for zz=1:L0_N
#         ur(:,:,zz)=u2rho_2d_mw(u2(:,:,zz))
#         vr(:,:,zz)=v2rho_2d_mw(v2(:,:,zz))
#     end
#     clear u2 v2
#
#     % Compute Northward and Eastward velocities
#
#     for zz=1:L0_N
#         vel(:,:)=ur(:,:,zz)+vr(:,:,zz).*sqrt(-1)
#         vel=vel .* exp ( sqrt(-1) * L0_angler)
#         ur(:,:,zz)=real(vel)
#         vr(:,:,zz)=imag(vel)
#     end
#
#     ur2=griddata(L0_lonr_z,L0_latr_z,L0_zr,ur,lon_rho_romsz,lat_rho_romsz,gn.z_r)
#     datax = smooth3(ur2,'box',5)
#     ur2 = datax
#     clear datax
#
#     vr2=griddata(L0_lonr_z,L0_latr_z,L0_zr,vr,lon_rho_romsz,lat_rho_romsz,gn.z_r)
#     datax = smooth3(vr2,'box',5)
#     vr2 = datax
#     clear datax
#
#     clear ur vr
#
#     % Rotate velocities to ROMS grid, important!
#
#     for zz=1:L1_N
#         ur(:,:)=ur2(:,:,zz).*cos(gn.angle)+vr2(:,:,zz).*sin(gn.angle)
#         vr(:,:)=vr2(:,:,zz).*cos(gn.angle)-ur2(:,:,zz).*sin(gn.angle)
#         L1_u(:,:,zz,tt)=rho2u_2d_mw(ur(:,:))  % defined at u points
#         L1_v(:,:,zz,tt)=rho2v_2d_mw(vr(:,:))  % defined at v points
#     end
#
#     clear vel ur vr ur2 vr2 v2 u2
#
#     % Remove possible Nan values
#
#     for zz=1:L1_N
#         aa(:,:)=squeeze(L1_temp(:,:,zz,tt))
#         aa=maplev(aa)
#         L1_temp(:,:,zz,tt)=double(aa)
#         clear aa
#
#         aa(:,:)=squeeze(L1_salt(:,:,zz,tt))
#         aa=maplev(aa)
#         L1_salt(:,:,zz,tt)=double(aa)
#         clear aa
#
#         aa(:,:)=squeeze(L1_u(:,:,zz,tt))
#         aa=maplev(aa)
#         L1_u(:,:,zz,tt)=double(aa)
#         clear aa
#
#         aa(:,:)=squeeze(L1_v(:,:,zz,tt))
#         aa=maplev(aa)
#         L1_v(:,:,zz,tt)=double(aa)
#         clear aa
#     end
# end
#
# toc
# % save('E:\FLORIDA_GRIDS\L1\L0_L1_GOMSAB.mat','time','time4d','L1_temp','L1_salt','L1_u','L1_v','L1_ubar_4d','L1_vbar_4d','L1_zeta_4d','L1_ubar','L1_vbar','L1_zeta','-v7.3')
#


# %% CREATE INITIAL CONDITION
# init_time=time4d(1,1)-datenum(1858,11,17)
# create_roms_init_from_coawst(modelgrid,init_file,init_time,...
#     Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,...
#     Sinp.N,L1_u(:,:,:,1),L1_v(:,:,:,1),L1_ubar(:,:,1),L1_vbar(:,:,1),...
#     L1_temp(:,:,:,1),L1_salt(:,:,:,1),L1_zeta(:,:,1))
#
# ncwrite(init_file,'hc',20.0)
#
# %
# % clear zeta_coawst ubar_coawst vbar_coawst u_coawst v_coawst temp_coawst salt_coawst
# %
# %% CREATE CLIMATOLOGY
# clm_time=time4d-datenum(1858,11,17)
# nt_ts=length(init_time)
#
# create_roms_clm_from_coawst(modelgrid,clm_file,clm_time,...
#     Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,...
#     Sinp.N,L1_u,L1_v,L1_ubar_4d,L1_vbar_4d,L1_temp,L1_salt,L1_zeta_4d)
#
# %% CREATE BOUNDARY CONDITION
# time_bry=time-datenum(1858,11,17)
# time4d_bry=time4d-datenum(1858,11,17)

zeta_north  = L1_zeta[:,-1,:]
ubar_north  = L1_ubar[:,-1,:]
vbar_north  = L1_vbar[:,-1,:]
u_north     = L1_u   [:,-1,:,:]
v_north     = L1_v   [:,-1,:,:]
salt_north  = L1_salt[:,-1,:,:]
temp_north  = L1_temp[:,-1,:,:]

zeta_south  = L1_zeta[:,0,:]
ubar_south  = L1_ubar[:,0,:]
vbar_south  = L1_vbar[:,0,:]
u_sout      = L1_u   [:,0,:,:]
v_sout      = L1_v   [:,0,:,:]
salt_south  = L1_salt[:,0,:,:]
temp_south  = L1_temp[:,0,:,:]

zeta_east   = L1_zeta[-1,:,:]
ubar_east   = L1_ubar[-1,:,:]
vbar_east   = L1_vbar[-1,:,:]
u_east      = L1_u   [-1,:,:,:]
v_east      = L1_v   [-1,:,:,:]
salt_east   = L1_salt[-1,:,:,:]
temp_east   = L1_temp[-1,:,:,:]

zeta_west   = L1_zeta[0,:,:]
ubar_west   = L1_ubar[0,:,:]
vbar_west   = L1_vbar[0,:,:]
u_west      = L1_u   [0,:,:,:]
v_west      = L1_v   [0,:,:,:]
salt_west   = L1_salt[0,:,:,:]
temp_west   = L1_temp[0,:,:,:]
#
# create_roms_bry_from_coawst(modelgrid,bry_file,time_bry,time4d_bry,...
#     Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,Sinp.N,...
#     zeta_north,ubar_north,vbar_north,u_north,v_north,salt_north,temp_north,...
#     ubar_south,vbar_south,zeta_south,u_south,v_south,salt_south,temp_south,...
#     zeta_east,ubar_east,vbar_east,u_east,v_east,salt_east,temp_east,...
#     zeta_west,ubar_west,vbar_west,u_west,v_west,salt_west,temp_west)