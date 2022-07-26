# This routine creates boundary and initial condition files for ROMS:
# coawst_bdy.nc  coawst_ini.nc on a user-defined grid for a user-defined date.
#
#  This is currently set up to use L0 grid outputs to force L1 grids for the
#  NOPP-NHCI project, based on efforts by:
#  written by Maitane Olabarrieta 09/26/2021

import datetime
import numpy as np
import pyroms
# import pyromstools
from netCDF4 import Dataset


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

# Read the input grid.
gridIn = pyroms.grid.get_ROMS_grid('L1', L0_grid, L0_out)
L1_lat_rho = gridIn.hgrid.L1_lat_rho
L1_lat_u   = gridIn.hgrid.L1_lat_u
L1_lat_v   = gridIn.hgrid.L1_lat_v
L1_lon_rho = gridIn.hgrid.L1_lat_rho
L1_lon_u   = gridIn.hgrid.L1_lat_u
L1_lon_v   = gridIn.hgrid.L1_lat_v

nxr, nyr = L1_lat_rho.shape
nxu, nyu = L1_lat_u.  shape
nxv, nyv = L1_lat_v.  shape

L1_u       = np.zeros(nxu, nyu, L1_N, len(time4d))
L1_v       = np.zeros(nxv, nyv, L1_N, len(time4d))
L1_temp    = np.zeros(nxr, nyr, L1_N, len(time4d))
L1_salt    = np.zeros(nxr, nyr, L1_N, len(time4d))
L1_ubar_4d = np.zeros(nxu, nyu, len(time4d))
L1_vbar_4d = np.zeros(nxv, nyv, len(time4d))
L1_zeta_4d = np.zeros(nxr, nyr, len(time4d))
L1_ubar    = np.zeros(nxu, nyu, len(time))
L1_vbar    = np.zeros(nxv, nyv, len(time))
L1_zeta    = np.zeros(nxr, nyr, len(time))

# Computes the domain ranges (with 0.1 deg of margin).
L1_lon_max = max(max(L1_lon_rho)) + 0.1
L1_lon_min = min(min(L1_lon_rho)) - 0.1
L1_lat_max = max(max(L1_lat_rho)) + 0.1
L1_lat_min = min(min(L1_lat_rho)) - 0.1

# READ L0 GRID INFORMATION AND EXTRACT VARIABLES
gridIn = pyroms.grid.get_ROMS_grid('L1', L1_grid, L0_out)


#
#
# lonr=ncread(L0_grid,'lon_rho')
# latr=ncread(L0_grid,'lat_rho')
# [IIr,JJr]=find((lonr<=L1_lon_max & lonr>=L1_lon_min) & (latr<=L1_lat_max & latr>=L1_lat_min))
# L0_xinir=IIr(1)
# L0_xendr=IIr(end)
# L0_yinir=JJr(1)
# L0_yendr=JJr(end)
# L0_lonr=lonr(IIr(1):IIr(end),JJr(1):JJr(end))
# L0_latr=latr(IIr(1):IIr(end),JJr(1):JJr(end))
# [L0_nxr,L0_nyr]=size(L0_lonr)
#
# L0_xiniu=IIr(1)
# L0_xendu=IIr(end)-1
# L0_yiniu=JJr(1)
# L0_yendu=JJr(end)
# L0_xiniv=IIr(1)
# L0_xendv=IIr(end)
# L0_yiniv=JJr(1)
# L0_yendv=JJr(end)-1
#
# clear lonr latr IIr JJr
#
# L0_lonu=ncread(L0_grid,'lon_u',[L0_xinir,L0_yinir],[L0_nxr-1,L0_nyr])
# L0_latu=ncread(L0_grid,'lat_u',[L0_xinir,L0_yinir],[L0_nxr-1,L0_nyr])
# L0_lonv=ncread(L0_grid,'lon_v',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr-1])
# L0_latv=ncread(L0_grid,'lat_v',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr-1])
# [L0_nxu,L0_nyu]=size(L0_lonu)
# [L0_nxv,L0_nyv]=size(L0_lonv)
#
# L0_maskr=ncread(L0_grid,'mask_rho',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr])
# L0_angler=ncread(L0_grid,'angle',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr])
# L0_h=ncread(L0_grid,'h',[L0_xinir,L0_yinir],[L0_nxr,L0_nyr])
# L0_masku=ncread(L0_grid,'mask_u',[L0_xinir,L0_yinir],[L0_nxu,L0_nyu])
# L0_maskv=ncread(L0_grid,'mask_v',[L0_xinir,L0_yinir],[L0_nxv,L0_nyv])
# L0_hc=ncread(L0_out,'hc')
# [L0_N,pp] = size(ncread(L0_out,'s_rho'))
# clear pp
#
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
#
#
#
#
# for tt=1:L0_nt
#
#     disp(tt)
#
#     % Interpolate 3-dimensional variables
#
#     % Free surface elevation
#     L0_zeta=double(ncread(L0_out,'zeta',[L0_xinir,L0_yinir,L0_tini+tt-1],[L0_nxr,L0_nyr,1]))
#     aa(:,:)=squeeze(L0_zeta(:,:,1))
#     aa(L0_maskr==0)=nan
#     aa(:,:)=maplev(aa)
#     zz=griddata(L0_lonr,L0_latr,aa,gn.lon_rho,gn.lat_rho)
#     L1_zeta(:,:,tt)=zz(:,:)
#     clear aa zz
#     display('zeta:',num2str(tt))
#
#     L0_ubar=double(ncread(L0_out,'ubar',[L0_xiniu,L0_yiniu,L0_tini+tt-1],[L0_nxu,L0_nyu,1]))
#     au(:,:)=squeeze(L0_ubar(:,:,1))
#     au(L0_masku==0)=nan
#     au(:,:)=maplev(au)
#     aur=u2rho_2d_mw(au)
#     display('ubar:',num2str(tt))
#
#     L0_vbar=double(ncread(L0_out,'vbar',[L0_xiniv,L0_yiniv,L0_tini+tt-1],[L0_nxv,L0_nyv,1]))
#     av(:,:)=squeeze(L0_vbar(:,:,1))
#     av(L0_maskv==0)=nan
#     av(:,:)=maplev(av)
#     avr=v2rho_2d_mw(av)
#     display('vbar:',num2str(tt))
#     % Compute Northward and Eastward velocities, important!
#
#     vel=aur + avr.*sqrt(-1)
#     vel=vel .* exp ( sqrt(-1) * L0_angler)
#     velu=real(vel)
#     velv=imag(vel)
#     velu1=griddata(L0_lonr,L0_latr,velu,gn.lon_rho,gn.lat_rho)
#     velv1=griddata(L0_lonr,L0_latr,velv,gn.lon_rho,gn.lat_rho)
#
#     % Rotate velocities to ROMS grid, important!
#
#     ubar1(:,:)=velu1.*cos(gn.angle)+velv1.*sin(gn.angle)
#     vbar1(:,:)=velv1.*cos(gn.angle)-velu1.*sin(gn.angle)
#     L1_ubar(:,:,tt)=rho2u_2d_mw(ubar1)  % defined at u points
#     L1_vbar(:,:,tt)=rho2v_2d_mw(vbar1)  % defined at v points
#
#     clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
#     clear L0_ubar L0_vbar
#
# end
# toc
#
# tic
# tt=0
# for tt1=1:8:L0_nt
#     tt=tt+1
#     disp(tt)
#
# % Interpolate 3-dimensional variables
#
#     % Free surface elevation
#     L0_zeta=double(ncread(L0_out,'zeta',[L0_xinir,L0_yinir,L0_tini+tt1-1],[L0_nxr,L0_nyr,1]))
#     aa(:,:)=squeeze(L0_zeta(:,:,1))
#     aa(L0_maskr==0)=nan
#     aa(:,:)=maplev(aa)
#     zz=griddata(L0_lonr,L0_latr,aa,gn.lon_rho,gn.lat_rho)
#     L1_zeta_4d(:,:,tt)=zz(:,:)
#     clear aa zz
#
#     L0_u=double(ncread(L0_out,'u',[L0_xiniu,L0_yiniu,1,L0_tini+tt1-1],[L0_nxu,L0_nyu,inf,1]))
#     L0_v=double(ncread(L0_out,'v',[L0_xiniv,L0_yiniv,1,L0_tini+tt1-1],[L0_nxv,L0_nyv,inf,1]))
#
#     L0_z(:,:,:)=set_depth(L0_Vtransform,L0_Vstretching,L0_theta_s,L0_theta_b,L0_hc,L0_N, ...
#         5,L0_h,squeeze(L0_zeta(:,:,1)))
#
#     L0_z(isnan(L0_z)==1)=0
#
#     for zz=1:L0_N
#         L0_Hz(:,:,zz)=abs(L0_z(:,:,zz+1)-L0_z(:,:,zz))
#     end
#
#     [L0_ubar(:,:),L0_vbar(:,:)]=uv_barotropic(squeeze(L0_u(:,:,:,1)),squeeze(L0_v(:,:,:,1)),squeeze(L0_Hz(:,:,:,1)))
#
#     clear L0_Hz
#
#     au(:,:)=L0_ubar(:,:)
#     au(L0_masku==0)=nan
#     au(:,:)=maplev(au)
#     aur=u2rho_2d_mw(au)
#
#     av(:,:)=L0_vbar(:,:)
#     av(L0_maskv==0)=nan
#     av(:,:)=maplev(av)
#     avr=v2rho_2d_mw(av)
#
#     % Compute Northward and Eastward velocities, important!
#
#     vel=aur + avr.*sqrt(-1)
#     vel=vel .* exp ( sqrt(-1) * L0_angler)
#     velu=real(vel)
#     velv=imag(vel)
#     velu1=griddata(L0_lonr,L0_latr,velu,gn.lon_rho,gn.lat_rho)
#     velv1=griddata(L0_lonr,L0_latr,velv,gn.lon_rho,gn.lat_rho)
#
#     % Rotate velocities to ROMS grid, important!
#
#     ubar1(:,:)=velu1.*cos(gn.angle)+velv1.*sin(gn.angle)
#     vbar1(:,:)=velv1.*cos(gn.angle)-velu1.*sin(gn.angle)
#     L1_ubar_4d(:,:,tt)=rho2u_2d_mw(ubar1)  % defined at u points
#     L1_vbar_4d(:,:,tt)=rho2v_2d_mw(vbar1)  % defined at v points
#
#     clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
#     clear L0_ubar L0_vbar
#
#     % Interpolate 4-dimensional variables
#     L0_temp=double(ncread(L0_out,'temp',[L0_xinir,L0_yinir,1,L0_tini+tt1-1],[L0_nxr,L0_nyr,inf,1]))
#     L0_salt=double(ncread(L0_out,'salt',[L0_xinir,L0_yinir,1,L0_tini+tt1-1],[L0_nxr,L0_nyr,inf,1]))
#
#     temp1(:,:,:)=L0_temp(:,:,:)
#     salt1(:,:,:)=L0_salt(:,:,:)
#     u1(:,:,:)=L0_u(:,:,:)
#     v1(:,:,:)=L0_v(:,:,:)
#
#     for zz=1:L0_N
#         aa(:,:)=squeeze(temp1(:,:,zz))
#         aa(L0_maskr==0)=nan
#         aa=maplev(aa)
#         temp2(:,:,zz)=aa
#         clear aa
#
#         aa(:,:)=squeeze(salt1(:,:,zz))
#         aa(L0_maskr==0)=nan
#         aa=maplev(aa)
#         salt2(:,:,zz)=aa
#         clear aa
#
#         aa(:,:)=squeeze(u1(:,:,zz))
#         aa(L0_masku==0)=nan
#         aa=maplev(aa)
#         u2(:,:,zz)=aa
#         clear aa
#
#         aa(:,:)=squeeze(v1(:,:,zz))
#         aa(L0_maskv==0)=nan
#         aa=maplev(aa)
#         v2(:,:,zz)=aa
#         clear aa
#     end
#
#     clear temp1 salt1 u1 v1
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
# zeta_north=squeeze(L1_zeta(:,end,:))
# ubar_north=squeeze(L1_ubar(:,end,:))
# vbar_north=squeeze(L1_vbar(:,end,:))
# u_north=squeeze(L1_u(:,end,:,:))
# v_north=squeeze(L1_v(:,end,:,:))
# salt_north=squeeze(L1_salt(:,end,:,:))
# temp_north=squeeze(L1_temp(:,end,:,:))
#
# zeta_south=squeeze(L1_zeta(:,1,:))
# ubar_south=squeeze(L1_ubar(:,1,:))
# vbar_south=squeeze(L1_vbar(:,1,:))
# u_south=squeeze(L1_u(:,1,:,:))
# v_south=squeeze(L1_v(:,1,:,:))
# salt_south=squeeze(L1_salt(:,1,:,:))
# temp_south=squeeze(L1_temp(:,1,:,:))
# %
# zeta_east=squeeze(L1_zeta(end,:,:))
# ubar_east=squeeze(L1_ubar(end,:,:))
# vbar_east=squeeze(L1_vbar(end,:,:))
# u_east=squeeze(L1_u(end,:,:,:))
# v_east=squeeze(L1_v(end,:,:,:))
# salt_east=squeeze(L1_salt(end,:,:,:))
# temp_east=squeeze(L1_temp(end,:,:,:))
# %
# zeta_west=squeeze(L1_zeta(1,:,:))
# ubar_west=squeeze(L1_ubar(1,:,:))
# vbar_west=squeeze(L1_vbar(1,:,:))
# u_west=squeeze(L1_u(1,:,:,:))
# v_west=squeeze(L1_v(1,:,:,:))
# salt_west=squeeze(L1_salt(1,:,:,:))
# temp_west=squeeze(L1_temp(1,:,:,:))
#
# create_roms_bry_from_coawst(modelgrid,bry_file,time_bry,time4d_bry,...
#     Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,Sinp.N,...
#     zeta_north,ubar_north,vbar_north,u_north,v_north,salt_north,temp_north,...
#     ubar_south,vbar_south,zeta_south,u_south,v_south,salt_south,temp_south,...
#     zeta_east,ubar_east,vbar_east,u_east,v_east,salt_east,temp_east,...
#     zeta_west,ubar_west,vbar_west,u_west,v_west,salt_west,temp_west)