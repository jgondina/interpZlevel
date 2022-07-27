import xarray as xr
import xesmf

def regrid_GLBy(src_grd, dst_grd, fld, method='nearest_s2d'):
    coords = xr.open_dataset('/import/AKWATERS/kshedstrom/gridpak/Arctic2/grid_Arctic_4.nc')
    coords = coords.rename({'lon_rho': 'lon', 'lat_rho': 'lat'})
    gsource = xr.open_dataset('/import/AKWATERS/kshedstrom/HYCOM/Svalbard/data/HYCOM_GLBy0.08_2018_345.nc')

    srcCords = {'lat': src_grd.hgrid.lat_rho, 'lon': src_grd.hgrid.lon_rho}
    dstCords = {'lat': dst_grd.hgrid.lat_rho, 'lon': dst_grd.hgrid.lon_rho}

    regrid = xesmf.Regridder(
        srcCords, dstCoords,
        method=method,
        periodic=False,
        filename='regrid_t.nc',
        reuse_weights=True
    )
    tdest = regrid(fld)
    return tdest