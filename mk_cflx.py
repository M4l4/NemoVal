import netCDF4 as nc
import numpy as np

nc_in =nc.Dataset('adj_uncs_park-ccmp-takahashi_fgco2_1990-2009.nc', 'r')
lon = nc_in.variables['long']
lat = nc_in.variables['lat']
time = nc_in.variables['time']
cflx = nc_in.variables['Cflx']

out_nc = nc.Dataset('uncs_cflx_park_takahasi_1990-2009.nc', 'w')
out_nc.createDimension('lon', len(lon))
out_nc.createDimension('lat', len(lat))
out_nc.createDimension('time', None)

cflx_out = out_nc.createVariable('Cflx', float, ('time', 'lat', 'lon'), fill_value=cflx[0,0,0]*-1)
cflx_out.units = 'mol/m2/s'
cflx_out.long_name = cflx.long_name

lon_out = out_nc.createVariable('lon', float, ('lon'))
lon_out.units = lon.units
lon_out.axis = lon.axis
lon_out.standard_name = lon.standard_name

lat_out = out_nc.createVariable('lat', float, ('lat'))
lat_out.units = lat.units
lat_out.axis = lat.axis
lat_out.standard_name = lat.standard_name


time_out = out_nc.createVariable('time', float, ('time'))
time_out.units = time.units
time_out.standard_name = time.standard_name
time_out.calendar = time.calendar

lon_out[:] = lon[:]
lat_out[:] = lat[:]
time_out[:] = time[:]
cflx_out[:] = cflx[:]*-1
out_nc.close()
