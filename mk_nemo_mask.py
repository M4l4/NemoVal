import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

mesh_nc = nc.Dataset('/raid/ra40/data/ncs/nemo_out/nemo_3.4_orca2_mesh_mask.nc', 'r')
grid_nc = nc.Dataset('/raid/ra40/data/ncs/nemo_out/nue/mc_nue_1m_00010101_00011231_grid_t.nc.001', 'r')

tmask = mesh_nc.variables['tmask']

tmask = mesh_nc.variables['tmask']
x = mesh_nc.dimensions['x']
y = mesh_nc.dimensions['y']
z = mesh_nc.dimensions['z']

e1t = mesh_nc.variables['e1t']
e2t = mesh_nc.variables['e2t']
e3t = mesh_nc.variables['e3t']
gdept = mesh_nc.variables['gdept']

deptht = grid_nc.variables['deptht']
nav_lon = grid_nc.variables['nav_lon']
nav_lat = grid_nc.variables['nav_lat']


tmask_data = tmask[:]
tmask_data[tmask_data==0] = -999

nc4_bathy = nc.Dataset('orca2_tmask.nc', 'w')
nc4_bathy.createDimension('y', len(y))
nc4_bathy.createDimension('x', len(x))
nc4_bathy.createDimension('z', len(deptht))

deptht_out = nc4_bathy.createVariable(
    'deptht', float, ('z'))
deptht_out.units = 'm'

gdept_out = nc4_bathy.createVariable(
    'gdept', float, ('z', 'y', 'x'))
gdept_out.units = 'm'

nav_lat_out = nc4_bathy.createVariable(
    'nav_lat', float, ('y', 'x'))
nav_lat_out.units = 'degrees_north'

nav_lon_out = nc4_bathy.createVariable(
    'nav_lon', float, ('y', 'x'))
nav_lon_out.units = 'degrees_east'

tmask_out = nc4_bathy.createVariable(
    'tmask', float, ('z', 'y', 'x'), fill_value=-999)
tmask_out.coordinates = 'deptht nav_lat nav_lon'
tmask_out.set_auto_mask(True)

e1t_out = nc4_bathy.createVariable(
    'e1t', float, ('y', 'x'))
e1t_out.units = 'm^2'
e1t_out.coordinates = 'nav_lat nav_lon'

e2t_out = nc4_bathy.createVariable(
    'e2t', float, ('y', 'x'))
e2t_out.units = 'm^2'
e2t_out.coordinates = 'nav_lat nav_lon'

e3t_out = nc4_bathy.createVariable(
    'e3t', float, ('z', 'y', 'x'))
e3t_out.units = 'm^3'
e3t_out.coordinates = 'deptht nav_lat nav_lon'

tmask_out[:] = tmask_data.squeeze()
deptht_out[:] = deptht[:].squeeze()
gdept_out[:] = gdept[:].squeeze()
nav_lat_out[:] = nav_lat[:].squeeze()
nav_lon_out[:] = nav_lon[:].squeeze()
e1t_out[:] = e1t[:].squeeze()
e2t_out[:] = e2t[:].squeeze()
e3t_out[:] = e3t[:].squeeze()

grid_nc.close()
mesh_nc.close()
nc4_bathy.close()