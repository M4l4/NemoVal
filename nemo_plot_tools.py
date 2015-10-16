"""Basic library for quick plots of NEMO results. Dirty.

neil.swart@ec.gc.ca, 10/2015
"""
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.pyplot as plt
import matplotlib as mpl
import brewer2mpl
from discrete_cmap import discrete_cmap
#plt.ion()
plt.close('all')
font = {'size'   : 12}
plt.rc('font', **font)

import cmipdata as cd
import cdo; cdo = cdo.Cdo()
from colormaps import viridis

def default_pcolor_args(var):
    """Returns a dict with default pcolor params as key:value pairs"""

    d = {'vmin' : var.min(),
         'vmax' : var.max(),
         'cmap' : viridis,
         }

    return d

def load(ifile, varname, mask_file, maskname):
    """Time average ifile, interpolate to 1x1 degrees, and return var, lon, lat"""
    cdostr=('remapdis,r360x180 -timmean -setctomiss,0 ')
    var = cd.loadvar(ifile, varname, cdostr=cdostr)
    try:
        depth = cd.loadvar(ifile, 'deptht')
    except:
        depth = 0

    lon = np.linspace(0, 359, 360)
    lat = np.linspace(-90,90,180)
    return var, lon, lat, depth

def anom_cmap():
    """return a discrete blue-red cmap from colorbrewer"""
    ncols = 11
    cmap_anom = brewer2mpl.get_map('RdBu', 'diverging', ncols,
                                   reverse=True).mpl_colormap
    cmap_anom = discrete_cmap(ncols, cmap_anom)
    return cmap_anom

def global_map(lon, lat, var, ax=None, vmin=None, vmax=None, cmap=None,
               cblabel=None, title=None, **kwargs):
    """Pcolor a var in a global map, using ax if supplied"""
    # setup a basic global map

    if not ax:
        fig, ax = plt.subplots(1,1, figsize=(8,8))
    else:
        fig = plt.gcf()

    if not vmin:
        vmin = var.min()

    if not vmax:
        vmax = var.max()

    if not cmap:
        cmap = viridis

    if not cblabel:
        cblabel = ''

    if not title:
        title = ''

    m = Basemap(projection='kav7',lon_0=-180,resolution='c', ax=ax)
    lons, lats = np.meshgrid(lon, lat)
    x, y = m(lons, lats)

    cot = m.pcolor(x, y, var, vmin=vmin, vmax=vmax, cmap=cmap,
                   ax=ax, rasterized=True, **kwargs)

    ax.autoscale(enable=True, axis='both', tight=True)
    m.drawcoastlines(linewidth=1.25, ax=ax)
    m.fillcontinents(color='0.8',ax=ax, zorder=2)
    m.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0], linewidth=0,
                    ax=ax)
    #m.drawmeridians(np.arange(0,360,90),labels=[0,0,0,1],
                    #linewidth=0,yoffset=0.5e6, ax=ax)

    m.colorbar(mappable=cot, location='right', label=cblabel)
    ax.set_title(title)

def map_comparison(lon, lat, data1, data2, **kwargs):
    fig, (axl, axm, axr) = plt.subplots(3,1, figsize=(8,8))

    vmax = np.max([data1.max(), data2.max()])
    vmin = np.min([data1.min(), data2.min()])

    global_map(lon, lat, data1, ax=axl, vmin=vmin, vmax=vmax, **kwargs)
    global_map(lon, lat, data2, ax=axm, vmin=vmin, vmax=vmax, **kwargs)
    anom =  data2 - data1
    anom_max = abs(anom).max()
    cm = anom_cmap()
    global_map(lon, lat, anom, ax=axr, cmap=cm, vmin=-1*anom_max, vmax=anom_max,
               **kwargs)

def section(x, z, data, ax=None, ax_args=None, pcolor_args=None, cblabel=''):
    if not ax:
        fig, ax = plt.subplots(1,1, figsize=(8,8))
        fig.subplots_adjust(top=0.8, right=0.8)
    else:
        fig = plt.gcf()

    if not pcolor_args : pcolor_args = default_pcolor_args(data)
    for key, value in default_pcolor_args(data).iteritems():
        if key not in pcolor_args or (pcolor_args[key] is None):
            pcolor_args[key] = value

    cot = ax.pcolormesh(x, z, data, **pcolor_args)
    ax.invert_yaxis()
    ax.autoscale(True, axis='both', tight='both')
    plt.setp(ax, **ax_args)

    box = ax.get_position()
    tl = fig.add_axes([box.x1 + box.width * 0.05, box.y0, 0.02, box.height])
    fig.colorbar(cot, cax=tl, label=cblabel)

if __name__ == '__main__':
    #ifile = ('/raid/ra40/data/ncs/nemo_out/nue/' +
             #'mc_nue_1m_00010101_00011231_grid_t.nc.001')

    #var, lons, lats, depth = load(ifile, 'sosstsst')
    #global_map(lons, lats, var)

    ###
    ifile_y2000 = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                   'mc_nue_1m_20000101_20001231_ptrc_t.nc.001')
    ifile_y0001 = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                   'mc_nue_1m_00010101_00011231_ptrc_t.nc.001')

    no3_y2000, lon, lat, depth = load(ifile_y2000, 'NO3')
    #no3_y0001, lon, lat, depth = load(ifile_y0001, 'NO3')

    #fig, (axl, axm, axr) = plt.subplots(3,1, figsize=(8,8))

    #global_map(lon, lat, no3_y0001[0,:,:], ax=axl, cblabel=r'$\mu$mol l$^{-1}$')
    #global_map(lon, lat, no3_y2000[0,:,:], ax=axm, cblabel=r'$\mu$mol l$^{-1}$')

    #no3_anom =  no3_y2000[0,:,:] - no3_y0001[0,:,:]
    #cm = anom_cmap()
    #global_map(lon, lat, no3_anom, ax=axr, cmap=cm, vmin=-20, vmax=20,
               #cblabel=r'$\mu$mol l$^{-1}$')

    #fig.text(0.4,0.95, 'NUE: ORCA2_CMOC')
    #axl.set_title('Year 0001', fontsize=10)
    #axm.set_title('Year 2000', fontsize=10)
    #axr.set_title('Year 2000 - Year 0001', fontsize=10)
    #plt.savefig('no3_maps_nue.pdf', bbox_inches='tight')

    #map_comparison(lon, lat, no3_y0001[0,:,:], no3_y2000[0,:,:],
                   #cblabel=r'$\mu$mol l$^{-1}$')

    fig, (axl, axm, axr) = plt.subplots(3,1, figsize=(8,8))
    fig.subplots_adjust(right=0.6)

    section(lat, depth, no3_y2000.mean(axis=2), cblabel=r'$\mu$mol l$^{-1}$',
            ax=axl, ax_args={'title':'testing', 'xlim':[-70, 80]})
    plt.show()



