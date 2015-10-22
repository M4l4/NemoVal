"""Basic library for quick plots of NEMO results. Dirty.

neil.swart@ec.gc.ca, 10/2015
"""
import subprocess
import os
import glob
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.pyplot as plt
import matplotlib as mpl
import brewer2mpl
from discrete_cmap import discrete_cmap
plt.close('all')
font = {'size'   : 12}
plt.rc('font', **font)

from netCDF4 import Dataset
import cdo; cdo = cdo.Cdo()
from colormaps import viridis

def default_pcolor_args(data, anom=False):
    """Returns a dict with default pcolor params as key:value pairs"""

    # Set 3-std range for colorbar to exclude outliers.
    if anom:
        # For anomalies, center range around 0
        anom_max = abs(data).mean() + abs(data).std()*3.0
        vmin = -1*anom_max
        vmax = anom_max
        # Anomaly cmap
        cmap = anom_cmap()
    else:
        # otherwise, center around the mean
        vmin = data.mean() - data.std()*3.0
        vmax = data.mean() + data.std()*3.0
        # Use true min/max if they are closer to the mean than the 3-std spread.
        if vmax > data.max() : vmax = data.max()
        if vmin < data.min() : vmin = data.min()
        # New mpl, colorblind friendly, continuously varying, default cmap
        cmap = viridis

    d = {'vmin' : vmin,
         'vmax' : vmax,
         'cmap' : cmap,
         'rasterized' : True
         }

    return d

def load(ifile_fp, varname):
    """Time average ifile, interpolate to 1x1 degrees, and return var, lon, lat"""

    path, ifile = os.path.split(ifile_fp)
    if not os.path.isfile('remap_' + ifile):
        cdo.remapdis('r360x180', input='-timmean -setctomiss,0 '
                     + ifile_fp, output='remap_' + ifile)

    #var = cd.loadvar(ifile, varname)
    nc = Dataset('remap_' + ifile, 'r' )
    ncvar = nc.variables[varname]
    data = ncvar[:].squeeze()

    try:
        units = ncvar.units
    except:
        units = ''

    try:
        for dimension in ncvar.dimensions:
            if 'depth' in dimension.lower():
                    depth = nc.variables[dimension][:]
                    break
        else:
            depth = 0
    except:
        depth = 0

    lon = np.linspace(0, 359, 360)
    lat = np.linspace(-90,90,180)
    return data, units, lon, lat, depth

def anom_cmap():
    """return a discrete blue-red cmap from colorbrewer"""
    ncols = 11
    cmap_anom = brewer2mpl.get_map('RdBu', 'diverging', ncols,
                                   reverse=True).mpl_colormap
    cmap_anom = discrete_cmap(ncols, cmap_anom)
    return cmap_anom

def global_map(lon, lat, data, ax=None, ax_args=None, pcolor_args=None, cblabel='',
               anom=False):
    """Pcolor a var in a global map, using ax if supplied"""
    # setup a basic global map

    if not ax:
        fig, ax = plt.subplots(1,1, figsize=(8,8))
    else:
        fig = plt.gcf()

    if not pcolor_args : pcolor_args = default_pcolor_args(data, anom=anom)

    for key, value in default_pcolor_args(data, anom=anom).iteritems():
        if key not in pcolor_args or (pcolor_args[key] is None):
            pcolor_args[key] = value

    m = Basemap(projection='kav7',lon_0=-180,resolution='c', ax=ax)
    lons, lats = np.meshgrid(lon, lat)
    x, y = m(lons, lats)

    cot = m.pcolor(x, y, data, **pcolor_args)

    if ax_args:
        plt.setp(ax, **ax_args)

    ax.autoscale(enable=True, axis='both', tight=True)
    m.drawcoastlines(linewidth=1.25, ax=ax)
    m.fillcontinents(color='0.8',ax=ax, zorder=2)
    m.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0], linewidth=0,
                    ax=ax)
    #m.drawmeridians(np.arange(0,360,90),labels=[0,0,0,1],
                    #linewidth=0,yoffset=0.5e6, ax=ax)

    m.colorbar(mappable=cot, location='right', label=cblabel)
    vals = [data.min(), data.max(), data.mean()]
    snam = ['min: ', 'max: ', 'mean: ']
    vals = [s + str(np.round(v,1)) for s, v in zip(snam, vals)]
    x, y = m(10, -88)
    ax.text(x, y, '  '.join(vals), fontsize=8)

def map_comparison(lon, lat, data1, data2, cblabel='', **kwargs):

    # Three panels to plots the obs, data and anomaly
    fig, (axl, axm, axr) = plt.subplots(3,1, figsize=(8,8))

    # compute the anomaly
    if data1.shape != data2.shape:
        raise ValueError('data1.shape != data2.shape, cannot make anomaly')
    else:
        anom =  data2 - data1

    # Check for pcolor and ax args. Assign default pcolor args for anomaly if
    # none have been provided, and for all others set to empty dict.
    for dargs in ['data1_args', 'data2_args', 'anom_args']:
        if (dargs not in kwargs.keys()):
            kwargs[dargs] = {}

        if 'pcolor_args' not in kwargs[dargs].keys():
            kwargs[dargs]['pcolor_args'] =  {}

        if 'ax_args' not in kwargs[dargs].keys():
            kwargs[dargs]['ax_args'] = {}

    global_map(lon, lat, data1, ax=axl,
               pcolor_args=kwargs['data1_args']['pcolor_args'],
               ax_args=kwargs['data1_args']['ax_args'],
               cblabel=cblabel)

    global_map(lon, lat, data2, ax=axm,
               pcolor_args=kwargs['data2_args']['pcolor_args'],
               ax_args=kwargs['data2_args']['ax_args'],
               cblabel=cblabel)

    global_map(lon, lat, anom, ax=axr,
               pcolor_args=kwargs['anom_args']['pcolor_args'],
               ax_args=kwargs['anom_args']['ax_args'],
               cblabel=cblabel, anom=True)

def section(x, z, data, ax=None, ax_args=None, pcolor_args=None, cblabel='',
            contour=True, anom=False):
    if not ax:
        fig, ax = plt.subplots(1,1, figsize=(8,8))
        fig.subplots_adjust(top=0.8, right=0.8)
    else:
        fig = plt.gcf()

    if not pcolor_args : pcolor_args = default_pcolor_args(data, anom=anom)
    for key, value in default_pcolor_args(data, anom=anom).iteritems():
        if key not in pcolor_args or (pcolor_args[key] is None):
            pcolor_args[key] = value

    cot = ax.pcolormesh(x, z, data, **pcolor_args)

    if contour:
        ax.contour(x, z, data, colors=['k'], vmin=pcolor_args['vmin'],
                   vmax=pcolor_args['vmax'])
    ax.invert_yaxis()
    ax.autoscale(True, axis='both', tight='both')
    if ax_args:
        plt.setp(ax, **ax_args)

    box = ax.get_position()
    tl = fig.add_axes([box.x1 + box.width * 0.05, box.y0, 0.02, box.height])
    fig.colorbar(cot, cax=tl, label=cblabel)

def section_comparison(x, z, data1, data2, cblabel='', **kwargs):

    # Three panels to plots the obs, data and anomaly
    fig, (axl, axm, axr) = plt.subplots(3,1, figsize=(8,8), sharex=True,
                                        sharey=True)

    # compute the anomaly
    if data1.shape != data2.shape:
        raise ValueError('data1.shape != data2.shape, cannot make anomaly')
    else:
        anom =  data2 - data1

    # Check for pcolor and ax args. Assign default pcolor args for anomaly if
    # none have been provided, and for all others set to empty dict.
    for dargs in ['data1_args', 'data2_args', 'anom_args']:
        if (dargs not in kwargs.keys()):
            kwargs[dargs] = {}

        if 'pcolor_args' not in kwargs[dargs].keys():
            kwargs[dargs]['pcolor_args'] =  {}

        if 'ax_args' not in kwargs[dargs].keys():
            kwargs[dargs]['ax_args'] = {}

    section(x, z, data1, ax=axl,
               pcolor_args=kwargs['data1_args']['pcolor_args'],
               ax_args=kwargs['data1_args']['ax_args'],
               cblabel=cblabel)

    section(x, z, data2, ax=axm,
               pcolor_args=kwargs['data2_args']['pcolor_args'],
               ax_args=kwargs['data2_args']['ax_args'],
               cblabel=cblabel)

    section(x, z, anom, ax=axr,
               pcolor_args=kwargs['anom_args']['pcolor_args'],
               ax_args=kwargs['anom_args']['ax_args'],
               cblabel=cblabel, anom=True)

    axl.set_xlabel('')
    axm.set_xlabel('')

def proc_plots(plots, obs4comp):
    """Process a list of 'plots'

       Each 'plot' is a dict of key value pairs, defining:

       'ifile' : [str] specifying the input file to load data from

       'variables' : [list] of variables to plot, where each var is a [str]

       'plot_type' : 'section' or 'global_map'

       'pcolor_args' : [dict] option set of key-value pairs passed to pcolor

       'ax_args' : [dict] optional set of key-value pairs used to set axis
                   attributes.

        Output pdfs are saved to plots/ and into a merged pdf called joined.pdf
    """
    plots_out = []
    old_plots = glob.glob('plots/*.pdf')
    for f in old_plots:
        os.remove(f)

    for p in plots:

        # Loop over each variable in the plot, and save a plot for it.
        for v in p['variables']:
            print
            data, units, lon, lat, depth = load(p['ifile'], v)

            if 'pcolor_args' not in p.keys():
                p['pcolor_args'] = None

            #if 'cblabel' not in p.keys():
            p['cblabel'] = units

            # Check is ax_args exists, and if it has a title, if no title set to
            # the varname
            if 'ax_args' not in p.keys():
                p['ax_args'] = {'title' : v}
            else :
                p['ax_args']['title'] = v

            if p['plot_type'] == 'section':
                print 'plotting section of ' + v
                print

                p['ax_args']['xlabel'] = 'Latitude'
                p['ax_args']['xticks'] = np.arange(-80, 81, 20)
                p['ax_args']['ylabel'] = 'Depth'

                try:
                    if data.ndim == 3:
                        zonmean = data.mean(axis=2)
                    elif data.ndim == 2:
                        zonmean = data.mean(axis=1)
                except:
                    'proc_plot cannot zonal mean dfor section ' + ifile + ' ' + v

                section(lat, depth, zonmean, ax_args=p['ax_args'],
                            pcolor_args=p['pcolor_args'], cblabel=p['cblabel'])

                plot_name = 'plots/' + v + '_section.pdf'
                plt.savefig(plot_name, bbox_inches='tight')
                plots_out.append(plot_name)

            if p['plot_type'] == 'global_map':
                print 'plotting global map of ' + v
                print

                if data.ndim > 2:
                    if 'depth' not in p.keys():
                        p['plot_depth'] = np.round(depth.min())
                        print('global_map: plot_depth not specified for ' +
                               v + ', using ' + str(p['plot_depth']))
                    else:
                        p['plot_depth'] = np.round(p['plot_depth'])

                    try:
                        depth_ind = np.where(np.round(depth) == p['plot_depth'])[0][0]
                    except:
                        print('Failed to extract depth ' +  p['plot_depth'] +
                              ' for ' + v)

                    data = data[depth_ind, :, :]

                global_map(lon, lat, data, ax_args=p['ax_args'],
                               pcolor_args=p['pcolor_args'], cblabel=p['cblabel'])
                plot_name = 'plots/' + v + '_map.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

            if (p['plot_type'] == 'global_map_comp') or (p['plot_type'] == 'section_comp'):

                for dargs, atit in zip(['data1_args', 'data2_args', 'anom_args'],
                                       ['Observations', 'NEMO', 'NEMO - Obs.']):
                     if 'kwargs' not in p.keys():
                         p['kwargs'] = {}
                     if (dargs not in p['kwargs'].keys()):
                         p['kwargs'][dargs] = {}
                     if 'ax_args' not in p['kwargs'][dargs]:
                         p['kwargs'][dargs]['ax_args'] = {}

                     p['kwargs'][dargs]['ax_args']['title'] = atit + '(' + v + ')'
                                                        
                     if (p['plot_type'] == 'section_comp'):
                         p['kwargs'][dargs]['ax_args']['xlabel'] = 'Latitude'
                         p['kwargs'][dargs]['ax_args']['xticks'] = np.arange(-80, 81, 20)
                         p['kwargs'][dargs]['ax_args']['ylabel'] = 'Depth'

                if v not in obs4comp.keys():
                    print(v + ' not provided on obs4comp dict, cannot compare')
                    break
                else:
                    (obs_data, obs_units,
                        obs_lon, obs_lat, obs_depth) = load(obs4comp[v], v)

            if (p['plot_type'] == 'global_map_comp'):
                print('plotting global map comparison of ' + v + ' with obs')

                if data.ndim > 2:
                    if 'depth' not in p.keys():
                        p['plot_depth'] = np.round(depth.min())
                        print('global_map: plot_depth not specified for ' +
                               v + ', using ' + str(p['plot_depth']))
                    else:
                        p['plot_depth'] = np.round(p['plot_depth'])

                    try:
                        depth_ind = np.where(np.round(depth) == p['plot_depth'])[0][0]
                        data = data[depth_ind, :, :]
                    except:
                        print('Failed to extract depth ' +  p['plot_depth'] +
                              ' for ' + v)

                if obs_data.ndim == 3:
                     try:
                         ind_obs = np.where(
                         np.round(obs_depth)==p['plot_depth'])[0][0]
                         obs_data = obs_data[ind_obs, :, :]
                     except:
                         print('Failed to extract depth ' + str(p['plot_depth'])
                               + ' for observed ' + v)

                map_comparison(lon, lat, obs_data, data, cblabel=p['cblabel'],
                               **p['kwargs'])
                plot_name = 'plots/' + v + '_map-comp.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

            if (p['plot_type'] == 'section_comp'):
                print('plotting latitudinal section comparison of ' + v + ' with obs')

                if data.ndim > 2:
                    data = data.mean(axis=2)

                if obs_data.ndim > 2:
                    obs_data = obs_data.mean(axis=2)

                section_comparison(lat, depth, obs_data, data,
                                   cblabel=p['cblabel'], **p['kwargs'])
                plot_name = 'plots/' + v + '_section-comp.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

# END of plot loop

    subprocess.Popen(('pdfunite ' + ' '.join(plots_out) +
                      ' plots/joined.pdf'), shell=True).wait()

if __name__ == '__main__':
    ifile_y2000 = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                   'mc_nue_1m_20000101_20001231_ptrc_t.nc.001')
    ifile_y0001 = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                   'mc_nue_1m_00010101_00011231_ptrc_t.nc.001')

    no3_y2000, units, lon, lat, depth = load(ifile_y2000, 'NO3')
    no3_y0001, units, lon, lat, depth = load(ifile_y0001, 'NO3')

    map_comparison(lon, lat, no3_y0001[0,:,:], no3_y2000[0,:,:],
                   cblabel=r'$\mu$mol l$^{-1}$')

    plt.show()



