""" NemoView: A Basic library for quick and dirty plots of NEMO results.

version 0.1

neil.swart@canada.ca, 10/2015
"""
import subprocess
import os
import shutil
import glob
import numpy as np
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.pyplot as plt
import matplotlib as mpl
import brewer2mpl
plt.close('all')
font = {'size'   : 12}
plt.rc('font', **font)

from netCDF4 import Dataset
import cdo; cdo = cdo.Cdo()
from colormaps import viridis
import copy
import taylor
import scipy as sp
from scipy import stats

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    # return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    # https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

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

def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z

def round2SignifFigs(vals,n=0):
    """http://stackoverflow.com/questions/18915378/rounding-to-significant-figures-in-numpy
    """
    import numpy as np
    if np.all(np.isfinite(vals)) and np.all(np.isreal((vals))):
        eset = np.seterr(all='ignore')
        mags = 10.0**np.floor(np.log10(np.abs(vals)))  # omag's
        vals = np.round(vals/mags,n)*mags             # round(val/omag)*omag
        np.seterr(**eset)
    else:
        raise IOError('Input must be real and finite')
    return vals

def load2(ifile, varname, grid_file):
    from scipy.spatial import cKDTree
    lon_new = np.linspace(-179.5,179.5,360)
    lat_new = np.linspace(-89.5,89.5,180)
    tmask_data = tmask[0,0,:,:]
    lon2d, lat2d = np.meshgrid(lon_new, lat_new)
    xn, yn, zn =  lon_lat_to_cartesian(lon2d.flatten(), lat2d.flatten())
    xo, yo, zo =  lon_lat_to_cartesian(nav_lon[:].flatten(), nav_lat[:].flatten())
    tree = cKDTree(zip(xo, yo, zo))
    d, inds = tree.query(zip(xn, yn, zn), k = 1)
    tmask_new = tmask_data.flatten()[inds].reshape(lon2d.shape)

def load(ifile_fp, varname, mask_0=True):
    """Time average ifile, interpolate to 1x1 degrees, and return var, lon, lat"""

    path, ifile = os.path.split(ifile_fp)
    if not os.path.isfile('remap_' + ifile):
        if mask_0:
            cdo.remapdis('r360x180', input='-timmean -setctomiss,0 '
                         + ifile_fp, output='remap_' + ifile)
        else:
            cdo.remapdis('r360x180', input='-timmean '
                         + ifile_fp, output='remap_' + ifile)

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

def corr(obs, data, weights=None):
    """Compute a weighted correlation coefficient and std dev for obs and model

    Returns:
       r : correlation coefficient
       ovar : weighted stddev of obs
       dvar : weighted stddev of data

    """

    if not weights:
        weights = np.ones(obs.shape)

    obsf = obs.flatten()
    dataf = data.flatten()
    weightsf = weights.flatten()

    obar = np.ma.average(obsf, weights=weightsf)
    dbar = np.ma.average(dataf, weights=weightsf)
    ovar = np.sqrt(np.ma.average((obsf-obar)**2, weights=weightsf))
    dvar = np.sqrt(np.ma.average((dataf-dbar)**2, weights=weightsf))

    r = 1.0 / np.nansum(weightsf) * np.nansum( ( (obsf-obar)*(dataf-dbar)*weightsf )
                                        / (ovar*dvar) )

    return r, ovar, dvar


def taylor_plot(data, obs, weights=None, ax_args=None):

    fig = plt.figure(figsize=(8,4))

    print "Computing Taylor diagram statistics..."

    #if not pcolor_args : pcolor_args = default_pcolor_args(data, anom=anom)

    #for key, value in default_pcolor_args(data, anom=anom).iteritems():
    #    if key not in pcolor_args or (pcolor_args[key] is None):
    #        pcolor_args[key] = value

    corrcoef, refstd, stddev = corr(obs, data, weights)
    print "R:", corrcoef
    print "Std, obs, model:", refstd, stddev

    # Taylor diagram
    dia = taylor.TaylorDiagram(refstd, fig=fig, rect=122, label="Obs.")

    #colors = plt.matplotlib.cm.jet(np.linspace(0,1,len(samples)))

    # Add the models to Taylor diagram
    dia.add_sample(stddev, corrcoef, marker='o', ms=10, ls='',
                   mfc='r', mec='r', label="NEMO")

    # Add grid
    dia.add_grid()

    # Add RMS contours, and label them
    contours = dia.add_contours(colors='0.5')
    plt.clabel(contours, inline=1, fontsize=10)

    # Add a figure legend
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ],
               numpoints=1, prop=dict(size='small'), loc='upper right')

    ax = plt.gca()

    if ax_args:
        plt.setp(ax, **ax_args)

def global_map(lon, lat, data, ax=None, ax_args=None, pcolor_args=None, cblabel='',
               level='', anom=False):
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
    if anom:
       vals = [data.min(), data.max(), np.sqrt(np.mean(data**2))]
       snam = ['min: ', 'max: ', 'rmse: ']
    else:
       vals = [data.min(), data.max(), data.mean()]
       snam = ['min: ', 'max: ', 'mean: ']

    vals = [s + str(round2SignifFigs(va,1)) for s, va in zip(snam, vals)]
    x, y = m(10, -88)
    ax.text(x, y, '  '.join(vals), fontsize=8)

    if level:
        x,y = m(0, 90)
        ax.text(x*1.12, y*1.03, level + ' m', fontsize=10)


def map_comparison(lon, lat, data1, data2, cblabel='', level='', **kwargs_in):

    kwargs = copy.deepcopy(kwargs_in)

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

    # If neither of the plots pcolor_args was specified, set the same vmin/max
    if ('vmin' not in kwargs['data1_args']['pcolor_args'] and
        'vmin' not in kwargs['data2_args']['pcolor_args']):
            d1pca = default_pcolor_args(data1)
            d2pca = default_pcolor_args(data2)

            vmin = np.min([d1pca['vmin'], d2pca['vmin']])
            vmax = np.max([d1pca['vmax'], d2pca['vmax']])

            d1pca['vmin'] = vmin
            d1pca['vmax'] = vmax

            kwargs['data1_args']['pcolor_args'] = d1pca
            kwargs['data2_args']['pcolor_args'] = d1pca

    global_map(lon, lat, data1, ax=axl,
               pcolor_args=kwargs['data1_args']['pcolor_args'],
               ax_args=kwargs['data1_args']['ax_args'],
               cblabel=cblabel, level=level)

    global_map(lon, lat, data2, ax=axm,
               pcolor_args=kwargs['data2_args']['pcolor_args'],
               ax_args=kwargs['data2_args']['ax_args'],
               cblabel=cblabel, level=level)

    global_map(lon, lat, anom, ax=axr,
               pcolor_args=kwargs['anom_args']['pcolor_args'],
               ax_args=kwargs['anom_args']['ax_args'],
               cblabel=cblabel, level=level, anom=True)

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

    if anom:
       vals = [data.min(), data.max(), np.sqrt(np.mean(data**2))]
       snam = ['min: ', 'max: ', 'rmse: ']
    else:
       vals = [data.min(), data.max(), data.mean()]
       snam = ['min: ', 'max: ', 'mean: ']

    vals = [s + str(round2SignifFigs(va,1)) for s, va in zip(snam, vals)]
    ylims = ax.get_ylim()
    dy = max(ylims) - min(ylims)
    ypos = max(ylims) - 0.05*dy
    xlims = ax.get_xlim()
    dx = max(xlims) - min(xlims)
    xpos = min(xlims) + 0.05*dx

    ax.text(xpos, ypos, '  '.join(vals),
            fontsize=8, bbox=dict(facecolor='white', alpha=0.75))

    box = ax.get_position()
    tl = fig.add_axes([box.x1 + box.width * 0.05, box.y0, 0.02, box.height])
    fig.colorbar(cot, cax=tl, label=cblabel)

def section_comparison(x, z, data1, data2, cblabel='', **kwargs):

    # Three panels to plots the obs, data and anomaly
    fig, (axl, axm, axr) = plt.subplots(3,1, figsize=(8,8), sharex=True,
                                        sharey=True)

    fig.subplots_adjust(right=0.6)

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

def savefig(fig, plot_name, pdf=True, png=True, **kwargs):

    if pdf:
        fig.savefig(plot_name, bbox_inches='tight', format='pdf', **kwargs)
    if png:
        fig.savefig(plot_name, bbox_inches='tight', format='png', **kwargs)


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

    if os.path.isdir('./plots'):
        print("plots directory exists...overwriting")
        shutil.rmtree('./plots')
    else:
        print("Creating ./plots/ ...")
        print()

    os.mkdir('./plots')

    for p in plots:

        if (p['plot_type'] == 'global_map_comp') or (p['plot_type'] == 'section_comp'):

            for dargs, atit in zip(['data1_args', 'data2_args', 'anom_args'],
                                    ['Observations', 'NEMO', 'NEMO - Obs.']):
                    if 'kwargs' not in p.keys():
                        p['kwargs'] = {}
                    if (dargs not in p['kwargs'].keys()):
                        p['kwargs'][dargs] = {}
                    if 'ax_args' not in p['kwargs'][dargs]:
                        p['kwargs'][dargs]['ax_args'] = {}

                    if (p['plot_type'] == 'section_comp'):
                        p['kwargs'][dargs]['ax_args']['xlabel'] = 'Latitude'
                        p['kwargs'][dargs]['ax_args']['xticks'] = np.arange(-80, 81, 20)
                        p['kwargs'][dargs]['ax_args']['ylabel'] = 'Depth'

        # Loop over each variable in the plot, and save a plot for it.
        for v in p['variables']:

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

                if data.ndim > 2:
                    if 'plot_depth' not in p.keys():
                        p['plot_depth'] = np.round(depth.min())
                        print('global_map: plot_depth not specified for ' +
                               v + ', using ' + str(p['plot_depth']))
                    else:
                        p['plot_depth'] = np.round(p['plot_depth'])

                    try:
                        if p['plot_depth'] in np.round(depth):
                            depth_ind = np.where(np.round(depth) == p['plot_depth'])[0][0]
                        else:
                            print str(p['plot_depth']) + ' not a depth in the NEMO file...looking for closest match instead...'
                            anom = abs(np.round(depth) - p['plot_depth'])
                            depth_ind = np.where(anom == anom.min() )[0][0]
                            p['plot_depth'] = np.round(depth[depth_ind])
                            print '...using ' + str(p['plot_depth']) + ' m \n'

                        data = data[depth_ind, :, :]

                    except:
                        print('Failed to extract depth ' +  p['plot_depth'] +
                              ' for ' + v)

                if 'plot_depth' not in p.keys():
                    p['plot_depth'] = 0

                global_map(lon, lat, data, ax_args=p['ax_args'],
                               pcolor_args=p['pcolor_args'], cblabel=p['cblabel'],
                               level=str(p['plot_depth']))
                plot_name = 'plots/' + v + '_' + str(p['plot_depth']) +'_map.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

            if ( (p['plot_type'] == 'global_map_comp') or
                 (p['plot_type'] == 'section_comp') or
                 (p['plot_type'] == 'taylor_plot')        ):

                if v not in obs4comp.keys():
                    print(v + ' not provided on obs4comp dict, cannot compare')
                    break
                else:
                    (obs_data, obs_units,
                        obs_lon, obs_lat, obs_depth) = load(obs4comp[v], v, mask_0=False)

                if (p['plot_type'] != 'taylor_plot'):
                    for dargs, atit in zip(
                                       ['data1_args','data2_args','anom_args'],
                                       ['Observations', 'NEMO', 'NEMO - Obs.']):
                        p['kwargs'][dargs]['ax_args']['title'] = (atit + ' (' + v +')')
                else:
                    p['ax_args']['title'] = v


            if (p['plot_type'] == 'global_map_comp'):
                print('plotting global map comparison of ' + v + ' with obs')

                if data.ndim > 2:
                    if 'plot_depth' not in p.keys():
                        p['plot_depth'] = np.round(depth.min())
                        print('global_map: plot_depth not specified for ' +
                               v + ', using ' + str(p['plot_depth']))
                    else:
                        p['plot_depth'] = np.round(p['plot_depth'])

                    try:
                        if p['plot_depth'] in np.round(depth):
                            depth_ind = np.where(np.round(depth) == p['plot_depth'])[0][0]
                        else:
                            print 'Specified "plot_depth" of ' + str(p['plot_depth']) + ' m not a depth in the NEMO file...looking for closest match instead...'
                            anom = np.round(depth) - p['plot_depth']
                            depth_ind = np.where( abs(anom) == abs(anom).min())[0][0]
                            p['plot_depth'] = np.round(depth[depth_ind])

                        data = data[depth_ind, :, :]
                    except:
                       print('Failed to extract depth ' +  str(p['plot_depth']) +
                             ' for ' + v)

                elif 'plot_depth' not in p.keys():
                    p['plot_depth'] = 0

                if obs_data.ndim == 3:
                     try:
                         ind_obs = np.where(
                             np.round(obs_depth)==p['plot_depth'])[0][0]
                         obs_data = obs_data[ind_obs, :, :]
                     except:
                         print('Failed to extract depth ' + str(p['plot_depth'])
                               + ' for observed ' + v)

                map_comparison(lon, lat, obs_data, data, cblabel=p['cblabel'],
                               level=str(p['plot_depth']), **p['kwargs'])
                plot_name = 'plots/' + v + '_' + str(p['plot_depth']) + '_map-comp.pdf'
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

            if (p['plot_type'] == 'taylor_plot'):
                taylor_plot(data, obs_data, ax_args = p['ax_args'])
                plot_name = 'plots/' + v + '_taylor.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

# END of plot loop

    subprocess.Popen(('pdfunite ' + ' '.join(plots_out) +
                      ' plots/NemoView_plots.pdf'), shell=True).wait()

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



