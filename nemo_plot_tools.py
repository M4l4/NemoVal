""" NemoView: A Basic library for quick and dirty plots of NEMO results.

version 0.1

neil.swart@canada.ca, 10/2015
"""
import copy
import os
import shutil
import subprocess

import cdo
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

import taylor

cdo = cdo.Cdo()

plt.close('all')
font = {'size': 12}
plt.rc('font', **font)


def default_pcolor_args(data, anom=False):
    """Returns a dict with default pcolor params as key:value pairs"""

    # Set 3-std range for colorbar to exclude outliers.
    if anom:
        # For anomalies, center range around 0
        anom_max = abs(data).mean() + abs(data).std() * 3.0
        vmin = -1 * anom_max
        vmax = anom_max
        # Anomaly cmap
        cmap = cm.get_cmap("RdBu")
    else:
        # otherwise, center around the mean
        vmin = data.mean() - data.std() * 3.0
        vmax = data.mean() + data.std() * 3.0
        # Use true min/max if they are closer to the mean than the 3-std spread.
        if vmax > data.max():
            vmax = data.max()
        if vmin < data.min():
            vmin = data.min()
        # New mpl, colorblind friendly, continuously varying, default cmap
        cmap = cm.get_cmap("viridis")

    d = {'vmin': vmin,
         'vmax': vmax,
         'cmap': cmap,
         'rasterized': True
         }

    return d


def lon_lat_to_cartesian(lon, lat, r=1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius r
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = r * np.cos(lat_r) * np.cos(lon_r)
    y = r * np.cos(lat_r) * np.sin(lon_r)
    z = r * np.sin(lat_r)
    return x, y, z


def round_to_signif_figs(x, sigfigs):
    """
    Rounds the value(s) in x to the number of significant figures in sigfigs.
    Return value has the same type as x.
    Restrictions:
    sigfigs must be an integer type and store a positive value.
    x must be a real value or an array like object containing only real values.

    http://stackoverflow.com/questions/18915378/rounding-to-significant-figures-in-numpy
    """

    import decimal as decim
    decim.getcontext().prec = 64
    __logBase10of2_decim = decim.Decimal(2).log10()
    __logBase10of2 = float(__logBase10of2_decim)

    if not (type(sigfigs) is int or type(sigfigs) is long or isinstance(sigfigs, np.int)):
        raise TypeError("RoundToSigFigs: sigfigs must be an integer.")

    if sigfigs <= 0:
        raise ValueError("RoundtoSigFigs: sigfigs must be positive.")

    if not np.all(np.isreal(x)):
        raise TypeError("RoundToSigFigs: all x must be real.")

    xsgn = np.sign(x)
    absx = xsgn * x
    mantissas, binary_exponents = np.frexp(absx)

    decimal_exponents = __logBase10of2 * binary_exponents
    omags = np.floor(decimal_exponents)

    mantissas *= 10.0 ** (decimal_exponents - omags)

    if type(mantissas) is float or isinstance(mantissas, np.float):
        if mantissas < 1.0:
            mantissas *= 10.0
            omags -= 1.0

    else:
        fixmsk = mantissas < 1.0
        mantissas[fixmsk] *= 10.0
        omags[fixmsk] -= 1.0

    return xsgn * np.around(mantissas, decimals=sigfigs - 1) * 10.0 ** omags


def load(ifile_fp, varname, mask_0=True):
    """Time average ifile, interpolate to 1x1 degrees, and return var, lon, lat"""

    path, ifile = os.path.split(ifile_fp)
    if not os.path.isfile('remap_new_' + ifile):
        if mask_0:
            cdo.remapdis('r1440x720', input='-timmean -setctomiss,0 ' + ifile_fp, output='remap_new_' + ifile,
                         options='-L')
        else:
            cdo.remapdis('r1440x720', input='-timmean ' + ifile_fp, output='remap_new_' + ifile, options='-L')

    nc = Dataset('remap_new_' + ifile, 'r')
    ncvar = nc.variables[varname]
    data = ncvar[:].squeeze()
    masked_data = ma.masked_values(data, 0, atol=5e-8)

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

    lon = np.linspace(0, 359, 1440)
    lat = np.linspace(-90, 90, 720)
    return masked_data, units, lon, lat, depth


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
    ovar = np.sqrt(np.ma.average((obsf - obar) ** 2, weights=weightsf))
    dvar = np.sqrt(np.ma.average((dataf - dbar) ** 2, weights=weightsf))

    r = 1.0 / np.nansum(weightsf) * np.nansum(((obsf - obar) * (dataf - dbar) * weightsf)
                                              / (ovar * dvar))

    return r, ovar, dvar


def taylor_plot(data, obs, weights=None, ax_args=None):
    fig = plt.figure(figsize=(8, 4))

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
               [p.get_label() for p in dia.samplePoints],
               numpoints=1, prop=dict(size='small'), loc='upper right')

    ax = plt.gca()

    if ax_args:
        plt.setp(ax, **ax_args)


def npolar_map(lon, lat, data, ax=None, ax_args=None, pcolor_args=None, cblabel='', depth='', anom=False):
    """Pcolor a var in a polar map, using ax if supplied"""
    # setup a basic polar map

    if not ax:
        _, ax = plt.subplots(1, 1, figsize=(8, 8))

    if not pcolor_args:
        pcolor_args = default_pcolor_args(data, anom=anom)

    for key, value in default_pcolor_args(data, anom=anom).iteritems():
        if key not in pcolor_args or (pcolor_args[key] is None):
            pcolor_args[key] = value

    m = Basemap(projection='npstere', boundinglat=55, lon_0=0, resolution='l', round=True, ax=ax)

    lons, lats = np.meshgrid(lon, lat)
    x, y = m(lons, lats)

    cot = m.pcolor(x, y, data, **pcolor_args)

    if ax_args:
        plt.setp(ax, **ax_args)

    ax.autoscale(enable=True, axis='both', tight=True)
    m.drawcoastlines(linewidth=1.25, ax=ax)
    m.fillcontinents(color='0.8', ax=ax, zorder=2)
    m.drawparallels(np.arange(40, 80, 10))
    m.drawmeridians(np.arange(0, 360, 30), labels=[1, 0, 0, 1], fontsize=8)

    m.colorbar(mappable=cot, location='right', label=cblabel)
    if anom:
        vals = [data.min(), data.max(), np.sqrt(np.mean(data ** 2))]
        snam = ['min: ', 'max: ', '\nrmse: ']
    else:
        vals = [data.min(), data.max(), data.mean()]
        snam = ['min: ', 'max: ', '\nmean: ']

    vals = [s + str(round_to_signif_figs(va, 2)) for s, va in zip(snam, vals)]
    ax.text(0, 0, '  '.join(vals), fontsize=8, transform=ax.transAxes)

    if depth:
        ax.text(.7, 0, 'Depth: ' + depth + ' m', fontsize=8, transform=ax.transAxes)


def map_comparison(lon, lat, data1, data2, cblabel='', level='', **kwargs_in):
    kwargs = copy.deepcopy(kwargs_in)

    # Three panels to plots the obs, data and anomaly
    fig, (axl, axm, axr) = plt.subplots(3, 1, figsize=(8, 8))

    # compute the anomaly
    try:
        anom = data2 - data1
    except ValueError:
        print 'data1.shape != data2.shape, cannot make anomaly'
        raise

    # Check for pcolor and ax args. Assign default pcolor args for anomaly if
    # none have been provided, and for all others set to empty dict.
    for dargs in ['data1_args', 'data2_args', 'anom_args']:
        if dargs not in kwargs.keys():
            kwargs[dargs] = {}
        if 'pcolor_args' not in kwargs[dargs].keys():
            kwargs[dargs]['pcolor_args'] = {}
        if 'ax_args' not in kwargs[dargs].keys():
            kwargs[dargs]['ax_args'] = {}

    # If neither of the plots pcolor_args was specified, set the same vmin/max
    if 'vmin' not in kwargs['data1_args']['pcolor_args'] and 'vmin' not in kwargs['data2_args']['pcolor_args']:
        d1pca = default_pcolor_args(data1)
        d2pca = default_pcolor_args(data2)

        vmin = np.min([d1pca['vmin'], d2pca['vmin']])
        vmax = np.max([d1pca['vmax'], d2pca['vmax']])

        d1pca['vmin'] = vmin
        d1pca['vmax'] = vmax

        kwargs['data1_args']['pcolor_args'] = d1pca
        kwargs['data2_args']['pcolor_args'] = d1pca

    npolar_map(lon, lat, data1, ax=axl, ax_args=kwargs['data1_args']['ax_args'],
               pcolor_args=kwargs['data1_args']['pcolor_args'], cblabel=cblabel, depth=level)

    npolar_map(lon, lat, data2, ax=axm, ax_args=kwargs['data2_args']['ax_args'],
               pcolor_args=kwargs['data2_args']['pcolor_args'], cblabel=cblabel, depth=level)

    npolar_map(lon, lat, anom, ax=axr, ax_args=kwargs['anom_args']['ax_args'],
               pcolor_args=kwargs['anom_args']['pcolor_args'], cblabel=cblabel, depth=level, anom=True)


def section(x, z, data, ax=None, ax_args=None, pcolor_args=None, cblabel='', contour=True, anom=False):
    if not ax:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        fig.subplots_adjust(top=0.8, right=0.8)
    else:
        fig = plt.gcf()

    if not pcolor_args:
        pcolor_args = default_pcolor_args(data, anom=anom)
    for key, value in default_pcolor_args(data, anom=anom).iteritems():
        if key not in pcolor_args or (pcolor_args[key] is None):
            pcolor_args[key] = value

    cot = ax.pcolormesh(x, z, data, **pcolor_args)

    if contour:
        ax.contour(x, z, data, colors=['k'], vmin=pcolor_args['vmin'], vmax=pcolor_args['vmax'])
    ax.invert_yaxis()
    ax.autoscale(True, axis='both', tight='both')
    if ax_args:
        plt.setp(ax, **ax_args)

    if anom:
        vals = [data.min(), data.max(), np.sqrt(np.mean(data ** 2))]
        snam = ['min: ', 'max: ', 'rmse: ']
    else:
        vals = [data.min(), data.max(), data.mean()]
        snam = ['min: ', 'max: ', 'mean: ']

    vals = [s + str(round_to_signif_figs(va, 2)) for s, va in zip(snam, vals)]
    ylims = ax.get_ylim()
    dy = max(ylims) - min(ylims)
    ypos = max(ylims) - 0.05 * dy
    xlims = ax.get_xlim()
    dx = max(xlims) - min(xlims)
    xpos = min(xlims) + 0.05 * dx

    ax.text(xpos, ypos, '  '.join(vals), fontsize=8, bbox=dict(facecolor='white', alpha=0.75))

    box = ax.get_position()
    tl = fig.add_axes([box.x1 + box.width * 0.05, box.y0, 0.02, box.height])
    fig.colorbar(cot, cax=tl, label=cblabel)


def section_comparison(x, z, data1, data2, cblabel='', **kwargs):
    # Three panels to plots the obs, data and anomaly
    fig, (axl, axm, axr) = plt.subplots(3, 1, figsize=(8, 8), sharex=True, sharey=True)

    fig.subplots_adjust(right=0.6)

    # compute the anomaly
    if data1.shape != data2.shape:
        raise ValueError('data1.shape != data2.shape, cannot make anomaly')
    else:
        anom = data2 - data1

    # Check for pcolor and ax args. Assign default pcolor args for anomaly if
    # none have been provided, and for all others set to empty dict.
    for dargs in ['data1_args', 'data2_args', 'anom_args']:
        if dargs not in kwargs.keys():
            kwargs[dargs] = {}
        if 'pcolor_args' not in kwargs[dargs].keys():
            kwargs[dargs]['pcolor_args'] = {}
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
       'plot_type' : 'section' or 'npolar_map'
       'pcolor_args' : [dict] option set of key-value pairs passed to pcolor
       'ax_args' : [dict] optional set of key-value pairs used to set axis attributes.

        Output pdfs are saved to plots/ and into a merged pdf called joined.pdf
    """
    plots_out = []

    if os.path.isdir('./plots'):
        print("plots directory exists...overwriting")
        shutil.rmtree('./plots')
    else:
        print("Creating ./plots/ ...")

    os.mkdir('./plots')

    for p in plots:
        if (p['plot_type'] == 'npolar_map_comp') or (p['plot_type'] == 'section_comp'):
            for dargs in ['data1_args', 'data2_args', 'anom_args']:
                if 'kwargs' not in p.keys():
                    p['kwargs'] = {}
                if dargs not in p['kwargs'].keys():
                    p['kwargs'][dargs] = {}
                if 'ax_args' not in p['kwargs'][dargs]:
                    p['kwargs'][dargs]['ax_args'] = {}

                if p['plot_type'] == 'section_comp':
                    p['kwargs'][dargs]['ax_args']['xlabel'] = 'Latitude'
                    p['kwargs'][dargs]['ax_args']['xticks'] = np.arange(-80, 81, 20)
                    p['kwargs'][dargs]['ax_args']['ylabel'] = 'Depth'

        # Loop over each variable in the plot, and save a plot for it.
        for x, v in enumerate(p['variables']):
            data, units, lon, lat, depth = load(p['ifile'], v)

            if 'pcolor_args' not in p.keys():
                p['pcolor_args'] = None

            p['cblabel'] = units

            # Check is ax_args exists, and if it has a title, if no title set to the varname
            if 'ax_args' not in p.keys():
                p['ax_args'] = {'title': v}
            else:
                p['ax_args']['title'] = v

            if p['plot_type'] == 'section':
                print 'plotting section of ' + v

                p['ax_args']['xlabel'] = 'Latitude'
                p['ax_args']['xticks'] = np.arange(-80, 81, 20)
                p['ax_args']['ylabel'] = 'Depth'

                if data.ndim == 3:
                    zonmean = data.mean(axis=2)
                elif data.ndim == 2:
                    zonmean = data.mean(axis=1)
                else:
                    print 'proc_plot cannot zonal mean dfor section '  # + ifile + ' ' + v
                    raise ValueError

                section(lat, depth, zonmean, ax_args=p['ax_args'],
                        pcolor_args=p['pcolor_args'], cblabel=p['cblabel'])

                plot_name = 'plots/' + v + '_section.pdf'
                plt.savefig(plot_name, bbox_inches='tight')
                plots_out.append(plot_name)

            if p['plot_type'] == 'npolar_map':
                print 'plotting npolar map of ' + v

                if data.ndim > 2:
                    if 'plot_depth' not in p.keys():
                        p['plot_depth'] = np.round(depth.min())
                        print('global_map: plot_depth not specified for ' + v + ', using ' + str(p['plot_depth']))
                    else:
                        p['plot_depth'] = np.round(p['plot_depth'])

                    try:
                        if p['plot_depth'] in np.round(depth):
                            depth_ind = np.where(np.round(depth) == p['plot_depth'])[0][0]
                        else:
                            print str(p['plot_depth']) + ' not a depth in the NEMO file...looking for closest match' \
                                                         ' instead...'
                            anom = abs(np.round(depth) - p['plot_depth'])
                            depth_ind = np.where(anom == anom.min())[0][0]
                            p['plot_depth'] = np.round(depth[depth_ind])
                            print '...using ' + str(p['plot_depth']) + ' m \n'

                        data = data[depth_ind, :, :]

                    except:
                        print('Failed to extract depth ' + p['plot_depth'] + ' for ' + v)
                        raise

                if 'plot_depth' not in p.keys():
                    p['plot_depth'] = 0

                npolar_map(lon, lat, data, ax_args=p['ax_args'], pcolor_args=p['pcolor_args'], cblabel=p['cblabel'],
                           depth=str(p['plot_depth']))
                plot_name = 'plots/' + v + '_' + str(p['plot_depth']) + '_map.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

            if ((p['plot_type'] == 'npolar_map_comp') or (p['plot_type'] == 'section_comp') or
                    (p['plot_type'] == 'taylor_plot')):
                try:
                    comp_var = p['compare_to'][x]
                except IndexError:
                    print 'No comparison variable provided for ' + v
                    raise
                try:
                    obs_data, _, _, _, obs_depth = load(obs4comp[comp_var], comp_var, mask_0=False)
                except ValueError:
                    print(v + ' not provided on obs4comp dict, cannot compare')
                    raise

                if p['plot_type'] != 'taylor_plot':
                    for dargs, atit in zip(
                            ['data1_args', 'data2_args', 'anom_args'],
                            ['Observations', 'NEMO', 'NEMO - Obs.']):
                        p['kwargs'][dargs]['ax_args']['title'] = (atit + ' (' + v + ')')
                else:
                    p['ax_args']['title'] = v

            if p['plot_type'] == 'npolar_map_comp':
                print('plotting north polar map comparison of ' + v + ' with obs')

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
                            print 'Specified "plot_depth" of ' + str(p['plot_depth']) + ' m not a depth in the NEMO ' \
                                                                                        'file...looking for closest ' \
                                                                                        'match instead...'
                            anom = np.round(depth) - p['plot_depth']
                            depth_ind = np.where(abs(anom) == abs(anom).min())[0][0]
                            p['plot_depth'] = np.round(depth[depth_ind])

                        data = data[depth_ind, :, :]
                    except:
                        print('Failed to extract depth ' + str(p['plot_depth']) + ' for ' + v)
                        raise

                elif 'plot_depth' not in p.keys():
                    p['plot_depth'] = 0

                if obs_data.ndim == 3:
                    try:
                        ind_obs = np.where(
                            np.round(obs_depth) == p['plot_depth'])[0][0]
                        obs_data = obs_data[ind_obs, :, :]
                    except:
                        print('Failed to extract depth ' + str(p['plot_depth']) + ' for observed ' + v)
                        raise

                map_comparison(lon, lat, obs_data, data, cblabel=p['cblabel'],
                               level=str(p['plot_depth']), **p['kwargs'])
                plot_name = 'plots/' + v + '_' + str(p['plot_depth']) + '_map-comp.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

            if p['plot_type'] == 'section_comp':
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

            if p['plot_type'] == 'taylor_plot':
                taylor_plot(data, obs_data, ax_args=p['ax_args'])
                plot_name = 'plots/' + v + '_taylor.pdf'
                plots_out.append(plot_name)
                plt.savefig(plot_name, bbox_inches='tight')

                # END of plot loop

    subprocess.Popen(('pdfunite ' + ' '.join(plots_out) + ' plots/NemoView_plots.pdf'), shell=True).wait()


if __name__ == '__main__':
    ifile_y2000 = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                   'mc_nue_1m_20000101_20001231_ptrc_t.nc.001')
    ifile_y0001 = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                   'mc_nue_1m_00010101_00011231_ptrc_t.nc.001')

    no3_y2000, _, longitude, latitude, depths = load(ifile_y2000, 'NO3')
    no3_y0001, _, longitude, latitude, depths = load(ifile_y0001, 'NO3')

    map_comparison(longitude, latitude, no3_y0001[0, :, :], no3_y2000[0, :, :],
                   cblabel=r'$\mu$mol l$^{-1}$')

    plt.show()
