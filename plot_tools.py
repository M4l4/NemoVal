# coding=utf-8
""" PlotTools: A Basic library for quick and dirty plots of model results.

version 0.2

malachitehmiller@gmail.com, 03/2017
"""
import copy
import os
import shutil
import subprocess

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset, num2date

import taylor

plt.close('all')
font = {'size': 12}
plt.rc('font', **font)
mpl.rcParams['image.cmap'] = 'viridis'


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


def load(ifile_fp, varname):

    path, ifile = os.path.split(ifile_fp)

    if not os.path.isfile('remapped/intlev_{}'.format(ifile)):
        try:
            import cdo
            cdo = cdo.Cdo()
            cdo.env = {'REMAP EXTRAPOLATE': 'off'}
        except ImportError:
            print 'CDO must be installed to preproces files.'
            raise

        print 'Preprocessing', ifile

        with open('parameters/default_depths') as f:
            default_depths = f.read().replace('\n', ',')

        try:
            with Dataset('parameters/masks.nc', 'r') as mask:
                with Dataset(ifile_fp, 'r') as data:
                    if np.all(np.array(mask['nav_lat'][:]) == np.array(data['nav_lat'][:])) and \
                       np.all(np.array(mask['nav_lon'][:]) == np.array(data['nav_lon'][:])):
                        cdo.yearmean(input=ifile_fp, output='temp_' + ifile)
                    else:
                        raise AttributeError
                with Dataset('temp_' + ifile, 'a') as to_mask:
                    for var in to_mask.variables:
                        if len(to_mask[var].shape) == 4:
                            for i in range(0, to_mask[var].shape[0]):
                                to_mask[var][i, :, :, :] = ma.masked_where(
                                    np.logical_not(np.array(mask['tmask'][0, :, :, :], dtype=bool)),
                                    np.array(to_mask[var][i, :, :, :]))[:]
            cdo.intlevel(default_depths, input='-remapdis,parameters/default_grid temp_' + ifile,
                         output='remapped/intlev_' + ifile, options='-L', force=False)
            os.remove('temp_' + ifile)
        except (AttributeError, IndexError):
            print 'Data from {} is not on the NAA grid, no land mask applied.'.format(ifile)
            # Year average ifile, and interpolate to 1/2x1/2 degrees
            cdo.intlevel(default_depths,
                         input='-remapdis,parameters/default_grid -yearmean -sellonlatbox,-180,180,40,90 ' + ifile_fp,
                         output='remapped/intlev_' + ifile, options='-L', force=False)
        print ifile, 'preprocessed.'

    remaped_file = 'remapped/intlev_' + ifile
    nc = Dataset(remaped_file, 'r')
    ncvar = nc.variables[varname]
    data = ncvar[:].squeeze()
    masked_data = ma.masked_values(data, 0, atol=5e-8)

    try:
        units = ncvar.units
    except AttributeError:
        print 'Units not given for {} in {}, leaving empty.'.format(varname, ifile)
        units = ''

    dimensions = ncvar.dimensions

    try:
        for dimension in dimensions:
            if 'depth' in dimension.lower():
                depth = nc.variables[dimension][:]
                break
        else:
            raise IndexError
    except IndexError:
        print '\nDepths not given for {} in {}.'.format(varname, ifile)
        raise

    lon = np.linspace(0, 360, 721)
    lat = np.linspace(40, 90, 101)

    dates = num2date(nc['time_counter'][:], units=nc['time_counter'].units, calendar=nc['time_counter'].calendar)
    years = [x.year for x in dates]

    return masked_data, units, lon, lat, depth, dimensions, years


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

    corrcoef, refstd, stddev = corr(obs, data, weights)
    print "R:", corrcoef
    print "Std, obs, model:", refstd, stddev

    # Taylor diagram
    dia = taylor.TaylorDiagram(refstd, fig=fig, rect=122, label="Obs.")

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


def npolar_map(lon, lat, data, ax=None, ax_args=None, pcolor_args=None, cblabel='', depth=None, anom=False):
    """Pcolor a var in a polar map, using ax if supplied"""
    # setup a basic polar map

    if data.max() is ma.masked:
        raise ValueError

    if not ax:
        _, ax = plt.subplots(1, 1, figsize=(8, 8))

    if not pcolor_args:
        pcolor_args = default_pcolor_args(data, anom=anom)

    # for key, value in default_pcolor_args(data, anom=anom).iteritems():
    #     if key not in pcolor_args or (pcolor_args[key] is None):
    #         pcolor_args[key] = value

    m = Basemap(projection='npstere', boundinglat=55, lon_0=0, resolution='l', round=True, ax=ax)

    lons, lats = np.meshgrid(lon, lat)
    x, y = m(lons, lats)
    graphed_data = m.pcolormesh(x, y, data, **pcolor_args)

    if ax_args:
        plt.setp(ax, **ax_args)

    ax.autoscale(enable=True, axis='both', tight=True)
    m.drawcoastlines(linewidth=.25, ax=ax)
    m.drawparallels(np.arange(60, 80, 10))
    m.drawmeridians(np.arange(0, 360, 30), labels=[1, 0, 0, 1], fontsize=12)
    x, y = m(90, 60)
    plt.text(x-25000, y+30000, s=u"60°N", fontsize=12, ha='right')
    x, y = m(90, 70)
    plt.text(x-25000, y+30000, s=u"70°N", fontsize=12, ha='right')

    text_kwargs = {
        'fontsize': 14,
        'transform': ax.transAxes,
    }

    cb = m.colorbar(mappable=graphed_data, location='right', extend='both', format='%.3g')
    cb.set_label(cblabel, fontsize=text_kwargs['fontsize'])
    cb.ax.tick_params(labelsize=text_kwargs['fontsize'])

    x_coord = -0.103  # Lines up just about right with the 90° label, in three_model at least
    if anom:
        ax.text(x_coord, .1, 'min: {:.2g}'.format(data.min()), text_kwargs)
        ax.text(x_coord, .05, 'max: {:.2g}'.format(data.max()), text_kwargs)
        ax.text(x_coord, 0, 'rmse: {:.2g}'.format(np.sqrt(np.mean(data ** 2))), text_kwargs)
    else:
        ax.text(x_coord, .1, 'min: {:.2g}'.format(data.min()), text_kwargs)
        ax.text(x_coord, .05, 'max: {:.2g}'.format(data.max()), text_kwargs)
        ax.text(x_coord, 0, 'mean: {:.2g}'.format(data.mean()), text_kwargs)

    if depth:
        ax.text(.7, 0, 'Depth: {:.2f}m'.format(depth), text_kwargs)


def map_comparison(lon, lat, data1, data2, cblabel='', depth='', kwargs=None):
    # Three panels to plots the obs, data and anomaly
    fig, (axl, axm, axr) = plt.subplots(3, 1, figsize=(8, 24))

    # compute the anomaly
    try:
        anom = data1 - data2
    except ValueError:
        print 'data1.shape != data2.shape, cannot make anomaly'
        raise

    data1_args = default_pcolor_args(data1)
    data2_args = default_pcolor_args(data2)
    anom_args = default_pcolor_args(anom, anom=True)

    # Check for pcolor and ax args. Assign default pcolor args for anomaly if
    # none have been provided, and for all others set to empty dict.
    manually_set = False
    for dargs, var in zip(['data1_args', 'data2_args', 'anom_args'], [data1_args, data2_args, anom_args]):
        if dargs in kwargs.keys():
            if 'pcolor_args' in kwargs[dargs].keys():
                for x in kwargs[dargs]['pcolor_args']:
                    var[x] = kwargs[dargs]['pcolor_args'][x]
                    manually_set = True
        if 'ax_args' not in kwargs[dargs].keys():
            kwargs[dargs]['ax_args'] = {}
    if not manually_set:
        data1_args['vmax'] = data2_args['vmax'] = max(data1_args['vmax'], data2_args['vmax'])
        data1_args['vmin'] = data2_args['vmin'] = min(data1_args['vmin'], data2_args['vmin'])

    if 'title' not in kwargs['data1_args']['ax_args'].keys():
        kwargs['data1_args']['ax_args']['title'] = 'Model'
    if 'title' not in kwargs['data2_args']['ax_args'].keys():
        kwargs['data2_args']['ax_args']['title'] = 'Observations'
    if 'title' not in kwargs['anom_args']['ax_args'].keys():
        kwargs['anom_args']['ax_args']['title'] = 'Anomaly'

    try:
        npolar_map(lon, lat, data1, ax=axl, ax_args=kwargs['data1_args']['ax_args'],
                   pcolor_args=data1_args, cblabel=cblabel, depth=depth)

        npolar_map(lon, lat, data2, ax=axm, ax_args=kwargs['data2_args']['ax_args'],
                   pcolor_args=data2_args, cblabel=cblabel, depth=depth)

        npolar_map(lon, lat, anom, ax=axr, ax_args=kwargs['anom_args']['ax_args'],
                   pcolor_args=anom_args, cblabel=cblabel, depth=depth, anom=True)
    except ValueError:
        raise


def proc_plots(plots, obs4comp):
    """Process a list of 'plots'

       Each 'plot' is a dict of key value pairs, defining:
       'ifile' : [str] specifying the input file to load data from
       'variables' : [list] of variables to plot, where each var is a [str]
       'plot_type' : 'section' or 'npolar_map'
       'pcolor_args' : [dict] option set of key-value pairs passed to pcolor
       'ax_args' : [dict] optional set of key-value pairs used to set axis attributes.

        Output pdfs are saved to plots/ and into a merged pdf called ArcticGraphs_plots.pdf
    """
    plots_out = []

    if os.path.isdir('./plots'):
        print("plots directory exists...overwriting")
        shutil.rmtree('./plots')
    else:
        print("Creating ./plots/ ...")

    os.mkdir('./plots')

    for i, p in enumerate(plots):
        os.mkdir('./plots/{}'.format(i))
        if p['plot_type'] == 'npolar_map_comp':
            for dargs in ['data1_args', 'data2_args', 'anom_args']:
                if 'kwargs' not in p.keys():
                    p['kwargs'] = {}
                if dargs not in p['kwargs'].keys():
                    p['kwargs'][dargs] = {}
                if 'ax_args' not in p['kwargs'][dargs]:
                    p['kwargs'][dargs]['ax_args'] = {}

        # Loop over each variable in the plot, and save a plot for it.
        for x, v in enumerate(p['variables']):
            data, units, lon, lat, depth, dimensions, years = load(p['ifile'], v)

            if 'pcolor_args' not in p.keys():
                p['pcolor_args'] = None
            elif not p['pcolor_args'] is None:
                p['pcolor_args']['rasterized'] = True

            p['cblabel'] = units

            # Check is ax_args exists, and if it has a title, if no title set to the varname
            if 'ax_args' not in p.keys():
                p['ax_args'] = {'title': v}
            else:
                p['ax_args']['title'] = v

            if p['plot_type'] == 'npolar_map':
                print 'plotting npolar map of', v

                if data.ndim > 2:
                    if 'plot_depth' not in p.keys():
                        p['plot_depth'] = np.round(depth.min())
                        print 'global_map: plot_depth not specified for {}, using {}'.format(v, p['plot_depth'])
                    else:
                        p['plot_depth'] = np.round(p['plot_depth'])

                    try:
                        depth_index = (np.abs(depth - p['plot_depth'])).argmin()  # Find the nearest depth to the
                        p['plot_depth'] = depth[depth_index]                      # one requested

                        if len(data.shape) == 4:
                            data = data[:, depth_index, :, :]
                        elif len(data.shape) == 3:
                            data = data[depth_index, :, :]

                    except:
                        print('Failed to extract depth {}  for {}'.format(p['plot_depth'], v))
                        raise

                if 'plot_depth' not in p.keys():
                    p['plot_depth'] = 0

                if len(data.shape) == 3:  # If there are multiple time steps,
                    for i in range(0, data.shape[0]):  # Plot each one
                        pcolor_args = default_pcolor_args(data)
                        try:
                            npolar_map(lon, lat, data[i, :, :], ax_args=p['ax_args'], pcolor_args=pcolor_args,
                                       cblabel=p['cblabel'], depth=p['plot_depth'])
                            plot_name = 'plots/{}/{}_{}_map_{}.pdf'.format(i, v, p['plot_depth'], years[i])

                            plots_out.append(plot_name)
                            plt.savefig(plot_name, bbox_inches='tight')
                        except ValueError:
                            print 'Graph number {} of {} has no data, skipping...'.format(i+1, v)

                else:
                    try:
                        npolar_map(lon, lat, data, ax_args=p['ax_args'], pcolor_args=p['pcolor_args'],
                                   cblabel=p['cblabel'], depth=p['plot_depth'])
                        plot_name = 'plots/{}/{}_{:.2f}_map.pdf'.format(i, v, p['plot_depth'])
                        plots_out.append(plot_name)
                        plt.savefig(plot_name, bbox_inches='tight')
                    except ValueError:
                        print 'The graph {} has no data, skipping...'.format(v)

            if (p['plot_type'] == 'npolar_map_comp') or (p['plot_type'] == 'taylor_plot'):
                try:
                    comp_var = p['compare_to'][x]
                except KeyError:
                    print 'No comparison variable provided for', v
                    raise
                try:
                    obs_data, _, _, _, obs_depth, dimensions, years = load(obs4comp[comp_var], comp_var)
                except ValueError:
                    print(v + ' not provided on obs4comp dict, cannot compare')
                    raise

                if p['plot_type'] != 'taylor_plot':
                    for dargs, atit in zip(
                            ['data1_args', 'data2_args', 'anom_args'],
                            ['Observations', 'NEMO', 'NEMO - Obs.']):
                        p['kwargs'][dargs]['ax_args']['title'] = ('{} ({})'.format(atit, v))
                else:
                    p['ax_args']['title'] = v

                if p['plot_type'] == 'npolar_map_comp':
                    print 'plotting north polar map comparison of {} with observations {}.'.format(v, comp_var)

                    if data.ndim > 2:
                        if 'plot_depth' not in p.keys():
                            p['plot_depth'] = np.round(depth.min())
                            print 'global_map: plot_depth not specified for {}, using {}'.format(v, p['plot_depth'])
                        else:
                            p['plot_depth'] = np.round(p['plot_depth'])

                        try:
                            if p['plot_depth'] in np.round(depth):
                                depth_ind = np.where(np.round(depth) == p['plot_depth'])[0][0]
                            else:
                                print 'Specified plot_depth of {}m not a depth in the NEMO file...looking for ' \
                                      'closest match instead...'.format(p['plot_depth'])
                                anom = np.round(depth) - p['plot_depth']
                                depth_ind = np.where(abs(anom) == abs(anom).min())[0][0]
                                p['plot_depth'] = np.round(depth[depth_ind])

                            data = data[depth_ind, :, :]
                        except:
                            print 'Failed to extract depth {} for {}'.format(p['plot_depth'], v)
                            raise

                    elif 'plot_depth' not in p.keys():
                        p['plot_depth'] = 0

                    print 'Plotting at {}m'.format(p['plot_depth'])

                    if obs_data.ndim == 3:
                        try:
                            ind_obs = np.where(np.round(obs_depth) == p['plot_depth'])[0][0]
                            obs_data = obs_data[ind_obs, :, :]
                        except:
                            print 'Failed to extract depth {} for observed {}'.format(p['plot_depth'], comp_var)
                            raise

                    try:
                        map_comparison(lon, lat, data, obs_data, cblabel=p['cblabel'], depth=p['plot_depth'], kwargs=p['kwargs'])
                        plot_name = 'plots/{}/{}_{}_map-comp.pdf'.format(i, v, p['plot_depth'])
                        plots_out.append(plot_name)
                        plt.savefig(plot_name, bbox_inches='tight')
                    except ValueError:
                        print 'The comparison of {} has no data in one o, skipping...'.format(v)

                if p['plot_type'] == 'taylor_plot':
                    taylor_plot(data, obs_data, ax_args=p['ax_args'])
                    plot_name = 'plots/{}/{}_taylor.pdf'.format(i, v)
                    plots_out.append(plot_name)
                    plt.savefig(plot_name, bbox_inches='tight')

        plt.close('all')
        # END of plot loop

    if len(plots_out):
        subprocess.Popen(('pdfunite ' + ' '.join(plots_out) + ' plots/ArcticGraphs_plots.pdf'), shell=True).wait()
    else:
        print '\nNo plots were generated.'
