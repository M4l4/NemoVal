# coding=utf-8
""" PlotTools: A Basic library for quick and dirty plots of model results.

version 0.2

malachitehmiller@gmail.com, 03/2017
"""
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
import json

import taylor

plt.close('all')
font = {'size': 12}
plt.rc('font', **font)
mpl.rcParams['image.cmap'] = 'viridis'
mpl.rcParams['mathtext.default'] = 'regular'


def default_pcolor_args(data, steps=None, anom=False):
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
        cmap = cm.get_cmap("viridis", steps)

    d = {'vmin': vmin,
         'vmax': vmax,
         'cmap': cmap,
         'rasterized': True
         }

    return d


def load(ifile_fp, varname, file_type=None):
    path, ifile = os.path.split(ifile_fp)

    if not os.path.isfile('remapped/yearmean_{}'.format(ifile)):
        try:
            import cdo
            cdo = cdo.Cdo()
            cdo.env = {'REMAP EXTRAPOLATE': 'off'}
        except ImportError:
            raise SystemExit('CDO must be installed to preproces files.')
        cdo.yearmean(input=ifile_fp, output='remapped/yearmean_' + ifile)
        with Dataset('parameters/masks.nc', 'r') as mask, Dataset('remapped/yearmean_' + ifile, 'a') as to_mask:
            if file_type:
                with open('parameters/scale_factors.json') as json_file:
                    var_defs = json.load(json_file)
                    for var in to_mask.variables:
                        if var in var_defs[file_type].keys():
                            to_mask[var].setncattr('scale', var_defs[file_type][var][0])
                            to_mask[var].units = var_defs[file_type][var][1]
                            to_mask[var].long_name = var_defs[file_type][var][2]
            for var in to_mask.variables:
                if len(to_mask[var].shape) == 4:
                    to_mask[var][:] = np.ma.masked_invalid(np.where(mask['tmask'][0, ], to_mask[var][:], np.NaN))

    with Dataset('remapped/yearmean_' + ifile, 'r') as nc:

        ncvar = nc.variables[varname]
        data = ncvar[:].squeeze()

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
            raise SystemExit('\nDepths not given for {} in {}.'.format(varname, ifile))

        try:
            dates = num2date(nc['time_counter'][:], units=nc['time_counter'].units,
                             calendar=nc['time_counter'].calendar)
            years = [x.year for x in dates]
        except IndexError:
            years = [0]

        lon = nc.variables['nav_lon']
        lon = lon[:].squeeze()
        lat = nc.variables['nav_lat']
        lat = lat[:].squeeze()

    return data, units, lon, lat, depth, dimensions, years


def load_remap(ifile_fp, varname, file_type=None):
    path, ifile = os.path.split(ifile_fp)

    if not os.path.isfile('remapped/intlev_{}'.format(ifile)):
        try:
            import cdo
            cdo = cdo.Cdo()
            cdo.env = {'REMAP EXTRAPOLATE': 'off'}
        except ImportError:
            print 'CDO must be installed to preproces files.'
            raise

        print 'preprocessing', ifile

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
                    if file_type:
                        with open('parameters/scale_factors.json') as json_file:
                            var_defs = json.load(json_file)
                            for var in to_mask.variables:
                                if var in var_defs[file_type].keys():
                                    print var
                                    to_mask[var].setncattr('scale', var_defs[file_type][var][0])
                                    to_mask[var].units = var_defs[file_type][var][1]
                                    to_mask[var].long_name = var_defs[file_type][var][2]
                    for var in to_mask.variables:
                        if len(to_mask[var].shape) == 4:
                            to_mask[var][:] = np.ma.masked_invalid(np.where(mask['tmask'][0, ], to_mask[var][:], np.NaN))
            cdo.intlevel(default_depths, input='-remapdis,parameters/default_grid temp_' + ifile,
                         output='remapped/intlev_' + ifile, options='-L', force=False)
            os.remove('temp_' + ifile)
        except (AttributeError, IndexError):
            print 'Data from {} is not on the NAA grid, no land mask applied.'.format(ifile)
            # Year average ifile, and interpolate to 1/2x1/2 degrees
            cdo.intlevel(default_depths, input='-remapdis,parameters/default_grid -yearmean ' + ifile_fp,
                         output='remapped/intlev_' + ifile, options='-L', force=False)
        print ifile, 'preprocessed.'

    remaped_file = 'remapped/intlev_' + ifile
    with Dataset(remaped_file, 'r') as nc:
        ncvar = nc.variables[varname]
        data = ncvar[:].squeeze()

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
            raise SystemExit('\nDepths not given for {} in {}.'.format(varname, ifile))

        try:
            dates = num2date(nc['time_counter'][:], units=nc['time_counter'].units,
                             calendar=nc['time_counter'].calendar)
            years = [x.year for x in dates]
        except IndexError:
            years = [0]

    lon = np.linspace(0, 360, 721)
    lat = np.linspace(40, 90, 101)

    return data, units, lon, lat, depth, dimensions, years


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


def npolar_map(lon, lat, data, ax_args=None, pcolor_args=None, cblabel='', depth=None, ax=None, anom=False, remap=True):
    """Pcolor a var in a polar map, using ax if supplied"""
    # setup a basic polar map

    if data.max() is ma.masked:
        raise ValueError

    if not ax:
        _, ax = plt.subplots(1, 1, figsize=(8, 8))

    if not pcolor_args:
        pcolor_args = default_pcolor_args(data, anom=anom)

    m = Basemap(projection='npstere', boundinglat=55, lon_0=0, resolution='l', round=True, ax=ax)

    if remap:
        lon_mesh, lat_mesh = np.meshgrid(lon, lat)
        x, y = m(lon_mesh, lat_mesh)
    else:
        x, y = m(lon, lat)
    graphed_data = m.pcolormesh(x, y, data, **pcolor_args)

    if ax_args:
        plt.setp(ax, **ax_args)

    ax.autoscale(enable=True, axis='both', tight=True)
    m.drawcoastlines(linewidth=.25, ax=ax)
    m.drawparallels(np.arange(60, 80, 10))
    m.drawmeridians(np.arange(0, 360, 30), labels=[1, 0, 0, 1], fontsize=12)
    x, y = m(90, 60)
    plt.text(x - 25000, y + 30000, s=u"60°N", fontsize=12, ha='right')
    x, y = m(90, 70)
    plt.text(x - 25000, y + 30000, s=u"70°N", fontsize=12, ha='right')

    text_kwargs = {
        'fontsize': 14,
        'transform': ax.transAxes,
    }

    cb = m.colorbar(mappable=graphed_data, location='right', extend='both', format='%.3g')
    cb.set_label(r'${}$'.format(cblabel), fontsize=text_kwargs['fontsize'])
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
        npolar_map(lon, lat, data1, kwargs['data1_args']['ax_args'], data1_args, cblabel, depth, axl)
        npolar_map(lon, lat, data2, kwargs['data2_args']['ax_args'], data2_args, cblabel, depth, axm)
        npolar_map(lon, lat, anom, kwargs['anom_args']['ax_args'], anom_args, cblabel, depth, axr, anom=True)
    except ValueError:
        raise


def proc_plots(plots, obs4comp):
    """Process a list of 'plots'

    See single_plots.py for a breakdown of the plot dictionary's keys.
    
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

        # Loop over each variable in the plot, and save a plot for it.
        for x, v in enumerate(p['variables']):
            try:
                if p['plot_type'] == 'npolar_map' and not p['remap']:
                    if 'file_type' not in p.keys():
                        p['file_type'] = None
                    data, units, lon, lat, depth, dimensions, years = load(p['ifile'], v, p['file_type'])
                else:
                    data, units, lon, lat, depth, dimensions, years = load_remap(p['ifile'], v)
            except KeyError:
                raise
                data, units, lon, lat, depth, dimensions, years = load_remap(p['ifile'], v)

            # Check is ax_args exists, and if it has a title, if no title set to the varname
            if 'ax_args' not in p.keys():
                p['ax_args'] = {'title': v}
            else:
                p['ax_args']['title'] = v

            try:
                depth_index = np.abs(depth - p['plot_depth']).argmin()  # Find the nearest depth
            except KeyError:
                raise SystemExit('Please specify a depth for the plot of {}.'.format(v))
            plot_depth = depth[depth_index]
            print 'plotting {} at {:.3g}m'.format(v, plot_depth)

            if p['plot_type'] == 'npolar_map':
                print 'plotting npolar map of', v

                try:
                    steps = p['color_bar_steps']
                except KeyError:
                    steps = None
                color_args = default_pcolor_args(data, steps)
                try:
                    remap = p['remap']
                except KeyError:
                    remap = True
                try:
                    color_args['vmax'] = p['pcolor_args']['vmax']
                    color_args['vmin'] = p['pcolor_args']['vmin']
                except (KeyError, TypeError):
                    pass

                try:
                    if len(data.shape) == 4:
                        data = data[:, depth_index, :, :]
                    elif len(data.shape) == 3:
                        data = data[depth_index, :, :]
                except:
                    print('Failed to extract depth {}  for {}'.format(p['plot_depth'], v))
                    raise

                if len(data.shape) == 3:  # If there are multiple time steps,
                    for j in range(0, data.shape[0]):  # Plot each one
                        try:
                            npolar_map(lon, lat, data[j, :, :], p['ax_args'], color_args, units, plot_depth,
                                       remap=remap)
                            plot_name = 'plots/{}/{}_{:.2f}_NPole_{}.pdf'.format(j, v, plot_depth, years[j])
                            plots_out.append(plot_name)
                            plt.savefig(plot_name, bbox_inches='tight')
                        except ValueError:
                            print 'Graph number {} of {} has no data, skipping...'.format(j + 1, v)

                else:
                    try:
                        npolar_map(lon, lat, data, p['ax_args'], color_args, units, plot_depth, remap=remap)
                        plot_name = 'plots/{}/{}_{:.2f}_map.pdf'.format(i, v, plot_depth)
                        plots_out.append(plot_name)
                        plt.savefig(plot_name, bbox_inches='tight')
                    except ValueError:
                        print 'The graph {} has no data, skipping...'.format(v)

            else:
                try:
                    comp_var = p['compare_to'][x]
                except KeyError:
                    print 'No comparison variable provided for', v
                    raise
                try:
                    comp_data, _, _, _, _, _, _ = load_remap(obs4comp[comp_var], comp_var)
                except ValueError:
                    print(v + ' not provided in obs4comp dict, cannot compare')
                    raise

                if p['plot_type'] == 'npolar_map_comp':
                    if 'kwargs' not in p.keys():
                        p['kwargs'] = None

                    if data.ndim > 2:
                        try:
                            if len(data.shape) == 4:
                                data = data[:, depth_index, :, :]
                            elif len(data.shape) == 3:
                                data = data[depth_index, :, :]
                        except:
                            print('Failed to extract depth {}  for {}'.format(p['plot_depth'], v))
                            raise

                    if comp_data.ndim == 3:
                        try:
                            comp_data = comp_data[depth_index, :, :]
                        except:
                            print 'Failed to extract depth {} for observed {}'.format(p['plot_depth'], comp_var)
                            raise

                    try:
                        map_comparison(lon, lat, data, comp_data, units, plot_depth, p['kwargs'])
                        plot_name = 'plots/{}/{}_{}_map-comp.pdf'.format(i, v, p['plot_depth'])
                        plots_out.append(plot_name)
                        plt.savefig(plot_name, bbox_inches='tight')
                    except ValueError:
                        print 'The comparison of {} has no data in one o, skipping...'.format(v)

                if p['plot_type'] == 'taylor_plot':
                    taylor_plot(data, comp_data, ax_args=p['ax_args'])
                    plot_name = 'plots/{}/{}_taylor.pdf'.format(i, v)
                    plots_out.append(plot_name)
                    plt.savefig(plot_name, bbox_inches='tight')

        plt.close('all')
        # END of plot loop

    if len(plots_out):
        subprocess.Popen(('pdfunite ' + ' '.join(plots_out) + ' plots/ArcticGraphs_plots.pdf'), shell=True).wait()
    else:
        print '\nNo plots were generated.'
