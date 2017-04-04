"""NemoView: Produce a set of plots to validate the NEMO BGC

version 0.1 / standard plots

Usage: 

    $ ipython  nemo_bgc_val.py

To define custom plots, in the upper section modify the "plots" list. To define
custom observational data to compare (must match CMIP5 variable names and units),
insert/change values in the 'obs4comp' dictionary in the lower section.

WARNING: **This code is under development. There _ will_ be bugs, and future
versions
           may not be backwards compatible**


neil.swart@canada.ca, 10/2015
"""
import os
import plot_tools as npt

ifile = 'models/CanOE2_NAA-EPM032_365h_19840101_19841231_ptrc_T.nc'

if not os.path.exists(ifile):
    raise SystemExit('File not found:' + ifile)

#===================================================================================
#      DEFINE OBSERVATIONS
#
#      obs4comp is a dict where the keys are (NEMO & observation file) variable
#      names, and the values are the full path to the observation file.
#===================================================================================

obs4comp = {
    # 'alk': 'models/PISCES_NAA-EPM032_365h_19840101_19841231_ptrc_T.nc'
}

#===================================================================================
#      DEFINE PLOT LIST
#===================================================================================
#
#       Each 'plot' is a dict of key value pairs, defining:
#
#       'ifile': str specifying the input file to load data from
#
#       'plot_type': 'npolar_map' or 'npolar_map_comp' or 'taylor_plot' (soon)
#
#       'remap': if ifile is on the NAA grid and 'plot_type' is 'npolar_map', set this to False to use the NAA grid
#
#       'depth': int specifying the depth to plot at
#
#       'variables': [str] of variables to plot
#
#       'compare_to': [str] of variables to make an anomaly with each item in 'variables'. There must be an equal number
#                     of items in 'variables' and 'compare_to' if 'plot_type' is 'npolar_map_comp' or 'taylor_plot'.
#
#       'color_bar_steps': [int] Number of discrete colors in the color bar
#
#       'pcolor_args': [dict] optional set of key-value pairs passed to pcolor
#
#       'comp_args': [dict] same as 'pcolor_args' but for the comparison data
#
#       'anom_args': [dict] same as 'pcolor_args' but for the anomaly
#
#       'ax_args': [dict] optional set of key-value pairs used to set axis attributes.
#
#===================================================================================
plots = [
    {'ifile': ifile,
     'plot_type': 'npolar_map',
     'file_type': 'canoe',
     'plot_depth': 0,
     'variables': ['PHYC', ],
     'remap': False,
     'color_bar_steps': 5,
     },
    # {'ifile': ifile,
    #  'plot_type': 'npolar_map',
    #  'plot_depth': 10,
    #  'variables': ['PHYC'],
    #  'pcolor_args': {'vmin': 0.3, 'vmax': .9}
    #  },
    # {'ifile': ifile,
    #  'plot_type': 'npolar_map_comp',
    #  'plot_depth': 0,
    #  'variables': ['alk'],
    #  'compare_to': ['alk'],
    #  },
    # {'ifile': ifile,
    #  'plot_type': 'npolar_map_comp',
    #  'variables': ['alk'],
    #  'compare_to': ['alk'],
    #  'pcolor_args': {'vmin': 0e-3, 'vmax': 2.4e-3},
    #  'comp_args': {'vmin': 1e-3, 'vmax': 2e-3},
    #  'anom_args': {'vmin': -5e-6, 'vmax': 5e-6}
    #  },
    # {'ifile': ifile,
    #  'variables': ['alk'],
    #  'compare_to': ['alk'],
    #  'plot_depth': 0,
    #  'plot_type': 'taylor_plot',
    #  },
]

npt.proc_plots(plots, obs4comp)
