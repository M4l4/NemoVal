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
import argparse
import os

import nemo_plot_tools as npt

# Get input files from the command line:
parser = argparse.ArgumentParser(prog='nemo_bgc_val.py', usage='%(prog)s [options]')
parser.add_argument('infiles', nargs=1)
args = parser.parse_args()

ifile = args.infiles[0]
ifile = 'NAA-EPM032_337h_20060101_20251231_ptrc_T.nc'

if not os.path.exists(ifile):
    raise SystemExit('File not found:' + ifile)

#===================================================================================
#      DEFINE OBSERVATIONS
#
#      obs4comp is a dict where the keys are (NEMO & observation file) variable
#      names, and the values are the full path to the observation file.
#===================================================================================

obs_root = '/home/mal/Artic/ftp.cccma.ec.gc.ca/pub/nswart/obs4comp/'

obs4comp = {  # 'NO3': obs_root + 'uncs_orca2_data_data_n_an_nomask.nc',
    # 'DIC': obs_root + 'uncs_orca2_data_data_TCO2_nomask.nc',
    # 'NCHL': obs_root + 'uncs_seawifs_mean_1998_2005.nc',
    'O2': obs_root + 'CMOCNAA-EPM032_1460h_20020101_20030101_ptrc_.nc',
    # 'Alkalini': obs_root + 'uncs_orca2_data_data_Alk_nomask.nc',
    # 'Cflx': obs_root + 'uncs_cflx_park_takahasi_1990-2009.nc',
    # 'EPC100': obs_root + 'mol-m2-s_AWI_export.nc',
    # 'salt': obs_root + 'phc3.0_monthly.nc',
    # 'temp': obs_root + 'phc3.0_monthly.nc',
    'sosaline': 'NAA-EPM032_337h_20060101_20251231_grid_T.nc',
}

#===================================================================================
#      DEFINE PLOT LIST
#===================================================================================
#
#       Each 'plot' is a dict of key value pairs, defining:
#
#       'ifile' : [str] specifying the input file to load data from
#
#       'variables' : [list] of variables to plot, where each var is a [str]
#
#       'plot_type' : 'section' or 'global_map' or 'global_map_comp' or
#                     'section_comp'
#
#       'pcolor_args' : [dict] option set of key-value pairs passed to pcolor
#
#       'ax_args' : [dict] optional set of key-value pairs used to set axis
#                   attributes.
#
#       'kwargs' : [dict] for 'comp' plots, a dict containing pcolor_args and
#                   ax_args dicts, for one or each of the keys 'data1_args',
#                   'data2_args', 'anom_args'.
#===================================================================================
plots = [
    {'ifile': ifile,
     'variables': ['PHY'],
     'plot_depth': 0,
     'plot_type': 'npolar_map',
     },
    # {'ifile': ifile,
    #  'variables': ['isssalin'],
    #  'compare_to': ['sosaline'],
    #  'plot_depth': 0,
    #  'plot_type': 'npolar_map_comp',
    #  },
    # {'ifile': ifile,
    #  'variables': ['isssalin'],
    #  'compare_to': ['sosaline'],
    #  'plot_depth': 0,
    #  'plot_type': 'taylor_plot',
    #  },

    # {'ifile': ifile_ptrc,
    #  'variables': ['NO3', 'O2'],
    #  'plot_type': 'global_map_comp',
    #  'plot_depth': 5,
    #  },
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['NO3'],
    #  'plot_type': 'global_map_comp',
    #  'plot_depth': 5,
    #  },
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['NO3', 'O2'],
    #  'plot_type': 'global_map_comp',
    #  'plot_depth': 2000,
    #  },
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['DIC'],
    #  'plot_type': 'global_map_comp',
    #  #'kwargs' : {'data1_args' : {'pcolor_args' :{'vmin' : 1800, 'vmax': 2200}},
    #  #'data2_args' : {'pcolor_args' :{'vmin' : 1800, 'vmax': 2200}},
    #  #}
    #  },
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['DIC'],
    #  'plot_type': 'global_map_comp',
    #  #'kwargs' : {'data1_args' : {'pcolor_args' :{'vmin' : 1800, 'vmax': 2200}},
    #  #'data2_args' : {'pcolor_args' :{'vmin' : 1800, 'vmax': 2200}},
    #  #}
    #  'plot_depth': 2000,
    #  },
    #
    # #{'ifile' : ifile_ptrc,
    # #'variables' : ['Alkalini'],
    # #'plot_type' : 'global_map_comp',
    # #},
    #
    #
    # {'ifile': ifile_diad,
    #  'variables': ['Cflx'],
    #  'plot_type': 'global_map_comp',
    #  'kwargs': {'data1_args': {'pcolor_args': {'vmin': -2e-7, 'vmax': 2e-7, 'cmap': npt.anom_cmap()}},
    #             'data2_args': {'pcolor_args': {'vmin': -2e-7, 'vmax': 2e-7, 'cmap': npt.anom_cmap()}},
    #             'anom_args': {'pcolor_args': {'vmin': -2e-7, 'vmax': 2e-7}}
    #             }
    #  },
    #
    # {'ifile': ifile_diad,
    #  'variables': ['EPC100'],
    #  'plot_type': 'global_map_comp',
    #  'kwargs': {'data1_args': {'pcolor_args': {'vmin': 0, 'vmax': 2e-7}},
    #             'data2_args': {'pcolor_args': {'vmin': 0, 'vmax': 2e-7}},
    #             }
    #  },
    #
    # {'ifile': ifile_diad,
    #  'variables': ['Oflx', 'Nfix', 'PPPHY', 'EPC100'],
    #  'plot_type': 'global_map',
    #  },
    #
    # #{'ifile' : ifile_diad,
    # #'variables' : ['LNFe', 'LNnut'],
    # #'plot_type' : 'global_map',
    # #},
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['PHY', 'ZOO'],
    #  'plot_type': 'global_map',
    #  },
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['NCHL'],
    #  'plot_type': 'global_map_comp',
    #  'kwargs': {'data1_args': {'pcolor_args': {'vmin': 0, 'vmax': 0.5}},
    #             'data2_args': {'pcolor_args': {'vmin': 0, 'vmax': 0.5}},
    #             'anom_args': {'pcolor_args': {'vmin': -0.25, 'vmax': 0.25}}
    #             }
    #  },
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['NO3', 'O2', 'DIC'],  #, 'Alkalini'],
    #  'plot_type': 'section_comp',
    #  },
    #
    # {'ifile': ifile_ptrc,
    #  'variables': ['PHY', 'ZOO'],
    #  'plot_type': 'section',
    #  'ax_args': {'ylim': [200, 0]}
    #  },
]

npt.proc_plots(plots, obs4comp)
