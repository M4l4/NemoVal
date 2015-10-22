"""Produce a set of plots to validate the NEMO BGC


neil.swart@ec.gc.ca, 10/2015
"""
import nemo_plot_tools as npt

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
#       'plot_type' : 'section' or 'global_map'
#
#       'pcolor_args' : [dict] option set of key-value pairs passed to pcolor
#
#       'ax_args' : [dict] optional set of key-value pairs used to set axis
#                   attributes.
#===================================================================================

# The ptrc input file
ifile_ptrc = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                  'mc_nue_1m_20000101_20001231_ptrc_t.nc.001')

ifile_diad = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                  'mc_nue_1m_20000101_20001231_diad_t.nc.001')

# The list of plots
plots = [
         {'ifile' : ifile_ptrc,
          'variables' : ['NO3', 'O2'],
          'plot_type' : 'global_map_comp',
         },

         {'ifile' : ifile_ptrc,
          'variables' : ['DIC', 'Alkalini'],
          'plot_type' : 'global_map_comp',
          'kwargs' : {'data1_args' : {'pcolor_args' :{'vmin' : 1800, 'vmax': 2200}},
                      'data2_args' : {'pcolor_args' :{'vmin' : 1800, 'vmax': 2200}},
                      }
         },

         {'ifile' : ifile_diad,
          'variables' : ['Cflx'],
          'plot_type' : 'global_map_comp',
         },

         {'ifile' : ifile_ptrc,
          'variables' : ['PHY', 'ZOO'],
          'plot_type' : 'global_map',
         },

         {'ifile' : ifile_ptrc,
          'variables' : ['NCHL'],
          'plot_type' : 'global_map_comp',
          'kwargs' : {'data1_args' : {'pcolor_args' :{'vmin' : 0, 'vmax': 0.5}},
                      'data2_args' : {'pcolor_args' :{'vmin' : 0, 'vmax': 0.5}},
                      'anom_args'  : {'pcolor_args' :{'vmin' : -0.25, 'vmax': 0.25}}
                      }
         },

         {'ifile' : ifile_ptrc,
          'variables' : ['NO3', 'O2', 'DIC', 'Alkalini'],
          'plot_type' : 'section_comp',
         },


         {'ifile' : ifile_ptrc,
          'variables' : ['PHY', 'ZOO'],
          'plot_type' : 'section',
          'ax_args' : {'ylim' : [200, 0]}
         },

        ]

#===================================================================================
#      DEFINE OBSERVATIONS
#
#      obs4comp is a dict where the keys are (NEMO & observation file) variable
#      names, and the values are the full path to the observation file.
#===================================================================================

obs_root = '/raid/ra40/data/ncs/nemo_out/obs4comp/'

obs4comp = {'NO3' : obs_root + 'uncs_orca2_data_data_n_an_nomask.nc',
            'DIC' : obs_root + 'uncs_orca2_data_data_TCO2_nomask.nc',
            'NCHL': obs_root + 'uncs_seawifs_mean_1998_2005.nc',
            'O2'  : obs_root + 'umol-l_uncs_orca2_data_data_o_an_nomask.nc',
            'Alkalini' : obs_root + 'uncs_orca2_data_data_Alk_nomask.nc',
            'Cflx' : obs_root + 'adj_uncs_park-ccmp-takahashi_fgco2_1990-2009.nc'
           }

npt.proc_plots(plots, obs4comp)