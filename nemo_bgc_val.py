"""Produce a set of plots to validate the NEMO BGC


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
import subprocess
import nemo_plot_tools as npt

def bgc_main(plot_vars):

    for ifile, variables in plot_vars.iteritems():
        for var, kwargs in variables.iteritems():
            print var
            data, lon, lat, depth = npt.load(ifile, var)
            npt.global_map(lon, lat, data[0,:,:], title=var, **kwargs)
            plt.savefig( var + '_map.pdf', bbox_inches='tight')

            npt.section(lat, depth, data.mean(axis=2), ax_args={'xlim':[-80, 80]})
            plt.savefig( var + '_section.pdf', bbox_inches='tight')

    subprocess.Popen('pdfunite *_map.pdf output_maps.pdf', shell=True).wait()


if __name__ == '__main__':
    vars_ptrc = ['NO3', 'NCHL', 'DIC', 'PHY', 'ZOO']
    ifile_ptrc = ('/raid/ra40/data/ncs/nemo_out/nue/' +
                  'mc_nue_1m_00010101_00011231_ptrc_t.nc.001')



    plot_vars = {ifile_ptrc : {'NO3' : {},
                               'NCHL' : {},
                               'DIC': {'vmin' : 1800, 'vmax' : 2300},
                               'PHY' : {},
                               'ZOO' : {}
                               }
                 }
    bgc_main(plot_vars)