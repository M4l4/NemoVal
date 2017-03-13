import nemo_plot_tools as npt
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

f1 = 'CanOE2_NAA-EPM032_365h_19840101_19841231_ptrc_T.nc'
var1 = ['dFe', 'no3']
f2 = 'PISCESS_NAA-EPM032_365h_19840101_19841231_ptrc_T.nc'
var2 = ['Fer', 'no3']
depth = 10

if os.path.isdir('./plots'):
    print("plots directory exists...overwriting")
    shutil.rmtree('./plots')
else:
    print("Creating ./plots/ ...")

os.mkdir('./plots')

for i in range(0, len(var1)):
    data1, units1, lon, lat, depths, dimensions, years = npt.load(f1, var1[i])
    data2, units2, _, _, _, _, _ = npt.load(f2, var2[i])

    pargs1 = npt.default_pcolor_args(data1)
    pargs2 = npt.default_pcolor_args(data2)

    vmin = min(pargs1['vmin'], pargs2['vmin'])
    vmax = max(pargs1['vmax'], pargs2['vmax'])
    pargs1['vmin'] = vmin
    pargs1['vmax'] = vmax

    depth_index = (np.abs(depths-depth)).argmin()
    plot_depth = depths[depth_index]
    print 'ploting at {}m'.format(plot_depth)

    if len(data1.shape) == 4:
        data1 = data1[:, depth_index, :, :]
        data2 = data2[:, depth_index, :, :]
    else:
        data1 = data1[depth_index, :, :]
        data2 = data2[depth_index, :, :]

    print 'Plotting graph of {}...'.format(var1[i])
    if len(data1.shape) == 3:  # If there are multiple time steps,
        for j in range(0, data1.shape[0]):  # Plot each one
            print 'Plotting graph {} of {}...'.format(j + 1, data1.shape[0])

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 16))
            plt.suptitle('\n'.join([var1[i], years[j]]), fontsize=30)
            ax_args1 = {'title': f1.split('_')[0],
                        }
            ax_args2 = {'title': f2.split('_')[0],
                        }
            try:
                npt.npolar_map(lon, lat, data1[j, :, :], ax=ax1, ax_args=ax_args1, pcolor_args=pargs1, cblabel=units1,
                               depth=plot_depth)
                npt.npolar_map(lon, lat, data2[j, :, :], ax=ax2, ax_args=ax_args2, pcolor_args=pargs1, cblabel=units1,
                               depth=plot_depth)

                plot_name = 'plots/{}_{}_map_{}.pdf'.format(var1[i], plot_depth, j + 1)
                plt.savefig(plot_name, bbox_inches='tight')
            except ValueError:
                print 'Graph number {} of {} has no data, skipping...'.format(j + 1, 'ah why')
            plt.close(fig)
    else:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 16))
        plt.suptitle('\n'.join([var1[i], years[0]]), fontsize=30)
        ax_args1 = {'title': f1.split('_')[0],
                    }
        ax_args2 = {'title': f2.split('_')[0],
                    }
        try:
            npt.npolar_map(lon, lat, data1, ax=ax1, ax_args=ax_args1, pcolor_args=pargs1, cblabel=units1,
                           depth=plot_depth)
            npt.npolar_map(lon, lat, data2, ax=ax2, ax_args=ax_args2, pcolor_args=pargs1, cblabel=units1,
                           depth=plot_depth)

            plot_name = 'plots/{}_{}_map.pdf'.format(var1[i], plot_depth)
            plt.savefig(plot_name, bbox_inches='tight')
        except ValueError:
            print 'Graph number of {} has no data, skipping...'.format(var1[i])
        plt.close(fig)
