import plot_tools as pt
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

CAA_TRANSECT = {
    'lon': [46,  69,  108, 139, 162, 168, 173, 176, 180, 192, 201, 210, 221, 225, 241, 254, 260, 310, 466],
    'lat': [205, 205, 177, 115, 90,  72,  67,  60,  56,  56,  59,  50,  50,  52,  82,  83,  105, 119, 45],
    'ylim_max': 38,
    'title': 'CAA Transect',
    'ftype': 'CAA_Transect',
}

# Transect on 140W, from 70N to 80N
SHELF_1_TRANSECT = {
    'lon': [113, 181],
    'lat': [127, 184],
    'ylim_max': 38,
    'title': 'Shelf Transect 1',
    'ftype': 'Shelf_1_Transect',
}

# Transect on 150W, from 70N to 80N
SHELF_2_TRANSECT = {
    'lon': [97, 173],
    'lat': [153, 196],
    'ylim_max': 38,
    'title': 'Shelf Transect 1',
    'ftype': 'Shelf_2_Transect',
}

input_file = 'models/CanOE2_NAA-EPM032_365h_19840101_19841231_ptrc_T.nc'
var = ['no3']
transect_type = CAA_TRANSECT  # One of CAA_TRANSECT, SHELF_1_TRANSECT, SHELF_2_TRANSECT
file_type = 'canoe'  # Must be one of 'canoe', 'pisces', 'cmoc', or None. Determines the scaling applied to the data

# Number of discrete colors in the color bar, default None
color_bar_steps = 5

if os.path.isdir('./plots'):
    print("plots directory exists...overwriting")
    shutil.rmtree('./plots')
else:
    print("Creating ./plots/ ...")

os.mkdir('./plots')

for i in range(0, len(var)):
    data, units, lon, lat, depths, dimensions, years = pt.load(input_file, var[i], file_type)

    index_lon = np.array([], dtype=int)
    index_lat = np.array([], dtype=int)

    pcolor_args = pt.default_pcolor_args(data, color_bar_steps)
    for j in range(0, len(transect_type['lon'])-1):
        length = int(np.hypot(transect_type['lon'][j + 1] - transect_type['lon'][j],
                              transect_type['lat'][j + 1] - transect_type['lat'][j]))
        index_x = np.linspace(transect_type['lon'][j], transect_type['lon'][j + 1], length, dtype=int)
        index_y = np.linspace(transect_type['lat'][j], transect_type['lat'][j + 1], length, dtype=int)
        index_lon = np.append(index_lon, index_x)
        index_lat = np.append(index_lat, index_y)

    for j in range(0, len(years)):
        plt.close()
        fig, ax = plt.subplots()
        if len(data.shape) == 4:
            transect_data = data[j, :, index_lat, index_lon]
        else:
            transect_data = data[:, index_lat, index_lon]

        graphed_data = plt.pcolormesh(transect_data, **pcolor_args)

        plt.title('{} {}'.format(var[i], transect_type['title']))

        color_bar = plt.colorbar(graphed_data, extend='both', format='%.3g')
        color_bar.set_label(r'${}$'.format(units))

        if transect_type == CAA_TRANSECT:
            plt.xlim([0, 550])
            plt.xticks([0, 23, 71, 140, 210, 253, 291, 304, 326, 377, 550],
                       ['H', 'I', 'J', 'M', 'N', 'O', 'P', 'Q', 'R', 'T(OSD)', 'W/OSB'])
        elif transect_type == SHELF_1_TRANSECT:
            plt.xlim([0, 87])
            plt.xticks([0, 45, 87], ['70N', '75N', '80N'])
        else:
            plt.xlim([0, 86])
            plt.xticks([0, 44, 86], ['70N', '75N', '80N'])

        plt.ylim([transect_type['ylim_max'], 0])
        plt.ylabel('Depth (m)')
        yticks = ax.get_yticks()
        ax.set_yticklabels([int(depths[int(yticks[k])]) for k in range(0, len(yticks))], va='top')

        plot_name = 'plots/{}_{}_{}.pdf'.format(var[i], transect_type['ftype'], years[j])
        plt.savefig(plot_name, bbox_inches='tight')
