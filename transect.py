import plot_tools as pt
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

CAA_TRANSECT = {
    'lon': [46, 69, 108, 139, 162, 180, 202, 221, 241, 254, 260, 310, 466],
    'lat': [205, 205, 177, 115, 90, 57, 59, 50, 83, 83, 105, 119, 45],
    'str': 'CAA_Transect',
    'ylim_max': 38
}

ARC_TRANSECT = {}

input_file = 'models/CanOE2_NAA-EPM032_365h_19840101_19841231_ptrc_T.nc'
var = ['no3']
transect_type = CAA_TRANSECT  # Must be one of CAA_TRANSECT, ARC_TRANSECT
file_type = ['canoe']  # Must be one of ['canoe', 'pisces', 'cmoc'] or []. Determines the scaling applied to the data

CANOE_SCALE_FACTORS = {
    'PHYC': 1,
    'ZOO': 1,
    'PHY2C': 1,
    'ZOO2': 1,
    'no3': 1,
    'votemper': 1,
    'vosaline': 1,
    'ileadfra': 1,
    'iicethic': 1,
    'dic': 1e6,
    'PO4': 1e6/122,
    'Si': 1e6,
    'dFe': 1e-3,
    'vovecrtz': 1,
    'vozocrtx': 1,
    'vomecrty': 1,
    'alk': 1e6,
    'PPPHY': 1e9,
    'PPPHY2': 1e9,
    'NCHL': 1,
    'DCHL': 1,
    'heup': 1,
    'PAR': 1,
    'PH': 1,
    'soshfldo': 1,
    'somxl010': 1,
    'DOC': 1e6,
    'O2': 1e6,
    'NH4': 1,
}

PISCES_CMOC_SCALE_FACTORS = {
    'PHY': 1e6,
    'ZOO': 1e6,
    'PHY2': 1e6,
    'ZOO2': 1e6,
    'no3': 1/7.625e-6,
    'votemper': 1,
    'vosaline': 1,
    'ileadfra': 1,
    'iicethic': 1,
    'dic': 1e6,
    'PO4': 1e6/122,
    'Si': 1e6,
    'Fer': 1e9,
    'vovecrtz': 1,
    'vozocrtx': 1,
    'vomecrty': 1,
    'alk': 1e6,
    'PPPHY': 1e9,
    'PPPHY2': 1e9,
    'NCHL': 1e6,
    'DCHL': 1e6,
    'heup': 1,
    'PAR': 1,
    'PH': 1,
    'soshfldo': 1,
    'somxl010': 1,
}

if os.path.isdir('./plots'):
    print("plots directory exists...overwriting")
    shutil.rmtree('./plots')
else:
    print("Creating ./plots/ ...")

os.mkdir('./plots')

for i in range(0, len(var)):
    data, units, lon, lat, depths, dimensions, years = pt.load(input_file, var[i])

    if file_type[0] == 'canoe':
        data *= CANOE_SCALE_FACTORS[var[i]]
        print 'CanOE scaling applied to data.'
    elif file_type[0] == ('pisces', 'cmoc'):
        data *= PISCES_CMOC_SCALE_FACTORS[var[i]]
        print 'PISCES/CMOC scaling applied to data.'
    else:
        print 'No scaling applied to data.'

    index_lon = np.array([], dtype=int)
    index_lat = np.array([], dtype=int)
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
        graphed_data = plt.pcolormesh(transect_data, cmap='viridis')
        color_bar = plt.colorbar(graphed_data, extend='both', format='%.3g')
        color_bar.set_label(units)
        plt.xticks([0, 23, 71, 140, 210, 253, 291, 304, 326, 377, 549],
                   ['H', 'I', 'J', 'M', 'N', 'O', 'P', 'Q', 'R', 'T(OSD)', 'W/OSB'])
        plt.xlim([0, 549])
        plt.ylim([transect_type['ylim_max'], 0])
        plt.ylabel('Depth (m)')
        yticks = ax.get_yticks()
        for k in range(0, len(yticks)):
            yticks[k] = '{:.2f}'.format(depths[int(yticks[k])])
        ax.set_yticklabels(yticks, va='top')
        plot_name = 'plots/{}_{}_{}.pdf'.format(var[i], transect_type['str'], years[j])
        plt.savefig(plot_name, bbox_inches='tight')
