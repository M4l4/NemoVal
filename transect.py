import plot_tools as pt
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

CAA_TRANSECT = {
    'lon': [240.00, 240.60, 241.22, 241.84, 242.47, 243.10, 243.75, 244.39, 245.05, 245.71, 246.37, 247.04, 247.72,
            248.40, 249.09, 249.78, 250.48, 251.19, 251.90, 252.61, 253.33, 254.05, 254.78, 255.51, 256.24, 256.98,
            257.72, 258.47, 259.22, 259.97, 260.72, 261.48, 262.23, 262.99, 263.75, 264.51, 265.28, 266.04, 266.80,
            267.57, 268.33, 269.09, 269.85, 270.61, 271.37, 272.13, 272.89, 273.64, 274.40, 275.14, 275.89, 276.64,
            277.38, 278.11, 278.85, 279.58, 280.30, 281.02, 281.74, 282.45, 283.16, 283.86, 284.56, 285.25, 285.93,
            286.61, 287.29, 287.96, 288.62, 289.28, 289.93, 290.57, 291.21, 291.84, 292.47, 293.09, 293.70, 294.31,
            294.91, 295.50, 296.09, 296.67, 297.24, 297.81, 298.37, 298.92, 299.47, 300.01, 300.54, 301.07, 301.59,
            302.11, 302.61, 303.12, 303.61, 304.10, 304.58, 305.06, 305.53, 306.00],
    'lat': [72.50, 72.59, 72.68, 72.77, 72.85, 72.94, 73.02, 73.10, 73.17, 73.25, 73.32, 73.40, 73.46, 73.53, 73.59,
            73.66, 73.72, 73.77, 73.83, 73.88, 73.93, 73.98, 74.02, 74.07, 74.10, 74.14, 74.18, 74.21, 74.24, 74.26,
            74.29, 74.31, 74.33, 74.34, 74.36, 74.37, 74.37, 74.38, 74.38, 74.38, 74.38, 74.37, 74.36, 74.35, 74.34,
            74.32, 74.30, 74.28, 74.25, 74.23, 74.20, 74.16, 74.13, 74.09, 74.05, 74.01, 73.96, 73.91, 73.86, 73.81,
            73.75, 73.69, 73.63, 73.57, 73.50, 73.44, 73.37, 73.00, 73.22, 73.14, 73.06, 72.98, 72.90, 72.82, 72.73,
            72.64, 72.55, 72.46, 72.36, 72.27, 72.17, 72.07, 71.97, 71.86, 71.76, 71.65, 71.54, 71.43, 71.32, 71.21,
            71.09, 70.97, 70.86, 70.74, 70.62, 70.50, 70.37, 70.25, 70.12, 70.00],
    'str': 'CAA_Transect',
    'ylim_max': 36
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

    rounded_lat = np.round(np.multiply(transect_type['lat'], 2))/2
    index_lat = []
    for x in rounded_lat:
        index_lat.append(np.argwhere(lat == x)[0][0])

    rounded_lon = np.round(np.multiply(transect_type['lon'], 2))/2
    index_lon = []
    for x in rounded_lon:
        index_lon.append(np.argwhere(lon == x)[0][0])

    for j in range(0, len(years)):
        plt.close()
        fig, ax = plt.subplots()
        if len(data.shape) == 4:
            transect_data = data[j, :, index_lat, index_lon]
        else:
            transect_data = data[:, index_lat, index_lon]
        graphed_data = plt.pcolormesh(transect_data, cmap='viridis')
        color_bar = plt.colorbar(graphed_data)
        color_bar.set_label(units)
        plt.xticks([], [])
        plt.ylim([transect_type['ylim_max'], 0])
        yticks = ax.get_yticks()
        for k in range(0, len(yticks)):
            yticks[k] = '{:.2f}'.format(depths[int(yticks[k])])
        ax.set_yticklabels(yticks)
        plot_name = 'plots/{}_{}_{}.pdf'.format(var[i], transect_type['str'], years[j])
        plt.savefig(plot_name, bbox_inches='tight')

    # Don't mind this code, it's to figure out the lat/lon points of the transect
    #
    # trans_x, trans_y = np.linspace(65, 60, 100, dtype=int), np.linspace(480, 612, 100, dtype=int)
    # x, y = m(lon[trans_y], lat[trans_x])
    # trans_x, trans_y = np.linspace(x[0], x[-1], 100), np.linspace(y[0], y[-1], 100)
    # x, y = m(trans_x, trans_y, inverse=True)
    # print 360+x
    # print '-'*10
    # print y
    # plt.plot(trans_x, trans_y, 'g-',)
