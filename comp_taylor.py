import os
import taylor
import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmipdata as cd

plt.close('all')
plt.ion()
font = {'size': 12}
plt.rc('font', **font)

runid = 'nrb'
year = '1000'
##The ptrc input file
ifile_ptrc = ('/raid/ra40/data/ncs/nemo_out/' + runid + '/mc_' + runid + '_1m_' + year + '0101_' +
              year + '1231_ptrc_t.nc.001')

ifile_diad = ('/raid/ra40/data/ncs/nemo_out/' + runid + '/mc_' + runid + '_1m_' + year + '0101_' +
              year + '1231_diad_t.nc.001')

ifile_ptrc2 = '/raid/ra40/data/ncs/nemo_out/oda/mc_oda-o2c-cv2_1m_17870101_17871231_ptrc_t.nc.001'

ifile_diad2 = '/raid/ra40/data/ncs/nemo_out/oda/mc_oda-o2c-cv2_1m_17870101_17871231_diad_t.nc.001'

obs_root = '/raid/ra40/data/ncs/nemo_out/obs4comp/'

obs4comp = {'NO3': obs_root + 'uncs_orca2_data_data_n_an_nomask.nc',
            # 'NO3' : obs_root + 'data_NO3_nomask.nc',
            'DIC': obs_root + 'uncs_orca2_data_data_TCO2_nomask.nc',
            'NCHL': obs_root + 'orca2_seawifs_mean_1998_2005.nc',
            'O2': obs_root + 'umol-l_uncs_orca2_data_data_o_an_nomask.nc',
            'Alkalini': obs_root + 'uncs_orca2_data_data_Alk_nomask.nc',
            'Cflx': obs_root + 'orca2_landschuetzer_fgco2.nc',
            'EPC100': obs_root + 'orca2_mol-m2-s_AWI_export.nc',
            }


def corr(obs, data, weights=None):
    """Compute a weighted correlation coefficient and std dev for obs and model

    Returns:
       r : correlation coefficient
       ovar : weighted stddev of obs
       dvar : weighted stddev of data

    """

    if weights is None:
        weights = np.ones(obs.shape)

    obsf = obs.flatten()
    dataf = data.flatten()
    weightsf = weights.flatten()

    obar = np.ma.average(obsf, weights=weightsf)
    dbar = np.ma.average(dataf, weights=weightsf)
    ovar = np.sqrt(np.ma.average((obsf - obar) ** 2, weights=weightsf))
    dvar = np.sqrt(np.ma.average((dataf - dbar) ** 2, weights=weightsf))

    r = np.ma.average((obsf - obar) * (dataf - dbar), weights=weightsf) / (ovar * dvar)

    return r, ovar, dvar


def taylor_plot(datalist, obs, weights=None, fig=None, rect=111, collist=None,
                symlist=None, lablist=None, title=None, srange=(0, 1.5)):
    if not fig:
        fig = plt.figure(figsize=(8, 4))

    print "Computing Taylor diagram statistics..."

    corrcoef, refstd, stddev = corr(obs, obs, weights)

    # Taylor diagram
    dia = taylor.TaylorDiagram(refstd, fig=fig, rect=rect, label="Obs.",
                               title=title, srange=srange)

    if not collist:
        colors = plt.matplotlib.cm.jet(np.linspace(0, 1, len(datalist)))

    if not lablist:
        lablist = ['' for i in datalist]

    if not symlist:
        symlist = ['o' for i in datalist]

    for i, data in enumerate(datalist):
        corrcoef, refstd, stddev = corr(obs, data, weights)
        print "R:", corrcoef
        print "Std, obs, model:", refstd, stddev

        # Add the models to Taylor diagram
        dia.add_sample(stddev, corrcoef, marker=symlist[i], ms=10, ls='',
                       mfc=collist[i], mec=collist[i], label=lablist[i])

    # Add grid
    dia.add_grid(True, axis='x', linestyle='--', alpha=0.6,
                 color=[0.5, 0.5, 0.5], zorder=1)

    # Add RMS contours, and label them
    contours = dia.add_contours(colors='0.5')
    plt.clabel(contours, inline=1, fontsize=10)

    # Add a figure legend
    fig.legend(dia.samplePoints,
               [p.get_label() for p in dia.samplePoints],
               numpoints=1, prop=dict(size='small'), loc='upper right',
               frameon=False)


def taylor_plot_depths(datalist, obs, weights=None, fig=None, rect=111,
                       collist=None, lablist=None, title=None, srange=(0, 1.5)):
    if not fig:
        fig = plt.figure(figsize=(8, 4))

    print "Computing Taylor diagram statistics..."

    # Taylor diagram
    dia = taylor.TaylorDiagram(1, fig=fig, rect=rect, label="Obs.",
                               title=title, srange=srange)

    if not collist:
        collist = plt.matplotlib.cm.jet(np.linspace(0, 1, len(datalist)))

    if not lablist:
        lablist = ['' for i in datalist]

    for i, data in enumerate(datalist):
        print obs.shape, data.shape

        for k in range(data.shape[0]):
            corrcoef, refstd, stddev = corr(obs[k, :, :], data[k, :, :], weights)

            # Add the models to Taylor diagram
            dia.add_sample(stddev, corrcoef, marker='$%d$' % (k + 1),
                           normalize=refstd,
                           ms=10, ls='', mfc=collist[i], mec=collist[i],
                           label=lablist[i])

    # Add grid
    dia.add_grid(True, axis='x', linestyle='--', alpha=0.6,
                 color=[0.5, 0.5, 0.5], zorder=1)

    # Add RMS contours, and label them
    contours = dia.add_contours(colors='0.5')
    plt.clabel(contours, inline=1, fontsize=10)


def load(ifile_fp, varname, gridfile):
    """Time average ifile  and return var, dims"""

    path, ifile = os.path.split(ifile_fp)
    data = cd.loadvar(ifile_fp, varname, cdostr='timmean')
    tmask = cd.loadvar(gridfile, 'tmask')

    if data.ndim == 3:
        data = data * tmask
    else:
        data = data * tmask[0, :, :].squeeze()

    return np.ma.masked_equal(data, 0)


def genweights(gridfile):
    """Return weights based on grid area (2D) and
       volume (3D).
    """
    # Load grid dimensions and tmask
    e1t = cd.loadvar(gridfile, 'e1t')
    e2t = cd.loadvar(gridfile, 'e2t')
    e3t = cd.loadvar(gridfile, 'e3t')
    tmask = cd.loadvar(gridfile, 'tmask')

    # tile e1 and e2 areas to 3D
    e1t3 = np.tile(e1t, [e3t.shape[0], 1, 1])
    e2t3 = np.tile(e2t, [e3t.shape[0], 1, 1])

    vol = e1t3 * e2t3 * e3t * tmask
    vol = np.ma.masked_equal(vol, 0)
    weights3d = vol / vol.sum()

    area = e1t.squeeze() * e2t.squeeze() * tmask[0, :, :]
    area = np.ma.masked_equal(area, 0)
    weights2d = area / area.sum()

    return weights2d, weights3d


#def driver(varlist, ifilelist, obs4comp):
#"""Calls taylor_plot for each variable in varlist,
#and for each variable plots a point for obs in
#obs4comp and model in ifilelist
#"""
#for var in varlist:
#obs = load(obs4comp[var], var)

gridfile = '/home/ncs/ra40/nemo_out/nemo_3.4_orca2_mesh_mask.nc'
weights2d, weights3d = genweights(gridfile)

fig = plt.figure(figsize=(8, 8))

canoe_no3 = load(ifile_ptrc, 'NO3', gridfile)
cmoc_no3 = load(ifile_ptrc2, 'NO3', gridfile)
obs_no3 = load(obs4comp['NO3'], 'NO3', gridfile)
taylor_plot([canoe_no3, cmoc_no3], obs_no3, weights=weights3d,
            fig=fig, rect=221, lablist=['CanOE', 'CMOC'], symlist=['o', 'o'],
            title='NO3', collist=['r', 'b'])

canoe_dic = load(ifile_ptrc, 'DIC', gridfile)
cmoc_dic = load(ifile_ptrc2, 'DIC', gridfile)
obs_dic = load(obs4comp['DIC'], 'DIC', gridfile)
taylor_plot([canoe_dic, cmoc_dic], obs_dic, weights=weights3d,
            fig=fig, rect=222, lablist=['CanOE', 'CMOC'], symlist=['o', 'o'],
            title='DIC', collist=['r', 'b'])

canoe_NCHL = load(ifile_ptrc, 'NCHL', gridfile)
canoe_DCHL = load(ifile_ptrc, 'DCHL', gridfile)
canoe_chl = canoe_NCHL + canoe_DCHL
cmoc_NCHL = load(ifile_ptrc2, 'NCHL', gridfile)
obs_NCHL = load(obs4comp['NCHL'], 'NCHL', gridfile)
taylor_plot([canoe_chl[0, :, :], cmoc_NCHL[0, :, :]], obs_NCHL,
            weights=weights2d, fig=fig, rect=223, lablist=['CanOE', 'CMOC'],
            symlist=['o', 'o'], title='CHL', srange=(0, 2), collist=['r', 'b'])

canoe_cflx = load(ifile_diad, 'Cflx', gridfile) * 1e8
cmoc_cflx = load(ifile_diad2, 'Cflx', gridfile) * 1e8
obs_cflx = load(obs4comp['Cflx'], 'Cflx', gridfile) * 1e8
taylor_plot([canoe_cflx, cmoc_cflx], obs_cflx, weights=weights2d,
            fig=fig, rect=224, lablist=['CanOE', 'CMOC'], symlist=['o', 'o'],
            title='Cflx', collist=['r', 'b'])

plt.subplots_adjust(hspace=0.45)
plt.savefig('nemo_Taylor.png', dpi=300)
##############3
fig = plt.figure(figsize=(8, 8))

taylor_plot_depths([canoe_no3, cmoc_no3], obs_no3, weights=weights2d, fig=fig,
                   rect=121, collist=['r', 'b'], lablist=['CanOE', 'CMOC'],
                   title='NO3', srange=(0, 2))

taylor_plot_depths([canoe_dic, cmoc_dic], obs_dic, weights=weights2d, fig=fig,
                   rect=122, collist=['r', 'b'], lablist=['CanOE', 'CMOC'],
                   title='DIC', srange=(0, 2))

#taylor_plot_depths([canoe_no3, cmoc_no3], obs_no3, weights=weights2d, fig=fig,
#rect=223, collist=['r', 'b'], lablist=['CanOE', 'CMOC'],
#title='NO3', srange=(0,4))

#taylor_plot_depths([canoe_no3, cmoc_no3], obs_no3, weights=weights2d, fig=fig,
#rect=224, collist=['r', 'b'], lablist=['CanOE', 'CMOC'],
#title='NO3', srange=(0,4))

plt.savefig('nemo_Taylor_depths.png', dpi=300)

#################3

runid = 'nue'
year = '0001'
##The ptrc input file
ifile_ptrc2 = ('/raid/ra40/data/ncs/nemo_out/' + runid + '/'
                                                         'mc_' + runid + '_1m_' + year + '0101_' +
               year + '1231_ptrc_t.nc.001')

cmoc_no3 = load(ifile_ptrc2, 'NO3', gridfile)

fig = plt.figure(figsize=(8, 8))
taylor_plot_depths([cmoc_no3], obs_no3, weights=weights2d, fig=fig,
                   rect=111, collist=['b'], lablist=['CMOC'],
                   title='NO3', srange=(0, 2))
