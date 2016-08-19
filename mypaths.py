# -*- coding: utf-8 -*-
"""
Paths to data
"""
import glob
import os

opj = os.path.join


def _sg(path):
    return sorted(glob.glob(path))

#
# Common directory
#
homedir = opj(os.getenv('HOME'), 'UEA', 'ACCACIA')
extdir = opj('/media', os.getenv('USER'), 'Elements', 'ACCACIA')

#
# Output directories
#
localfig = 'figures'

#
# UM
#
um_dir = opj(homedir, 'UM', 'exp_results')
um_dir_ext = opj(extdir, 'UM', 'exp_results')

#
# Satellite data
#
# AVHRR
avhrr_dir = opj(homedir, 'satellite', 'all_r8')
avhrr_flist = _sg(opj(avhrr_dir, '*ch4.r8.tif'))
# Polar lows locations (from AVHRR imagery)
trajdir = opj(homedir, 'satellite')
trajf = opj(trajdir, 'pmc_loc_time_ch4_20Mar-02Apr.txt')
# ASCAT
ascat_dir = opj(homedir, 'satellite', 'ASCAT')
ascat_12km_flist = _sg(opj(ascat_dir, '12.5km', '*.nc'))
ascat_25km_flist = _sg(opj(ascat_dir, '25.0km', '*.nc'))
# Cloudsat
cloudsat_dir = opj(homedir, 'satellite', 'CLOUDSAT')
cloudsat_flist = _sg(opj(cloudsat_dir, '*.h5'))
cloudsat_geoprof_flist = _sg(opj(cloudsat_dir, '*CS_2B-GEOPROF_GRANULE*.h5'))
cloudsat_cwc_flist = _sg(opj(cloudsat_dir, '*CS_2B-CWC-RO_GRANULE*.h5'))

#
# In-situ observations: dropsondes and on-board aircraft probes
#
# "Core" data
faamdir = opj(homedir, 'core_processed')
faam_fname_pref = 'core_faam_20130326_v004_r0_b763'
faam_fname_post = '_1hz'
faamf = opj(faamdir, faam_fname_pref + faam_fname_post + '.nc')
# Cloud data
cloud_data_dir = opj(homedir, 'CloudData')
# Dropsondes
ds_pref = 'faam-dropsonde_faam_20130326'
ds_suff = '_r0_b763_proc.nc'
ds_tmark = ['112103', '112518', '112926', '113330', '113737', '114746',
            '115252', '115741', '120230', '120726', '121210']
ds_flist = sorted([opj(faamdir, ds_pref+i+ds_suff) for i in ds_tmark])
# Flux data base
fluxf = opj(faamdir, 'flux_data_32flights.nc')

#
# ERA-Interim
#
erai_in = opj(homedir, 'ERA-Interim', 'ERAI_05')
erai_plev = opj(erai_in, 'erai-20130320-20130403-05deg-plev.nc')
