# -*- coding: utf-8 -*-
import datetime
import numpy as np
from scipy import signal
import scipy.interpolate as scinter
# My modules
from phys_meteo import r_earth
#
# TODO: doc strings!
#


def minmax(x):
    return np.array((np.nanmin(x), np.nanmax(x)))


def scl2latex(x):
    print('DEPRECATED; use unit_format!')
    return r'$\times10^{{{0}}}$'.format(-int(np.log10(x)))


def unit_format(value, unit='1'):
    if value == 1:
        if unit == '1':
            string = ''
        else:
            string = unit
    else:
        exp = np.floor(np.log10(value))
        base = value / 10**exp
        if exp == 0 or exp == 1:
            string = r'${0}$'.format(value)
        elif exp == -1:
            string = r'${0:0.1f}$'.format(value)
        else:
            if int(base) == 1:
                string = r'$10^{{{0:d}}}$'.format(int(exp))
            else:
                string = r'${0:d}\times10^{{{1:d}}}$'.format(int(base),
                                                             int(exp))
        if not unit == '1':
            string += r' ${}$'.format(unit)

    return string


def domain_corners_ll(lon, lat):
    def _get_corners(a):
        return a[::a.shape[0]-1, ::a.shape[1]-1]
    model_dom_full = [(i, j) for i, j in zip(_get_corners(lon).flatten(),
                                             _get_corners(lat).flatten())]
    model_dom_full = model_dom_full + [model_dom_full[0]]
    model_dom_full[2], model_dom_full[3] = model_dom_full[3], model_dom_full[2]
    x = [i[0] for i in model_dom_full]
    y = [i[1] for i in model_dom_full]
    return x, y


def ticks_format(value, index):
    """
    get the value and returns the value as:
       integer: [0,99]
       1 digit float: [0.1, 0.99]
       n*10^m: otherwise
    To have all the number of the same size they are all returned
    as latex strings
    """
    exp = np.floor(np.log10(value))
    base = value/10**exp
    if exp == 0 or exp == 1:
        return '${0}$'.format(value)
    if exp == -1:
        return '${0:0.1f}$'.format(value)
    else:
        if int(base) == 1:
            return '$10^{{{0:d}}}$'.format(int(exp))
        else:
            return '${0:d}\\times10^{{{1:d}}}$'.format(int(base), int(exp))


def interp_nan(x):
    def _indices(z):
        return z.nonzero()[0]
    nans = np.isnan(x)
    x[nans] = np.interp(_indices(nans), _indices(~nans), x[~nans])
    return x


def butter_lowpass(cutoff, sample_rate, order):
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def blf(data, cutoff=0.01, sample_rate=1., order=5., interpnan=True):
    b, a = butter_lowpass(cutoff, sample_rate, order)
    if interpnan:
        data = interp_nan(data)
    y = signal.filtfilt(b, a, data)
    return y


def timestr2datetime(tref, fmt='%Y-%m-%d%H:%M:%S'):
    """
    Parse a timestring and return a datetime object
    """
    spl = tref.rsplit()
    if spl[-1][0] == '+' or spl[-1] == 'UTC':
        date_and_time = spl[-3] + spl[-2]
    else:
        date_and_time = spl[-2] + spl[-1]
    dt = datetime.datetime.strptime(date_and_time, fmt)

    if spl[0].lower() == 'days':
        nsec = 86400
    elif spl[0].lower() == 'hours':
        nsec = 3600
    elif spl[0].lower() == 'minutes':
        nsec = 60
    elif spl[0].lower() == 'seconds':
        nsec = 1
    else:
        nsec = 0
        print("Warning: unrecognized time unit, returned 0 for timestep")

    return dt, nsec


def generic_regz(zz, data_list_of_dict):
    regz = np.zeros((len(zz), len(data_list_of_dict)))
    for j, i in enumerate(data_list_of_dict):
        regz[:, j] = scinter.interp1d(i['coord'],
                                      i['fld'], bounds_error=False)(zz)

    return regz


def generic_cross_sect(xgrid, zz, dist, data_list_of_dict):
    regz = generic_regz(zz, data_list_of_dict)

    cross = np.zeros((len(zz), len(xgrid)))

    for i in range(len(zz)):
        cross[i, :] = scinter.interp1d(dist, regz[i, :],
                                       bounds_error=False)(xgrid)

    return cross


def ds_regz(zz, ds_id, data, coordname, varname, filt_kw):
    j = 0
    regz = np.zeros((len(zz), len(ds_id)))
    for i in ds_id:
        arr = getattr(data[i], varname).fil[:]
        if isinstance(filt_kw, dict):
            arr = blf(arr, **filt_kw)
        regz[:, j] = scinter.interp1d(getattr(data[i], coordname).fil[:],
                                      arr,
                                      bounds_error=False)(zz)
        j += 1

    return regz


def ds_cross_sect(xgrid, zz, leg_dist, ds_id, data, coordname, varname,
                  filt_kw):
    regz = ds_regz(zz, ds_id, data, coordname, varname, filt_kw)

    cross = np.zeros((len(zz), len(xgrid)))

    for i in range(0, len(zz)):
        cross[i, :] = scinter.interp1d(leg_dist, regz[i, :],
                                       bounds_error=False)(xgrid)

    return cross


def get_gridded(xpoints, ypoints, values, xx, yy):
    gr = scinter.griddata((xpoints, ypoints), values, (xx, yy),
                          method='linear')

    return gr


def get_leggrid(ds_id, data, xstep=1):
    dslat0 = np.mean(data[ds_id[0]].lat.fil[:])
    dslon0 = np.mean(data[ds_id[0]].lon.fil[:])

    leg_dist = np.zeros((len(ds_id)))
    leg_dist[0] = 0.
    j = 1
    for i in ds_id[1::]:
        dslat = np.mean(data[i].lat.fil[:])
        dslon = np.mean(data[i].lon.fil[:])
        leg_dist[j] = r_earth/1e3*np.arccos(np.sin(np.radians(dslat0))
                                            * np.sin(np.radians(dslat))
                                            + np.cos(np.radians(dslat0))
                                            * np.cos(np.radians(dslat))
                                            * np.cos(np.radians(dslon0-dslon))
                                            )
        j += 1

    leg_grid = np.arange(np.min(leg_dist), np.max(leg_dist), xstep)

    return leg_dist, leg_grid
