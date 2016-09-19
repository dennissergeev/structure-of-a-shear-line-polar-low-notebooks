# -*- coding: utf-8 -*-
"""
Collection of useful plotting functions
"""
import LatLon23 as LL
import matplotlib as mpl
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import numpy as np
import string
# My modules
import map_plot_func as mymap
import plot_params as pp
import phys_meteo as met
import misc_utils as misc
# TODO: refactoring, docstrings and PEP8


def add_subplot_labels(axs, xy=(0.05, 0.95), label_type='alpha'):
    if label_type == 'alpha':
        label_list = list(string.ascii_lowercase)

    for iax, label in zip(axs, label_list):
        iax.annotate('('+label+')',
                     xy, xycoords='axes fraction',
                     fontsize=8, fontweight='bold',
                     va='center', ha='center')


def add_hrow(parent_ax, cols=[],
             shift_factor=0.5, height_factor=0.25,
             fc='k', bc='w', fontsize=12, linewidth=1):
    newax = parent_ax.twinx()
    newax.grid('off')

    newax.xaxis.set_visible(False)
    newax.yaxis.set_visible(False)
    newax.patch.set_facecolor(bc)
    newax.set_ylim(0, 1)
    newax.axhline(0, linewidth=linewidth, color=fc)
    newax.axhline(1, linewidth=linewidth, color=fc)
    pos = newax.get_position()
    pos2 = [pos.x0, pos.y0 + shift_factor * pos.height,
            pos.width, height_factor * pos.height]
    newax.set_position(pos2)

    x0, x1 = newax.get_xlim()
    for icol in cols:
        if 'x0' in icol:
            x0 = icol['x0']
            newax.axvline(x0, linewidth=linewidth, color=fc)
        if 'x1' in icol:
            x1 = icol['x1']
            newax.axvline(x1, linewidth=linewidth, color=fc)
        x = 0.5*(x0 + x1)
        y = 0.5
        newax.annotate(icol['text'], xy=(x, y),
                       color=fc, fontsize=fontsize,
                       va='center', ha='center', zorder=100)
    return newax


def horiz_sect_model_obs(pack, plttype='default', clevs=None,
                         av_slicing='tz', add_colorbar=True, fill=True,
                         add_labels=None):
    if 'mapbounds' not in pack and 'model' in pack:
        lon_min = np.floor(pack['model'][0].min())
        lon_max = np.ceil(pack['model'][0].max())
        lat_min = np.floor(pack['model'][1].min())
        lat_max = np.ceil(pack['model'][1].max())
    else:
        lon_min, lon_max, lat_min, lat_max = pack['mapbounds']

    if 'ax' in pack:
        if pack['ax'] is not None:
            ax = pack['ax']
    else:
        fig, ax = plt.subplots()

    bm = mymap.make_map(ax=ax,
                        lon1=lon_min, lon2=lon_max, lat1=lat_min, lat2=lat_max,
                        tick_incr=[2., 1.], fill=fill, scale=True)
    if add_labels is not None:
        for p in add_labels:
            ix, iy = bm(p['lon'], p['lat'])
            ax.annotate(p['name'], (ix, iy), size=18, ha="center",
                        path_effects=[PathEffects.withStroke(linewidth=2,
                                                             foreground="w")])
    limname = pack['lims'][0]
    limmin = pack['lims'][1]
    limmax = pack['lims'][2]
    limmean = 0.5*(limmin+limmax)
    # print(limname, limmin, limmax)

    # Dropsondes
    lon_ds, lat_ds, var_ds = [], [], []
    if 'ds_obs' in pack:
        ds_list, ds_var = pack['ds_obs'][0], pack['ds_obs'][1]
        if isinstance(ds_var, str):
            ds_varname = ds_var
            ds_var = []
            for i in ds_list:
                ds_var.append(getattr(i, ds_varname).fil)

        for j, i in enumerate(ds_list):
            inds = np.where((limmin <= getattr(i, limname).fil)
                            & (getattr(i, limname).fil <= limmax))[0]
            if len(inds) > 0:
                lon_ds.append(np.mean(i.lon.fil[inds]))
                lat_ds.append(np.mean(i.lat.fil[inds]))
                var_ds.append(np.mean(ds_var[j][inds]))

    # Aircraft data
    lon_faam, lat_faam, var_faam = [], [], []
    if 'faam_obs' in pack:
        faam = pack['faam_obs'][0]
        n = 20.  # Optimal number of points to show
        inds = np.where((limmin <= getattr(faam, limname).val)
                        & (getattr(faam, limname).val <= limmax)
                        & (~np.isnan(faam.lon.val))
                        & (~np.isnan(faam.lat.val))
                        & (faam.lat.val > 72.))[0]
        if len(inds) > 0:
            nn = int(np.ceil(len(inds)/np.float(n)))
            lon_faam = faam.lon.val[inds[::nn]]
            lat_faam = faam.lat.val[inds[::nn]]
            var_faam = pack['faam_obs'][1][inds[::nn]]

    # Slice model data
    if plttype == 'default':
        model_data = pack['model'][2]
        if model_data.ndim == 3:
            # Assume 'zyx' dimensions
            zcoord = pack['model'][3]
            ilev = np.argmin(abs(zcoord.points-limmean))
            if av_slicing == 'z':
                ilev1 = np.argmin(abs(zcoord.points-limmin))
                ilev2 = np.argmin(abs(zcoord.points-limmax))
                ilev1, ilev2 = sorted([ilev1, ilev2])
                model_data2d = np.mean(model_data[ilev1:ilev2+1, ...], 0)
            else:
                model_data2d = model_data[ilev, ...]
            model_data = model_data2d
        elif model_data.ndim == 4:
            zcoord = pack['model'][3]
            tcoord = pack['model'][4]
            ilev = np.argmin(abs(zcoord.points-limmean))
            it1 = np.argmin(abs(tcoord.points
                                - tcoord.units.date2num(pack['ds_obs'][0][0].time.fil[0])))
            it2 = np.argmin(abs(tcoord.points-tcoord.units.date2num(pack['ds_obs'][0][-1].time.fil[0])))
            if 'z' in av_slicing and 't' in av_slicing:
                ilev1 = np.argmin(abs(zcoord.points-limmin))
                ilev2 = np.argmin(abs(zcoord.points-limmax))
                ilev1, ilev2 = sorted([ilev1, ilev2])
                model_data2d = np.mean(model_data[it1:it2+1, ...], 0)
                model_data2d = np.mean(model_data2d[ilev1:ilev2+1, ...], 0)
            elif av_slicing == 't':
                model_data2d = np.mean(model_data[it1:it2+1, ...], 0)
                model_data2d = model_data2d[ilev, ...]
            elif av_slicing == 'z':
                ilev1 = np.argmin(abs(zcoord.points-limmin))
                ilev2 = np.argmin(abs(zcoord.points-limmax))
                ilev1, ilev2 = sorted([ilev1, ilev2])
                model_data2d = np.mean(model_data[int(round(0.5*(it1 + it2))),
                                                  ilev1:ilev2+1, ...], 0)
            else:
                model_data2d = model_data[int(round(0.5*(it1 + it2))),
                                          ilev, ...]
            model_data = model_data2d

        assert model_data.ndim == 2, 'Error: model_data array has to be of rank 2, but the shape is {0} instead'.format(model_data.shape)

        if not all(len(i)==0  for i in (lon_faam, lat_faam, var_faam, lon_ds, lat_ds, var_ds)):
            lon_all = np.append(lon_ds, lon_faam)
            lat_all = np.append(lat_ds, lat_faam)
            var_all = np.append(var_ds, var_faam)
            c = map_model_obs_scatter(ax, bm, pack['model'][0], pack['model'][1], model_data2d, lon_all, lat_all, var_all, clevs=clevs, add_colorbar=add_colorbar)
            if not add_colorbar:
                return c
        else:
            print('No observations sampled, plotting only model data')
            xx, yy = bm(pack['model'][0], pack['model'][1])
            c = bm.contourf(xx, yy, model_data2d, clevs)
            if add_colorbar:
                plt.colorbar(c, ax=ax)

    elif plttype == 'barbs':
        # Assume that 3rd element of 'model' entry is a list or tuple with u and v
        u_and_v = []
        for iuv in pack['model'][2]:
            model_data = iuv
            if len(model_data.shape) == 4:
                ilev = int(np.argmin(abs(zcoord.points-limmean)))
                it1 = int(np.argmin(abs(tcoord.points-tcoord.units.date2num(ds[0].time.fil[0]))))
                it2 = int(np.argmin(abs(tcoord.points-tcoord.units.date2num(ds[-1].time.fil[0]))))
                if av_slicing == 'tz':
                    model_data2d = np.mean(model_data[it1:it2+1, ...], 0)
                    model_data2d = np.mean(model_data2d, 0)
                elif av_slicing == 't':
                    model_data2d = np.mean(model_data[it1:it2+1, ...], 0)
                    model_data2d = model_data2d[ilev,...]
                elif av_slicing == 'z':
                    model_data2d = np.mean(model_data[int(round(0.5*(it1 + it2))), ilev-1:ilev+1, ...], 0)
                model_data = model_data2d
            u_and_v.append(model_data)


def map_model_obs_scatter(ax, bm, model_lons, model_lats, model_data, obs_lons, obs_lats, obs_data, clevs=None, cmap=plt.cm.Oranges, add_colorbar=True):
    if clevs is None:
        clevs = np.linspace(min(model_data.min(),obs_data.min()), max(model_data.max(),obs_data.max()),100)
    xx, yy = bm(model_lons, model_lats)
    c = bm.contourf(xx, yy, model_data, clevs, cmap=cmap)
    xx, yy = bm(obs_lons, obs_lats)
    sc = bm.scatter(xx, yy, c=obs_data, cmap=cmap, s=2**7, rasterized=True, vmin=clevs[0], vmax=clevs[-1], edgecolors='w')
    if add_colorbar:
        plt.colorbar(sc, ax=ax, shrink=0.75)
    else:
        return c


def map_model_obs_barbs(ax, bm, model_lons, model_lats, model_uv, obs_lons, obs_lats, obs_uv, s=50):
    if hasattr(model_uv,'__iter__') and len(model_uv)==2 and \
       hasattr(obs_uv,'__iter__') and len(obs_uv)==2:
        xx, yy = bm(lon2d, lat2d)
        b1 = bm.barbs(xx[::s,::s], yy[::s,::s], model_uv[0][::s,::s], model_uv[1][::s,::s], color='r')
        xx, yy = bm(obs_lons, obs_lats)
        b2 = bm.barbs(xx, yy, obs_uv[0], obs_uv[1], color=pp.almost_black)
    else:
        raise ValueError('Input variables should consist of u- and v-wind pair')


def pcol_n_tseries_vs_hgt(varlist, mask=None, lon=None, lat=None, add_dist=False, time=None, fs=None, llstep=60, invert_xax=False, segments=None):

    if segments is None:
        nrows = len(varlist)
    else:
        nrows = len(varlist)+1

    if fs is None:
        fs = (12, nrows*4) # +1 to add segments
    
    fig = plt.figure(figsize=fs)
    axgr = AxesGrid(fig, 111, nrows_ncols=(nrows, 1), axes_pad=0.4, aspect=False,
                    cbar_location='right', cbar_mode='each', cbar_pad=0.1)

    
    cbars = []
    sbplt_labels = list(string.ascii_lowercase)
    for n, (ax, ivar, lab) in enumerate(zip(axgr.axes_all[1:], varlist, sbplt_labels)): 
        
        x = ivar['x']
        z = ivar['z']
        data = ivar['data']
        kw = ivar['kw']
        im = ax.pcolormesh(x, z, data, **kw)

        ax.set(xlim=[x[0], x[-1]], ylim=[z[0], z[-1]])
        ax.set_ylabel('Height (km)',fontsize=18)
        ax.tick_params(axis='y', which='major', labelsize=16)
        ax.get_yaxis().set_label_coords(-0.025,0.5)

        at = AnchoredText('({lab}) {title}'.format(lab=lab, title=ivar['ttl']), prop=dict(size=18), frameon=True, loc=2)
        at.patch.set_boxstyle('round', pad=0., rounding_size=0.2)
        at.set_zorder(400)
        ax.add_artist(at)

        cb = plt.colorbar(im, cax=axgr.cbar_axes[n+1])
        cb.ax.tick_params(labelsize=18)
        if 'cbkw' in ivar:
            cb.ax.set_yticklabels(ivar['cbkw']['ticks'])
        cb.ax.set_title(ivar['units'], fontsize=18)
        cbars.append(cb)
        #cb.ax.tick_params(labelsize=22)
        #    cb.set_ticks(np.arange(0,len(ivar['cb_ticks'])))
        #    cb.set_ticklabels(ivar['cb_ticks'])
            
        if invert_xax:
            ax.invert_xaxis()
            
        #for tick in ax.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(16) 
            
        ax.xaxis.set_visible(False)
        ax.spines['bottom'].set_color('w')
        _ax = add_xaxis_below(ax, np.linspace(0,1,len(time[::llstep])), [], 0)

    # Additional axes
    if time is not None:
        time =  time[::llstep]
        tticks = np.linspace(0,1,len(time))
        xlabels_t = [i.strftime('%H:%M:%S') for i in time]
        newax_t = add_xaxis_below(ax, tticks, xlabels_t, 0)
        if invert_xax:
            newax_t.invert_xaxis()
    if lon is not None and lat is not None:
        lon, lat = lon[::llstep], lat[::llstep]
        llticks = np.linspace(0,1,len(lon))
        xlabels_lat = [LL.Latitude(i).to_string("d%$^\circ$%m%'%H") for i in lat]
        xlabels_lon = [LL.Longitude(i).to_string("d%$^\circ$%m%'%H") for i in lon]
        newax_lon = add_xaxis_below(ax, llticks, xlabels_lon, 30)
        newax_lat = add_xaxis_below(ax, llticks, xlabels_lat, 60)
        if invert_xax:
            newax_lon.invert_xaxis()
            newax_lat.invert_xaxis()

        if add_dist:
            dist = [0]
            for i in range(1, len(lon)):
                dist.append(
                            met.r_earth/1000.*np.arccos( \
                            np.sin(np.radians(lat[i-1]))*np.sin(np.radians(lat[i])) + \
                            np.cos(np.radians(lat[i-1])) \
                            *np.cos(np.radians(lat[i]))*np.cos(np.radians(lon[i-1]-lon[i]))
                                                       )
                            )
            newax_dist = add_xaxis_below(ax, llticks, ['{0:3.0f} km'.format(i) for i in np.cumsum(dist)], 90)
            if invert_xax:
                newax_dist.invert_xaxis()
                
    for n_ax in (newax_t, newax_lon, newax_lat, newax_dist):
        n_ax.tick_params(axis='x', which='major', labelsize=16)
                
    if segments is not None:
        axgr[0].set_ylim(0, 2)
        add_hrow(axgr[0], **segments)
        axgr[0].axis('off')
        axgr.cbar_axes[0].remove()
        
    return fig, axgr, cbars


def plot_n_zprofiles(varlist, fs=None, filt_keys=dict(sample_rate=1., cutoff=1/100., order=6), legend=True):
    
    if fs is None:
        fs = (30,10)

    fig, axs = plt.subplots(ncols=len(varlist), sharey=True, figsize=fs)
    if len(varlist) == 1:
        axs = [axs]
    
    
    for n, ytup in enumerate(varlist):
        lnkw_data = []
        lnkw_filt = []
        for i, idic in enumerate(ytup[:-1]): # The last ytup is a dict with plot parameters
            lnkw_data.append(dict(color=idic['c'] if 'c' in idic else pp.cvec[i], 
                                  linewidth=idic['lw'] if 'lw' in idic else 1.5, alpha=0.6))
            lnkw_filt.append(dict(color=idic['c'] if 'c' in idic else pp.cvec[i], 
                                  linewidth=idic['lw'] if 'lw' in idic else 2.5))
        try:
            if ytup[-1]['plttype'].lower() == 'barbs':
                fake_handles = []
                for i, idic in enumerate(ytup[:-1]):
                    s = idic['data'][-1] if isinstance(idic['data'][-1], int) else 1
                    
                    axs[n].barbs(idic['data'][0][::s], idic['yax'][::s],
                          misc.blf(idic['data'][1], **filt_keys)[::s],
                          misc.blf(idic['data'][2], **filt_keys)[::s],
                          label=idic['label'] if 'label' in idic else 'data{}'.format(0),
                          barb_increments=pp.barb_incr,
                          color=lnkw_filt[i]['color'])
                    fake_handles.append(axs[n].plot([], [],
                                        color=lnkw_filt[i]['color'],
                                        linewidth=idic['lw'] if 'lw' in idic else 1)[0])
                _, labels = axs[n].get_legend_handles_labels()
                if legend:
                    axs[n].legend(fake_handles, labels, loc=1)

        except KeyError:
            for i, idic in enumerate(ytup[:-1]):
                hgt = idic['yax']
                if 'apply_filt' in idic:
                    if idic['apply_filt']:
                        axs[n].plot(idic['data'], hgt, **lnkw_data[i])
                        axs[n].plot(misc.blf(idic['data'], **filt_keys), hgt,
                                    label=idic['label'] if 'label' in idic else 'data{}'.format(0),
                                    **lnkw_filt[i])
                    else:
                        axs[n].plot(idic['data'], hgt,
                                    label=idic['label'] if 'label' in idic else 'data{}'.format(0),
                                    **lnkw_filt[i])
                else:
                    axs[n].plot(idic['data'], hgt,
                                label=idic['label'] if 'label' in idic else 'data{}'.format(0),
                                **lnkw_filt[i])

            if legend:
                axs[n].legend(loc=1)

        axs[n].set_title(ytup[-1]['ttl'], fontsize=20)
        axs[n].set_xlabel(ytup[-1]['xlab'], fontsize=20)
        axs[n].set_ylim(0, 6) #np.nanmin(hgt), np.nanmax(hgt))
        axs[n].grid('on')
        axs[n].tick_params(axis='both', which='major', labelsize=20)
            
    axs[0].set_ylabel('Height, km', fontsize=20)
     
    fig.tight_layout()
        
    return fig, axs

        
def plot_n_tseries(varlist,
                   time=None, lon=None, lat=None, add_dist=False, llstep=60, segments=None,
                   fs=None, filt_keys=dict(sample_rate=1., cutoff=1/100., order=6)):
    
    if fs is None:
        fs = (27,(len(varlist))*3)

    fig, axs = plt.subplots(nrows=len(varlist), sharex=True, figsize=fs)
    if len(varlist) == 1:
        axs = [axs]
    
    for n, ytup in enumerate(varlist):
        ax = axs[n]
        
        lnkw_data = []
        lnkw_filt = []
        for i, idic in enumerate(ytup[:-1]): # The last ytup is a dict with plot parameters
            lnkw_data.append(dict(color=idic['c'] if 'c' in idic else pp.cvec[i], linewidth=1., alpha=0.6))
            lnkw_filt.append(dict(color=idic['c'] if 'c' in idic else pp.cvec[i], linewidth=2.))
        try:
            if ytup[-1]['plttype'].lower() == 'barbs':
                fake_handles = []
                for i, idic in enumerate(ytup[:-1]):
                    mdt = idic['xax']
                    s = idic['data'][-1] if isinstance(idic['data'][-1], int) else 1

                    ax.barbs(idic['xax'][::s], idic['data'][0][::s],
                          misc.blf(idic['data'][1], **filt_keys)[::s],
                          misc.blf(idic['data'][2], **filt_keys)[::s],
                          label=idic['label'] if 'label' in idic else 'data{}'.format(0),
                          barb_increments=pp.barb_incr,
                          color=lnkw_filt[i]['color'])
                    fake_handles.append(ax.plot([], [],
                                                color=lnkw_filt[i]['color'])[0])
                _, labels = ax.get_legend_handles_labels()
                #ax.legend(fake_handles, labels, loc=1)

        except KeyError:
            for i, idic in enumerate(ytup[:-1]):
                mdt = idic['xax']
                if 'apply_filt' in idic:
                    if idic['apply_filt']:
                        ax.plot(mdt, idic['data'], **lnkw_data[i])
                        ax.plot(mdt, misc.blf(idic['data'], **filt_keys),
                                label=idic['label'] if 'label' in idic else 'data{}'.format(0),
                                **lnkw_filt[i])
                    else:
                        ax.plot(mdt, idic['data'], label=idic['label'] if 'label' in idic else 'data{}'.format(0), **lnkw_filt[i])
                else:
                    ax.plot(mdt, idic['data'], label=idic['label'] if 'label' in idic else 'data{}'.format(0), **lnkw_filt[i])

            #ax.legend(loc=1)
    
        ax.set_title(ytup[-1]['ttl'], fontsize=18)
        ax.set_ylabel(ytup[-1]['ylab'], fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.get_yaxis().set_label_coords(-0.06, 0.5)
        ax.set_xlim(mdt[0],mdt[-1])
        axs[n].grid('on')
        if 'pres' in ytup[-1]['ttl'].lower():
            ax.invert_yaxis()

        ax.xaxis.set_visible(False)
        ax.spines['bottom'].set_color('w')
        _ax = add_xaxis_below(ax, np.linspace(0, 1, len(time[::llstep])), [], 0)

    # Additional axes
    if time is not None:
        time =  time[::llstep]
        tticks = np.linspace(0,1,len(time))
        xlabels_t = [i.strftime('%H:%M') for i in time]
        newax_t = add_xaxis_below(ax, tticks, xlabels_t, 0)

    if lon is not None and lat is not None:
        lon, lat = lon[::llstep], lat[::llstep]
        llticks = np.linspace(0,1,len(lon))
        xlabels_lat = [LL.Latitude(i).to_string("d%$^\circ$%m%'%H") for i in lat]
        xlabels_lon = [LL.Longitude(i).to_string("d%$^\circ$%m%'%H") for i in lon]
        newax_lon = add_xaxis_below(ax, llticks, xlabels_lon, 20)
        newax_lat = add_xaxis_below(ax, llticks, xlabels_lat, 40)

        if add_dist:
            dist = [0]
            for i in range(1, len(lon)):
                dist.append(
                            met.r_earth/1000.*np.arccos( \
                            np.sin(np.radians(lat[i-1]))*np.sin(np.radians(lat[i])) + \
                            np.cos(np.radians(lat[i-1])) \
                            *np.cos(np.radians(lat[i]))*np.cos(np.radians(lon[i-1]-lon[i]))
                                                       )
                            )
            newax_dist = add_xaxis_below(ax, llticks,
                                         ['{0:3.0f} km'.format(i) for i in np.cumsum(dist)], 60)

    if segments is not None:
        add_hrow(ax, **segments)

    return fig, axs


def add_xaxis_below(parent_ax,xtick_array,xlab_array,shift_down):
    newax = parent_ax.twiny()
    newax.set_xticks(xtick_array)
    newax.set_xticklabels(xlab_array)
    newax.spines['left'].set_visible(False)
    newax.spines['right'].set_visible(False)
    newax.set_frame_on(True)
    newax.patch.set_visible(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.spines['bottom'].set_position(('outward', shift_down))
    newax.tick_params(axis='both', which='major')#, labelsize=16)
    newax.grid('off')
    return newax