# -*- coding: utf-8 -*-
"""
"""
# Standard libraries
import matplotlib as mpl
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
# Local libraries
import phys_meteo as met

# Default map boundaries
LON1 = -5.
LON2 = 40.
LAT1 = 65.
LAT2 = 78.

# Names on map
TOPONYMS = [dict(name='Svalbard', lon=16., lat=79.),
            dict(name='Norway', lon=23., lat=70.),
            dict(name='Norwegian Sea', lon=0., lat=69.),
            ]


def label_map(ax, bm, toponyms=TOPONYMS):
    pe_kw = dict(path_effects=[PathEffects.withStroke(linewidth=3,
                                                      foreground='k')])
    for i in toponyms:
        ix, iy = bm(i['lon'], i['lat'])
        txt = ax.annotate(i['name'], (ix, iy), size=20, ha='center', color='w',
                          **pe_kw)
        txt.set_zorder(200)


def make_map(ax=None, lon1=LON1, lon2=LON2, lat1=LAT1, lat2=LAT2,
             proj='lcc', resolution='l', tick_incr=[5, 1],
             fill=False, add_grid=True, coast=True, scale=False):
    if ax is None:
        ax = plt.gca()
    bm = Basemap(ax=ax,
                 llcrnrlon=lon1,
                 llcrnrlat=lat1,
                 urcrnrlon=lon2,
                 urcrnrlat=lat2,
                 lat_1=lat1+2./3*(lat2-lat1), lat_2=lat2,
                 lon_0=0.5*(lon1+lon2),
                 projection=proj, resolution=resolution)

    if add_grid:
        ticklon = np.array(tick_incr)[0]
        try:
            ticklat = np.array(tick_incr)[1]
        except IndexError:
            ticklat = ticklon

        meridians = np.arange(round(lon1)-30, lon2+30, ticklon)
        parallels = np.arange(round(lat1)-30, lat2+30, ticklat)
        bm.drawmeridians(meridians, labels=[0, 0, 0, 1],
                         color='0.5', size=14, latmax=90)
        bm.drawparallels(parallels, labels=[1, 0, 0, 0],
                         color='0.5', size=14)
    if coast:
        bm.drawcoastlines(color='0.5')
    if fill:
        bm.fillcontinents(color='0.75', alpha=0.25)
    if scale:
        ms = bm.drawmapscale(lon2-3.25, lat1+0.4, 0.5*(lon1+lon2),
                             lat1+2./3*(lat2-lat1), 100,
                             barstyle='fancy', fontsize=10, fontcolor='k',
                             zorder=200)
        for i in ms:
            if isinstance(i, mpl.text.Text):
                i.set_fontweight('bold')
                i.set_path_effects([PathEffects.withStroke(linewidth=1,
                                                           foreground='w')])

    return bm


def draw_rect_domain(bm, grd, units='km', color='b',
                     centre=True, circle=False, returnxy=False):
    """
    Draw a rectangular domain
    using central longitude and latitude, size and step in x- or y-direction

    **Arguments:**

    *bm*
        Basemap instance
    *grd*
        Grid parameters of the domain.
        Its type should be a dictionary with keywords:
        lon0, lat0, nx, ny, dx, dy

    **Optional arguments:**

    *units*
        "km" (default) or "deg"

    *color*
        Is passed to the matplotlib.pyplot.plot function. Default is "b".

    *centre*
        Mark centre of the domain. Default is True.

    *circle*
        Use Basemap.drawgreatcircle instead of plot function.

    *returnxy*
        Return map coordinates of the corners of the domain

    **Example:**
    Draw a domain defined by 400x400 points with the step of 0.02 degree
        >>> dom_2km = dict(lon0=20, lat0=70, dx=0.02, dy=0.02, nx=400, ny=400)
        >>> draw_rect_domain(bm, dom_2km, units='deg', color='r')
    """

    nx05 = grd['nx'] / 2.
    ny05 = grd['ny'] / 2.

    if units.lower() == 'km':
        m_in_deg = 2 * np.pi * met.r_earth / 360.
        llat = grd['lat0'] - ny05 * grd['dy'] * 1000. / m_in_deg
        ulat = grd['lat0'] + ny05 * grd['dy'] * 1000. / m_in_deg

        l_deg_len = m_in_deg * np.cos(np.radians(llat))
        u_deg_len = m_in_deg * np.cos(np.radians(ulat))
        llcrnrlon = grd['lon0'] - nx05 * grd['dx'] * 1000. / l_deg_len
        ulcrnrlon = grd['lon0'] - nx05 * grd['dx'] * 1000. / u_deg_len
        urcrnrlon = grd['lon0'] + nx05 * grd['dx'] * 1000. / u_deg_len
        lrcrnrlon = grd['lon0'] + nx05 * grd['dx'] * 1000. / l_deg_len

        lons = [llcrnrlon, ulcrnrlon, urcrnrlon, lrcrnrlon, llcrnrlon]
        lats = [llat,      ulat,      ulat,      llat,      llat]

    elif units.lower() == 'deg':
        llcrnrlon = grd['lon0'] - nx05 * grd['dx']
        urcrnrlon = grd['lon0'] + nx05 * grd['dx']
        llcrnrlat = grd['lat0'] - ny05 * grd['dy']
        urcrnrlat = grd['lat0'] + ny05 * grd['dy']

        lons = [llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon]
        lats = [llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat]

    else:
        raise ValueError('Keyword units must be km or deg')

    x, y = bm(lons, lats)

    if circle:
        if units.lower() == 'km':
            lons, lats = bm(x, y, inverse=True)
        for i, j, k, l in zip(lons[:-1], lats[:-1], lons[1:], lats[1:]):
            bm.drawgreatcircle(i, j, k, l,
                               marker='o', color=color, linewidth=2)
    else:
        bm.plot(x, y, marker='o', color=color, linewidth=2)

    if centre:
        x0, y0 = bm(grd['lon0'], grd['lat0'])
        bm.plot(x0, y0, marker='o', color=color, linewidth=2)

    if returnxy:
        return x[:-1], y[:-1]
