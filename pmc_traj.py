# -*- coding: utf-8 -*-
"""
Collection of some useful functions to work with
text file containing polar lows trajectories
"""
# Standard imports
import datetime
import matplotlib.pyplot as plt
import numpy as np
# Local imports
from phys_meteo import r_earth
# TODO: docstrings!


class pmc:
    def __init__(self, n, dat):
        self.n = n
        self.dt = [datetime.datetime.strptime(str(int(i)), '%Y%m%d%H%M')
                   for i in dat[dat[:, 0] == n, 1]]
        self.lon = [i for i in dat[dat[:, 0] == n, 2]]
        self.lat = [i for i in dat[dat[:, 0] == n, 3]]
        self.dist = []
        self.vel = []
        if len(self.lon) > 1:
            self.calc_dist(r_earth)
            self.calc_vel()
            self.calc_met_dir()
            self.calc_vel_uv()

    def calc_dist(self, r):
        lon = np.radians(self.lon)
        lat = np.radians(self.lat)
        self.dist = []
        for i in range(len(lon)-1):
            if lat[i+1] == lat[i] and lon[i+1] == lon[i]:
                self.dist.append(0.)
            else:
                self.dist.append(np.arccos(np.sin(lat[i]) * np.sin(lat[i+1])
                                 + np.cos(lat[i]) * np.cos(lat[i+1])
                                 * np.cos(lon[i] - lon[i+1])))
        self.dist = r*np.array(self.dist)

    def calc_vel(self):
        self.mean_vel = self.dist.sum() / (self.dt[-1]
                                           - self.dt[0]).total_seconds()
        self.vel = []
        for i in range(len(self.dt)-1):
            if self.dt[i] == self.dt[i+1]:
                self.vel.append(None)
            else:
                self.vel.append(self.dist[i] / (self.dt[i+1]
                                                - self.dt[i]).total_seconds())
        self.vel = np.array(self.vel)

    def calc_met_dir(self):
        self.mean_prop_dir = 180 + 180./np.pi*np.arctan2(self.lon[-1]
                                                         - self.lon[0],
                                                         self.lat[-1]
                                                         - self.lat[0])
        self.prop_dir = []
        for i in range(len(self.dt)-1):
            if self.dt[i] == self.dt[i+1]:
                self.prop_dir.append(None)
            else:
                dlon = self.lon[i+1] - self.lon[i]
                dlat = self.lat[i+1] - self.lat[i]
                self.prop_dir.append(180 + 180./np.pi*np.arctan2(dlon, dlat))

    def calc_vel_uv(self):
        assert hasattr(self, 'vel'), self.calc_vel()
        assert hasattr(self, 'prop_dir'), self.calc_met_dir()

        self.prop_u = []
        self.prop_v = []
        for i in range(len(self.dt)-1):
            self.prop_u.append(-self.vel[i]
                               * np.sin(np.radians(self.prop_dir[i])))
            self.prop_v.append(-self.vel[i]
                               * np.cos(np.radians(self.prop_dir[i])))

    def expand_last(self):
        if len(self.dist) < len(self.dt):
            self.dist += [self.dist[-1]]
        if len(self.vel) < len(self.dt):
            self.vel += [self.vel[-1]]
        if len(self.prop_dir) < len(self.dt):
            self.prop_dir += [self.prop_dir[-1]]
        if len(self.prop_u) < len(self.dt):
            self.prop_u += [self.prop_u[-1]]
        if len(self.prop_v) < len(self.dt):
            self.prop_v += [self.prop_v[-1]]


def get_pmc_list(fname):
    dat = np.loadtxt(fname, delimiter="\t")
    ids = [int(i) for i in set(dat[:, 0])]
    pmc_list = []
    for j in ids:
        pmc_list.append(pmc(j, dat))
    return pmc_list


def plot_pmc_traj_bw(m, pmc_list):
    for i in pmc_list:
        x, y = m(i.lon, i.lat)
        lnw = 1.5
        if i.n == 12:
            clr = 'r'
            lnst = '-'
        else:
            clr = 'k'
            lnst = '--'
        traj = m.plot(x, y,
                      color=clr, linestyle=lnst, linewidth=lnw,
                      label='pmc_'+str(i.n))
        beg = m.plot(x[0], y[0], marker='o', mfc=clr, mec=clr, ms=10)
        fin = m.plot(x[-1], y[-1], marker='o', mfc='w', mec=clr, ms=10)
        for ii, iart in enumerate([traj, beg, fin]):
            iart[0].set_zorder(100+ii)


def plot_pmc_traj_great(m, pmc_list):
    for i in pmc_list:
        x, y = m(i.lon, i.lat)
        lnw = 1.5
        if i.n == 12:
            clr = 'm'
            lnst = '-'
        else:
            clr = 'b'
            lnst = '--'
        for j in range(len(i.lon)-1):
            traj = m.drawgreatcircle(i.lon[j], i.lat[j],
                                     i.lon[j+1], i.lat[j+1],
                                     color=clr, linestyle=lnst, linewidth=lnw,
                                     label='pmc_'+str(i.n))
            traj[0].set_zorder(100)
        beg = m.plot(x[0], y[0], marker='o', mfc=clr, mec=clr, ms=10)
        fin = m.plot(x[-1], y[-1], marker='o', mfc='w', mec=clr, ms=10)
        for ii, iart in enumerate([beg, fin]):
            iart[0].set_zorder(101+ii)


if __name__ == '__main__':
    import sys
    import mypaths
    fmt = 'pdf'
    svfigkw = dict(format=fmt, dpi=300, bbox_inches='tight')

    fname = mypaths.trajf
    pmc_list = get_pmc_list(fname)

    fig, ax = plt.subplots(figsize=(25, 15))

    try:
        import map_plot_func as mymap
        m = mymap.make_map(ax=ax)
    except ImportError:
        try:
            from mpl_toolkits.basemap import Basemap
            lon1 = -5
            lon2 = 50.
            lat1 = 65.
            lat2 = 80.
            m = Basemap(ax=ax, projection='lcc',
                        llcrnrlon=lon1, llcrnrlat=lat1,
                        urcrnrlon=lon2, urcrnrlat=lat2,
                        lat_1=lat1+2./3*(lat2-lat1), lat_2=lat2,
                        lon_0=0.5*(lon1+lon2),
                        resolution='c')
            m.drawcoastlines(color='0.5')
            m.fillcontinents(color='0.75')

            plot_pmc_traj_bw(m, pmc_list)
            fig.savefig('pmc_traj_accacia.{fmt}'.format(fmt=fmt), **svfigkw)
        except:
            print('Something went wrong...')
            sys.exit(1)
