# -*- coding: utf-8 -*-
import matplotlib as mpl

dark_grey = '#444444'
almost_black = '#222222'
barb_incr = {'half': 2.5, 'full': 5.0, 'flag': 25.0}

# Line plot styles
cvec = ([[0.25, 0.25, 0.25],
         [0.90, 0.10, 0.10],
         [0.90, 0.90, 0.10],
         [0.25, 0.75, 0.25],
         [1.00, 0.50, 0.00],
         [0.25, 0.75, 1.00],
         [0.25, 0.75, 0.25],
         [1.00, 0.00, 0.25],
         [0.75, 0.25, 0.25],
         [1.00, 0.75, 0.25],
         [0.25, 1.00, 0.00],
         [0.50, 0.25, 0.75],
         [0.50, 0.75, 0.75],
         [0.75, 0.25, 0.25],
         [0.25, 1.00, 0.25]])


def div_cmap(numcolors=11, name='bwr_div_cmap',
             mincol='blue', midcol='white', maxcol='red',
             under=None, midcol_alpha=0, over=None):
    """ Create a custom diverging colormap with three colors

    Default is blue to transparent white to red with 11 colors.
    Colors can be specified in any way understandable
    by matplotlib.colors.ColorConverter.to_rgb()
    """
    c_max = mpl.colors.colorConverter.to_rgba(maxcol)
    c_min = mpl.colors.colorConverter.to_rgba(mincol)
    c_mid = mpl.colors.colorConverter.to_rgba(midcol, alpha=midcol_alpha)
    colors = [c_min, c_mid, c_max]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(name=name,
                                                        colors=colors,
                                                        N=numcolors)

    if under is not None:
        cmap.set_under(under)
    if over is not None:
        cmap.set_over(over)

    return cmap
