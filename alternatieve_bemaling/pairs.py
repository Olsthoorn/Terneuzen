#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:20:25 2018

@author: Theo
"""

import numpy as np
import matplotlib.pyplot as plt
from fdm.mfgrid import Grid

AND = np.logical_and
OR  = np.logical_or
NOT = np.logical_not

# find points to the leftof a line or a polygon
def cell_pairs(gr, polygon):
    '''return cell pairs left and right of polygon contour like for HB package.

    parameters
    ----------
        gr: fdm.mfgrid.Grid
            grid object
        polygon: list of coordinate pairs or a XY array [n, 2]
            contour coordinates
    '''


    In = gr.inpoly(polygon)
    mask_west  = np.hstack((In, np.zeros((gr.ny, 1), dtype=bool)))
    mask_south = np.vstack((np.zeros((1, gr.nx), dtype=bool), In))

    west  = np.abs(np.diff(mask_west , axis=1)) == 1
    south = np.abs(np.diff(mask_south, axis=0)) == 1

    east = np.hstack((np.zeros((gr.ny, 1), dtype=bool), west[:,:-1]))
    north= np.vstack((south[1:,:], np.zeros((1, gr.nx), dtype=bool)))

    pairs = np.array([np.hstack((gr.NOD[0][west], gr.NOD[0][north])),
                      np.hstack((gr.NOD[0][east], gr.NOD[0][south]))]).T

    return pairs

if __name__ == '__main__':

    x = np.linspace(-100., 100., 21)
    y = np.linspace(-100., 100., 21)
    z = [0, -10, -20]

    gr = Grid(x, y, z)

    polygon = np.array([(23, 15), (45, 50.), (10., 81.), (-5., 78), (-61., 51.), (-31., 11.),
               (-6., -4.), (-42., -20.), (-50., -63.), (-7., -95.),
               (31., -80.), (60., -71.), (81., -31.), (5., -63.), (25., -15.), (95., 40.),
               (23, 15)])

    pairs = cell_pairs(gr, polygon)


    fig, ax = plt.subplots()
    ax.set_title('Node pairs for the hor. flow-barrier package of Modflow')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    gr.plot_grid(world=False, ax=ax)
    ax.plot(polygon[:,0], polygon[:,1])
    ax.plot(gr.Xm.ravel()[pairs[:,0]], gr.Ym.ravel()[pairs[:,0]], '.r', label='column 1')
    ax.plot(gr.Xm.ravel()[pairs[:,1]], gr.Ym.ravel()[pairs[:,1]], '.b', label='column 2')
    for pair in pairs:
        ax.plot(gr.Xm.ravel()[pair], gr.Ym.ravel()[pair], 'k-')
