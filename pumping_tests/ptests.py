#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 21:16:38 2018

@author: Theo
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mlu.mlu_xml import mlu2xml, Mluobj
from mlu.ptest import Pumptest, MLU_ptest
from fdm.mfgrid import Grid

datapath = './data'

ptestName  = 'gat_boomse_klei'
#ptestName  = 'zuidelijke_landtong'

mlu_file = os.path.join(datapath, ptestName + '.mlu')
xml_file = os.path.join(datapath, ptestName + '.xml')


#%% Convert mlu file to xmlfile
mlu2xml(mlu_file, fout=xml_file)


#%% Generate a Mluobj and plot the drawdowns (contained in the xml_file)

mu = Mluobj(xml_file)

mu.plot_drawdown(mu.obsNames, yscale='linear', xscale='log', marker='.')

mu.plot_section(xscale='log')
fig = plt.gcf()
fig.set_size_inches(4, 8)

#mu.toshape(ptestName)


mu.plot_drawdown('PB11-6', yscale='linear', xscale='log', ylim=(-1,4), marker='.')


'''
#%% simulate pumping test contained in mlu file

mptest = MLU_ptest(xml_file, Rw=0.13, tshift=0.)


#%% Simulate ptest by setting up the simulatoin from ground using the
# parameters in the MLU file.


# Firstly define the grid
Rw   = 0.13    # m, well borehole radius
Rmax = 10000. # m, extent of model (should be large enough)

r  = np.hstack((0, np.logspace(np.log10(Rw), np.log10(Rmax),
                               int(np.ceil(10 * np.log10(Rmax/Rw) + 1)))))

dz = np.array([9., 16., 0.01, 3.79, 3.7, 1., 3.7, 1., 3.8, 1., 1.3, 5., 3.4, 6.3])
z  = np.hstack((0, -np.cumsum(dz)))
gr = Grid(r, [-0.5, 0.5], z, axial=True)

# Set soil parameters
kD = np.array([48, 11.37, 0.13, 0.15, 0.16, 0.75, 31.5])
c  = np.array([1500, 1.e-2, 2.85, 12.3, 76., 26., 56.7])
Sat= np.array([0., 0., 0., 0., 0., 0., 0.,])
Saq =np.array([1.e-04, 1.e-4, 1.e-5, 1.e-5, 1.e-5, 1.e-5, 1.e-5])
top_aquitard = True  # whether the top layer is a top aquitard
Q = np.array([0., 0., 0., 0., 0., 40.77, 0.])       # extraction from each aquifer

# simulation time
t = np.logspace(-5, 1, 101) # times for simulation
# specify observation points [(anme, r, iaquifer), ...]
well = ('pp4'   , 45619.637, 372414.25, 5)

points = [('pb11_3', 45620.875, 372422.43, 2),
          ('pb11_4', 45622.344, 372423.27, 3),
          ('pb11-5', 45621.608, 372421.04, 4),
          ('pb11-6', 45619.263, 372421.36, 5),
          ('pp4'   , 45619.637, 372414.25, 5)]

# Generate points as [(name, r, z), (...), ...]
# Distance to well
xy = np.array([(p[1], p[2]) for p in points])
r_ow= np.sqrt((well[1] - xy[:,0])**2 + (well[2] - xy[:,1])**2)
# Make sure min distance equals well radius
r_ow = np.fmax(r_ow, Rw)

# Convert (zero-based) aquifer numbers to z of aquifer center
Iaq = np.array([p[3] for p in points], dtype=int) # aquifer number
Ilay = Iaq * 2 + 1 if top_aquitard else Iaq * 2   # fdm layer number
z_ow = (z[Ilay] + z[Ilay + 1]) / 2.                  # z of layer center

# generate obswells for ptest
obswells = [(p[0], r, zm_) for p, r, zm_ in zip(points, r_ow, z_ow)]


# Get pumping test simulation object, it runs immediately
pt = Pumptest( gr=gr, kD=kD, c=c, S=(Saq, Sat), top_aquitard=True,
               t=t, Q=Q, obswells=obswells)

# show the drawdown curves for the observation points
pt.show(xscale='log', grid=True)

'''