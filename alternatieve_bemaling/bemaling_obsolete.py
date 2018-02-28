#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 21:16:38 2018

@author: Theo
"""

import numpy as np
import matplotlib.pyplot as plt
from fdm.mfgrid import Grid

from fdm.fdm3t import fdm3t
from showdd import show, contour

AND = np.logical_and

#%% Laagopbouw


# set soil parameters
kD =  np.array([50., 1.5, 7.5, 20])
c  =  np.array([100., 2000., 2, 10])
Saq = np.array([0.001, 0.001,0.001, 0.001])
Sat = np.array([0., 0., 0., 0.])


top_aquitad = True

kwand = 1e-3  # m/d
rwand = (64.95, 65.05)
zwand = (0., -45.)
rwell = (rwand[0] - 0.5, rwand[0])
zwell = (-36., -42. )
hwell = -20.
#      c1  kD1   c2    kD2   c3    kD3   c4  kD4
z = [0, -2., -22., -36., -42, -42.5, -45, -50, -60]
y = [ -0.5, 0.5]
r = np.hstack((0, *rwell, rwand[0], rwand[1], *rwand, np.logspace(0, 4, 101)))

gr = Grid(r, y, z, axial=True)

wand = np.zeros(gr.shape, dtype=bool)
wand[AND( AND( AND(gr.XM > rwand[0], gr.XM < rwand[1]),
                    gr.ZM < zwand[0]),
                    gr.ZM > zwand[1])] = True

well = np.zeros(gr.shape, dtype=bool)
well[AND( AND( AND(gr.XM > rwell[0], gr.XM < rwell[1]),
                    gr.ZM < zwell[0]),
                    gr.ZM > zwell[1])] = True

vloer = np.zeros(gr.shape, dtype=bool)
vloer[AND(gr.XM < rwand[0], gr.ZM > -22)] = True


iaquif = np.zeros(gr.nz, dtype=bool)
iaquif[1::2] = True
iatard = np.logical_not(iaquif)

kh = np.zeros(gr.nz)
kv = np.zeros(gr.nz)
ss = np.zeros(gr.nz)

kh[iaquif] = kD / gr.dz[iaquif]
kv[iatard] = gr.dz[iatard] /c
kv[iaquif] = kh[iaquif] / 3
kh[iatard] = kv[iatard] * 3; kv[0] /= 2. # top layer
ss[iaquif] = Saq / gr.dz[iaquif]
ss[iatard] = Sat / gr.dz[iatard]

Kh = gr.const(kh);  Kh[wand] = kwand
Kv = gr.const(kv)
Ss = gr.const(ss)


IBOUND = gr.const(1, dtype=int)
IBOUND[0]    = -1
IBOUND[well] = -1
IBOUND[vloer] = 0

FQ = gr.const(0.)
HI = gr.const(0.)
HI[well] = hwell

# simulation time
t = np.logspace(-3, 3, 101) # times for simulation

# specify observation points [(anme, r, iaquifer), ...]
points = [('well',    0.5*(rwell[0] + rwell[1]), -40.),
          ('TT-001 ',   1., -40.),
          ('TT-050' ,  50., -40.),
          ('TT-070' ,  70., -40.),
          ('TT-120' , 120., -40.),
          ('TT-175' , 175., -40.),
          ('TT-350' , 350., -40.),
          ('TT-700' , 700., -40.),
          ('TO-001',    1., -50.),
          ('TO-050',   50., -50.),
          ('TO-080',   70., -50.),
          ('TO-120',  120., -50.),
          ('TO-175',  175., -50.),
          ('TO-350',  350., -50.),
          ('TO-700',  700., -50.),
          ('BO-001',   1.,  -15.)]

points = [points[i] for i in [0, 1, 2, 8]]

out = fdm3t(gr=gr, t=t,
             kxyz=(Kh, Kh, Kv),
             Ss=Ss, FQ=FQ, HI=HI,
             IBOUND=IBOUND, epsilon=1.0)


show(gr, t, out, points, xscale='log', grid=True)


C = contour(gr, out, dphi=0.25, dpsi=10, xlim=(0, 150))
