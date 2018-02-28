#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 23:25:00 2018

@author: Theo
"""
import numpy as np
from fdm.fdm3t import Fdm3t
from fdm.mfgrid import Grid
import pandas as pd
import matplotlib.patches as patches
import matplotlib.pyplot as plt

AND = np.logical_and



def get_params(workbook=None, sheetname=None, header=0):
    '''Return model parameters from workbook/sheetname.

    parameters
    ----------
        workbook : Excel workbook
            workbook with parameters
        sheetname : sheetname in excel workbook
            sheet with the table with parameters.

            Columns required are:
                name, ztop, zbot, D, nsub, kh, kv, ss
        header : int
            zero-based rownumber in sheet with the column headers

    '''


    cols=['name', 'ztop', 'zbot', 'D', 'nsub', 'kh', 'kv', 'ss', 'color', 'alpha']
    layers = pd.read_excel(workbook, sheetname='model', header=4)
    layers = layers[cols]

    ilay = 0
    layer = dict()
    for i in range(len(layers)):
        layer[ilay] = {col :layers[col].iloc[i] for col in layers.columns}
        nsub = layer[ilay]['nsub']
        if nsub > 1:
            layer[ilay]['D'] = layer[ilay]['D'] / nsub
            layer[ilay]['zbot'] = layer[ilay]['ztop'] - layer[ilay]['D']
            layer[ilay]['nsub'] = 1
            ilay += 1
            for j in range(nsub - 1):
                layer[ilay] = layer[ilay-1].copy()
                layer[ilay]['ztop'] = layer[ilay-1]['zbot']
                layer[ilay]['zbot'] = layer[ilay]['ztop'] - layer[ilay]['D']
                ilay += 1
        else:
            ilay += 1
    return layer


#if __name__ == '__main__':


#%% Laagopbouw

iv = 9 # variant Nr


workbook = 'paramsGatBoomseKlei.xls'

layer = get_params(workbook, 'model', header=4)

variant = dict(pd.read_excel(workbook, sheetname='varianten', header=2, index_col='variant').T)


#if iv < 9:
#    keys = set(layer.keys())
#    for key in layer:
#        if 0.5 * (layer[key]['ztop'] + layer[key]['zbot']) < -59:
#            keys.remove(key)
#    layer = {k:layer[k] for k in keys}

kwand  = 1e-10  # m/d
rwand  = (64.95, 65.05)
zBklei = (-22.0, -36.0)

# choosing which variant to simulate

hwell = -11.2

zwand = (0., variant[iv]['zwand'])
rwell = (variant[iv]['rwell1'], variant[iv]['rwell2'])
zwell = (variant[iv]['zwell1'], variant[iv]['zwell2'])



#%% Grid

z = [layer[0]['ztop']]
kh = []
kv = []
ss = []
for i in range(len(layer)):
    z.append(layer[i]['zbot'])
    kh.append(layer[i]['kh'])
    kv.append(layer[i]['kv'])
    ss.append(layer[i]['ss'])

y  = [ -0.5, 0.5]
r  = np.hstack((0, *rwell, *rwand, np.logspace(0, 4, 101)))
gr = Grid(r, y, z, axial=True)


wand = np.zeros(gr.shape, dtype=bool)
wand[AND( AND( AND(gr.XM > rwand[0], gr.XM < rwand[1]),
                    gr.ZM < zwand[0]),
                    gr.ZM > zwand[1])] = True

well = np.zeros(gr.shape, dtype=bool)
well[AND( AND( AND(gr.XM  > rwell[0], gr.XM < rwell[1]),
                    gr.ZM < zwell[0]),
                    gr.ZM > zwell[1])] = True

vloer = np.zeros(gr.shape, dtype=bool)
vloer[AND(gr.XM < rwand[0], gr.ZM > -21.9)] = True

#%% parameters
Kh = gr.const(np.array(kh));  Kh[wand] = kwand
Kv = gr.const(np.array(kv))
Ss = gr.const(np.array(ss))

#%% Boundary array and boundary conditions
IBOUND = gr.const(1, dtype=int)
IBOUND[0]    = -1
IBOUND[well] = -1
IBOUND[vloer] = 0

FQ = gr.const(0.)
HI = gr.const(0.)
HI[well] = hwell

#%% simulation time
t = np.logspace(-3, 3, 101) # times for simulation

# %%specify observation points [(anme, r, iaquifer), ...]
points = [('well',    0.5*(rwell[0] + rwell[1]), -40.),
          ('TT-001',    1., -40.),
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

rw = np.round(0.5 * (rwell[0] + rwell[1]), 1)
zw = np.round(0.5 * (zwell[0] + zwell[1]), 1)

points = [('well_{:.0f}'.format(rw),    rw, zw),
          ('P0001_40',    1., -40.),
          ('P0001_50',    1., -50.),
          ('P0001_70',    1., -70.),
          ('p0075_50',   75., -50.),
          ('p0250_50',  250., -50.),
          ('p0500_50',  500., -50.),
          ('p1000_50', 1000., -50.),
          ('p5000_50', 5000., -50.)]

fdmObj = Fdm3t(gr=gr, t=t, kxyz=(Kh, Kh, Kv), Ss=Ss, FQ=FQ, HI=HI,
            IBOUND=IBOUND, epsilon=1.0)


Q = fdmObj.out['Q'][-1][:, 0, :]

Qd = np.sum(Q[1:, :])
Qh = Qd/24.
print('Extraction flow Q = {:.2f} m3/h', Qh)


fdmObj.show(points, xscale='log', size_inches=(11.0, 3.1),
            title='Sluishoofd bemaling, variant {}, Q={:.1f} m3/h'.format(iv, Qh))

p = list()
for iL in layer:
    p.append(patches.Rectangle((  0., layer[iL]['zbot']),\
                    width=200., height=layer[iL]['ztop'] - layer[iL]['zbot'],\
                    color=layer[iL]['color'], alpha=layer[iL]['alpha']))

# Boomse Klei
p.append(patches.Rectangle((  0., zBklei[1]),\
                width=rwand[0], height=(zBklei[0] - zBklei[1]), color='darkgreen', zorder=100))

# Vloer
p.append(patches.Rectangle((  0., zBklei[0]),\
                width=rwand[0], height=2., color='gray', zorder=100))

# Openruimte in sluishoofd
p.append(patches.Rectangle((  0., zBklei[0] + 2.0),\
                width= rwand[0] , height=zwand[0] - zBklei[0] - 2.0, color='white', zorder=101))

# Diepwand
p.append(patches.Rectangle((rwand[0], zwand[1]),\
                width=1.0 , height=zwand[0] - zwand[1], color='black', zorder=102))

# Onttrekking
p.append(patches.Rectangle((rwell[0], zwell[1]),\
                width=rwell[1] - rwell[0], height=zwell[0] - zwell[1], color='white', zorder=102))


#ax = fdmObj.contour(dphi=1.0, dpsi=100., xlim=(0, 200.), xscale='linear', patches=p)
ax = fdmObj.contour(dphi=1.0, dpsi=100., xlim=(0, 200.), xscale='linear',
                    patches=p, colors='r', linestyles='-', size_inches=(11, 4.5),
                    title='Stijghoogten en stroomlijnen, variant {}, Q={:.1f} m3/h'.format(iv, Qh))

print('{:2s} {:26s} {:8s} {:8s} {:8s} {:8s} {:8s} {:8s}'\
          .format('nr', 'name', 'ztop', 'zbot', 'D',
                      'kh', 'kv', 'ss'))

for ilay in range(len(layer)):
    lay = layer[ilay]
    print('{:2d} {:22s} {:8.2f} {:8.2f} {:8.2f} {:8.3g} {:8.3g} {:8.3g}'\
          .format(ilay,
                  lay['name'], lay['ztop'], lay['zbot'], lay['D'],
                      lay['kh'], lay['kv'], lay['ss']))


#%% Vertaling onttrekking naar individuele putten

n = 16
r = 65.
x = [r * np.cos(i/n * 2 * np.pi) for i in range(n)]
y = [r * np.sin(i/n * 2 * np.pi) for i in range(n)]

rw= 0.25

x0 = r
y0 = rw
x1 = r * np.cos(0.5/n * 2 * np.pi)
y1 = r * np.sin(0.5/n * 2 * np.pi)

Q = 60 * 24 / n
kD = 72.
R0 = 1000.

s0 = np.sum([-Q/(2 * np.pi * kD) * np.log(np.sqrt((xx - x0)**2 + (yy-y0)**2) / R0) for xx, yy in zip(x, y)])
s1 = np.sum([-Q/(2 * np.pi * kD) * np.log(np.sqrt((xx - x1)**2 + (yy-y1)**2) / R0) for xx, yy in zip(x, y)])

print(s0, s1)

Q* n/ (2 * np.pi * kD) * np.log(65/R0)

#%% HK

c = 7000 # d