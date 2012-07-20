#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 09:09:13 2012

@author: Louis Rossi
"""

from pylab import *
from string import split,atof,atoi
from matplotlib.patches import Patch,Ellipse

name = 'expD'
n = 5

vtxname = name+'{0:04d}'.format(n)+'.vtx'
grdname = name+'{0:04d}'.format(n)+'.grd'

try:
    f = open('egrid.default','r')
    txt = f.read()
    f.close()
    nums = split(txt)
    num = map(atof,nums)
except IOError:
    print "No egrid.default file found."
# make these smaller to increase the resolution

gridn = int(num[4])

try:
    f = open(grdname,'r')
    txt = f.read()
    f.close()
    ws = split(txt)
    w = map(atof,ws)
    w = reshape(w,(gridn,gridn))
    grdfile = True
except IOError:
    print "No grd file found."
    grdfile = False

try:
    f = open(vtxname,'r')
    txt = f.read()
    f.close()
    vs = split(txt)
    vdata = map(atof,vs)
    vdata = reshape(vdata,(len(vdata)/6,6))
    vtxfile = True
except IOError:
    print "No vtx file found."
    vtxfile = False

num_vorts = len(vdata)

axis([-num[0],num[1],-num[2],num[3]])
if grdfile:
    x = linspace(num[0],num[2],gridn)
    y = linspace(num[1],num[3],gridn)
    
    X,Y = meshgrid(x, y)
    
    pcolor(X,Y,w,edgecolors='None',shading='faceted')
    colorbar()

if vtxfile:
    ax = gca()
    
    for k in range(0,num_vorts):
        circ = Ellipse((vdata[k,0],vdata[k,1]), \
        (vdata[k,3]*vdata[k,4])**0.5, \
        (vdata[k,3]/vdata[k,4])**0.5, \
        angle = vdata[k,5]*180/pi, \
        fill=False,ls='dotted')
    
        ax.add_patch(circ)

show()