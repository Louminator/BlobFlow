#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 09:09:13 2012

@author: Louis Rossi
"""

from scipy import *
from sys import argv
from pylab import *
from string import split,atof,atoi
from matplotlib.patches import Patch,Ellipse

if (len(argv) < 2):
    print "Usage: vtxplot.py <filename> <list of frames>"

frameList = []
vtx = False
grd = True
argInd = 1
while (argInd < len(argv)):
    if (argv[argInd][0] == '-'):
        if argv[argInd][1:] == 'v':
            vtx = True
        elif argv[argInd][1:] == 'nv':
            vtx = False
        elif argv[argInd][1:] == 'g':
            grd = True
        elif argv[argInd][1:] == 'ng':
            grd = False
    else:
        try:
            frameList.append(atoi(argv[argInd]))
        except ValueError:
            name = argv[argInd]
    argInd += 1
    
#name = '/home/rossi/Research/Oseen-explorations/expB/expB'
#n = 0


try:
    f = open('egrid.default','r')
    txt = f.read()
    f.close()
    nums = split(txt)
    num = map(atof,nums)
except IOError:
    print "No egrid.default file found."

gridn = int(num[4])

for n in frameList:
    
    figure(n)
    
    axis([num[0],num[2],num[1],num[3]])
    if grd:
        grdname = name+'{0:04d}'.format(n)+'.grd'
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
            
        if grdfile:
            x = linspace(num[0],num[2],gridn)
            y = linspace(num[1],num[3],gridn)
            
            X,Y = meshgrid(x, y)
            
            pcolor(X,Y,w,edgecolors='None',shading='faceted')
            colorbar()
    
    if vtx:
        vtxname = name+'{0:04d}'.format(n)+'.vtx'
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