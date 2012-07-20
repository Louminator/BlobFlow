# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 09:09:13 2012

@author: Louis Rossi
"""

from pylab import *
from string import split,atof,atoi
from matplotlib.patches import Patch,Ellipse

def func3(x,y):
    return (1- x/2 + x**5 + y**3)*exp(-x**2-y**2)
#    return(x**2+y**2/4)

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
    f = open('expA0010.grd','r')
    txt = f.read()
    f.close()
    ws = split(txt)
    w = map(atof,ws)
    w = reshape(w,(gridn,gridn))
except IOError:
    print "No vtx file found."

try:
    f = open('expA0010.vtx','r')
    txt = f.read()
    f.close()
    vs = split(txt)
    vdata = map(atof,vs)
    vdata = reshape(vdata,(len(vdata)/6,6))
except IOError:
    print "No vtx file found."

num_vorts = len(vdata)


    
x = linspace(num[0],num[2],gridn)
y = linspace(num[1],num[3],gridn)

X,Y = meshgrid(x, y)

Z = func3(X, Y)
#pcolor(X, Y, Z)
pcolor(X,Y,w,edgecolors=None,alpha=0.5)
colorbar()
axis([-num[0],num[1],-num[2],num[3]])

ax = gca()

for k in range(0,num_vorts):
    circ = Ellipse((vdata[k,0],vdata[k,1]), \
    (vdata[k,3]*vdata[k,4])**0.5, \
    (vdata[k,3]/vdata[k,4]/2)**0.5, \
    angle = vdata[k,5]*180/pi, \
    fill=False)

    ax.add_patch(circ)

show()