#!/usr/bin/python

from string import split,atof,atoi,atof

filein = open("split15asym.dat","r")

strs = split(filein.read())
filein.close()

fileout = open("split15new.dat","w")
data = map(atof, strs)

for i in range(0,len(data),7):
    fileout.write("{%14.7e, %14.7e, %14.7e, %14.7e, %14.7e, %14.7e}\n" % \
                  (data[i],data[i+1],data[i+2],data[i+3],data[i+4],data[i+5]))
fileout.close
