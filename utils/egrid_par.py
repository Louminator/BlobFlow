#!/usr/bin/python

from Tkinter import *
from string import split,atof,atoi
import re
from os import system
from sys import argv

class egrid_front:
    
    def __init__(self,master):

	frame = Frame(master)

	self.filename_label = Label(frame,text="Base file name:")
	self.filename_label.grid(row=0,column=0,sticky=W)

        self.filename_entry = StringVar()
	self.filename_entry = Entry(frame)
	self.filename_entry.grid(row=0,column=1)

	self.ll_label = Label(frame,text="Lower left:")
	self.ll_label.grid(row=1,column=0,sticky=W)

	self.ll_entry = Entry(frame)
	self.ll_entry.grid(row=1,column=1)

	self.ur_label = Label(frame,text="Upper right:")
	self.ur_label.grid(row=2,column=0,sticky=W)

	self.ur_entry = Entry(frame)
	self.ur_entry.grid(row=2,column=1)

	self.n_label = Label(frame,text="Grid size:")
	self.n_label.grid(row=3,column=0,sticky=W)

	self.n_entry = Entry(frame)
	self.n_entry.grid(row=3,column=1)

	self.cores_label = Label(frame,text="Cores:")
	self.cores_label.grid(row=3,column=2,sticky=W)

	self.cores_entry = Entry(frame)
	self.cores_entry.grid(row=3,column=3)
        self.cores_entry.insert(END,"1")

	self.start_label = Label(frame,text="Start index:")
	self.start_label.grid(row=0,column=2,sticky=W)

	self.start_entry = Entry(frame)
	self.start_entry.grid(row=0,column=3)

	self.stop_label = Label(frame,text="Stop index:")
	self.stop_label.grid(row=1,column=2,sticky=W)

	self.stop_entry = Entry(frame)
	self.stop_entry.grid(row=1,column=3)

	self.stride_label = Label(frame,text="Stride:")
	self.stride_label.grid(row=2,column=2,sticky=W)

	self.stride_entry = Entry(frame)
	self.stride_entry.grid(row=2,column=3)

	self.quiet = StringVar()
	self.quiet.set("")
	self.cb = Checkbutton(frame,text="quiet",variable=self.quiet, \
	onvalue="-q",offvalue="")
	self.cb.grid(row=4,column=0,sticky=W)

        self.xanti = StringVar()
        self.xanti.set("")
	self.cb = Checkbutton(frame,text="X-antisymmetry", \
	variable=self.xanti,onvalue="-x",offvalue="")
	self.cb.grid(row=4,column=1,sticky=W)

	self.button = Button(frame,text="Project",fg="green", \
	command=self.project)
	self.button.grid(row=5,column=0)

	self.button = Button(frame,text="Quit",fg="red",command=frame.quit)
	self.button.grid(row=5,column=1)

	try:
	    f = open('egrid.default','r')
	    txt = f.read()
	    f.close()
	    nums = split(txt)
	    self.ll_entry.insert(END,"("+nums[0]+","+nums[1]+")");
	    self.ur_entry.insert(END,"("+nums[2]+","+nums[3]+")");
	    self.n_entry.insert(END,nums[4]);
	except IOError:
	    print "No egrid.default file found."

        self.start_entry.insert(END,"0");
        self.stop_entry.insert(END,"0");
        self.stride_entry.insert(END,"1");
        if (len(argv) == 2):
            self.filename_entry.insert(END,argv[1])


	frame.pack()

    def project(self):
	tmp = self.ll_entry.get()

        patternA = re.compile('\([-,+]*\d*\.*\d*,')
        patternB = re.compile(',[-,+]*\d*\.*\d*\)')
        searchA  = patternA.search(tmp)
        searchB  = patternB.search(tmp)

        x0 = atof(tmp[searchA.start()+1:searchA.end()-1])
	y0 = atof(tmp[searchB.start()+1:searchB.end()-1])

	tmp = self.ur_entry.get()
        searchA  = patternA.search(tmp)
        searchB  = patternB.search(tmp)
	x1 = atof(tmp[searchA.start()+1:searchA.end()-1])
	y1 = atof(tmp[searchB.start()+1:searchB.end()-1])

	n = atoi(self.n_entry.get())

	f = open('egrid.default','w')

	f.write(str(x0) + " " + str(y0)+"\n")
	f.write(str(x1) + " " + str(y1)+"\n")
	f.write(str(n)+"\n")
	f.close()

	start = atoi(self.start_entry.get())
	stop = atoi(self.stop_entry.get())
	stride = atoi(self.stride_entry.get())
        
	for num in range(start,stop+1,stride):
	    numstr = "%04d" % num
	    cmd = "mpiexec -np " + self.cores_entry.get() + " egrid_par " + \
            self.xanti.get() + " " + self.quiet.get() + " " + \
	    self.filename_entry.get() + numstr
#	    print cmd
	    try:
		f = open(self.filename_entry.get()+numstr+".vtx","r")
		f.close()
		system(cmd)
	    except IOError:
		print "File not present: "+ \
		self.filename_entry.get()+numstr+".vtx"

root = Tk()

root.title("egrid")
app = egrid_front(root)

root.mainloop()
