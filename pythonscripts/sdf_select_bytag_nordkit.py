#!/bin/python
"""Jean-Remy Marchand 24th Oct. 2016
This script aims at extracting in a SDF file
generated by chemaxon only the tautomer with an
estimated occupancy > 24.9%
and renames molecules with 2+ tautomers
"""

import sys, os

if len(sys.argv)!=3:
	print "Usage = python script.py file.sdf output.sdf"

a=open(sys.argv[1], "r")
ori=a.readlines()
a.close()

out=""
tmp=""
flag=False
shouldI=False

c=0
name=""
names=[]

for line in ori:
	# Naming the tautomer
	c+=1
	if c==1:
		d=2
		if len(line) > 0 and line.split()[0]+"_tauto_1" in names:
			while line.split()[0]+"_tauto_"+str(d) in names:
				d+=1
			name=line.split()[0]+"_tauto_"+str(d)
		if len(line) > 0 and line.split()[0]+"_tauto_1" not in names:
			name=line.split()[0]+"_tauto_1"
		names.append(name)
		tmp+=name
		tmp+="\n"
	# Checking the occupancy
	if flag==True:
		if float(line.split()[0])>24.9:
			shouldI=True
		flag=False
		tmp+=line
	elif "tauto_occupancy" in line:
		flag=True
		tmp+=line
	elif c!=1:
		tmp+=line
	if "$$$" in line:
		if shouldI==True:
			out+=tmp
			shouldI=False
		tmp=""
		c=0
b=open(sys.argv[2], "w")
b.write(out)
b.close()
