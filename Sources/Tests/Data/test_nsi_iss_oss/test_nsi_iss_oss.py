#! /usr/bin/python

import numpy as np
import os
import math
import sys

istat = 0

# Open file
fileid = 'output_test_nsi_iss_oss_stokes.txt'
f = open(fileid,'r')

# Read file
data = f.readlines()

# Close file
f.close()

# Extract data
for datum in data:
    datum = datum.split(":")
    key = datum[0]
    if(key==" Current memory usage"):
        memory_usage = float(datum[1].split("\n")[0])
        if(memory_usage>0):
            print "test_nsi_iss_oss with Stokes problem FAILS!"
            print "ERROR: Current memory at the end of the execution > 0"
            istat = 1
    elif(key==" Velocity error L2-norm"):
        velocity_norm = float(datum[1].split("\n")[0])
        if(velocity_norm > 1.0e-10):
            print "test_nsi_iss_oss with Stokes problem FAILS!"
            print "ERROR: Velocity error L2-norm > 1.0e-10"
            istat = 1
    elif(key==" Pressure error L2-norm"):
        pressure_norm = float(datum[1].split("\n")[0])
        if(pressure_norm > 1.0e-10):
            print "test_nsi_iss_oss with Stokes problem FAILS!"
            print "ERROR: Pressure error L2-norm > 1.0e-10"
            istat = 1

# Open file
fileid = 'output_test_nsi_iss_oss_navierstokes.txt'
f = open(fileid,'r')

# Read file
data = f.readlines()

# Close file
f.close()

# Extract data
velocity_norm = []
pressure_norm = []
for datum in data:
    datum = datum.split(":")
    key = datum[0]
    if(key==" Current memory usage"):
        memory_usage = float(datum[1].split("\n")[0])
        if(memory_usage>0):
            print "test_nsi_iss_oss with Navier-Stokes problem FAILS!"
            print "ERROR: Current memory at the end of the execution > 0"
            istat = 1
    elif(key==" Velocity error L2-norm"):
        velocity_norm.append(float(datum[1].split("\n")[0]))
    elif(key==" Pressure error L2-norm"):
        pressure_norm.append(float(datum[1].split("\n")[0]))

velocity_rate = []
velocity_rate.append(math.log10(velocity_norm[1]/velocity_norm[0])/math.log10(2))
velocity_rate.append(math.log10(velocity_norm[2]/velocity_norm[1])/math.log10(2))
pressure_rate = []
pressure_rate.append(math.log10(pressure_norm[1]/pressure_norm[0])/math.log10(2))
pressure_rate.append(math.log10(pressure_norm[2]/pressure_norm[1])/math.log10(2))

for i in range(len(velocity_rate)):
    if(abs(3.0-abs(velocity_rate[i])) > 0.05):
        print "test_nsi_iss_oss with Stokes problem FAILS!"
        print "ERROR: Velocity error L2-norm convergence rate < 3"
    if(abs(2.0-abs(pressure_rate[i])) > 0.05):
        print "test_nsi_iss_oss with Stokes problem FAILS!"
        print "ERROR: Pressure error L2-norm convergence rate < 2"

print istat
