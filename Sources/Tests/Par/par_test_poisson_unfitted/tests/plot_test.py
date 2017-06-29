#!/usr/bin/python

import sys
##import math as m
import numpy as np
import matplotlib.pyplot as plt

#======================================================
def parse_result_file(filename):

   l2_norm                = 0
   h1_semi_norm           = 0
   error_l2_norm          = 0
   error_h1_semi_norm     = 0

   with open(filename) as fp:
       for line in fp:
           line = line.strip()
           line = line.split(':')
           if line[0] == "l2_norm":
               l2_norm = float(line[1])
           if line[0] == "h1_semi_norm":
               h1_semi_norm = float(line[1])
           if line[0] == "error_l2_norm":
               error_l2_norm = float(line[1])
           if line[0] == "error_h1_semi_norm":
               error_h1_semi_norm = float(line[1])

   return l2_norm, h1_semi_norm, error_l2_norm, error_h1_semi_norm

#======================================================
def parse_serie(num_files):

    el2n  = np.zeros(num_files)
    eh1sn = np.zeros(num_files)
    h     = np.zeros(num_files)
    ul2n  = 0.0
    uh1sn = 0.0

    nx0 = 10
    h0  = 2./nx0

    for i in range(num_files):
        filename = "out_mesh" +  str(i+1) + ".txt"
        print "Parsing file " + filename
        l2_norm, h1_semi_norm, error_l2_norm, error_h1_semi_norm = parse_result_file(filename)
        el2n [i] = error_l2_norm
        eh1sn[i] = error_h1_semi_norm
        h    [i] =  h0/(2**(i-1))

    ul2n  = l2_norm
    uh1sn = h1_semi_norm
    return el2n, eh1sn, h, ul2n, uh1sn

#======================================================
def data_print_slope(h,el2n,v,slope):

    dy = v[3]-v[2]
    x0 = np.log10(h[-1])
    x1 = np.log10(h[ 0])
    y0 = np.log10(el2n[-1]) - 0.05*dy
    y1 = y0 + slope*(x1-x0)
    px = np.array([x0,x1])
    py = np.array([y0,y1])

    return px, py

#======================================================
def main(num_files):

    mm2in = 1.0/25.4

    el2n, eh1sn, h, ul2n, uh1sn = parse_serie(num_files)

    el2n  = el2n /ul2n
    eh1sn = eh1sn/uh1sn


    fwidth  = 150*mm2in
    fheight = 200*mm2in
    plt.figure(figsize=(fwidth,fheight))
    plt.plot(np.log10(h),np.log10(el2n),'bo-', label='L2 norm')
    plt.plot(np.log10(h),np.log10(eh1sn),'rs-', label='H1 semi norm')
    plt.xlabel('log10(Element size)')
    plt.ylabel('log10(Relative error)')
    v = plt.axis()
    px, py = data_print_slope(h,el2n,v,2)
    plt.plot(px,py,'k--', label='Slope = 2')
    px, py = data_print_slope(h,eh1sn,v,1)
    plt.plot(px,py,'k-.', label='Slope = 1')
    plt.legend(loc='best')
    plt.savefig('fig_conv_study.pdf',format='pdf')
    plt.close()

    return 0

#======================================================
if __name__ == "__main__":

    # Default values
    num_files = 4

    # Parse arguments
    argv = sys.argv[1:]
    for i in range(len(argv)):
        if argv[i].strip() == "--num_meshes":
            num_files = int(argv[i+1])

    # Run
    exit(main(num_files))
#======================================================
