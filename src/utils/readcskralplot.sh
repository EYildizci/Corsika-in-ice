#!/bin/bash
#
# job_submit -p1 -cp -t360 -m1000 readcskralplot.sh
# 
# create file list and run `readcskralplot` program:
# ------------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk001234/
# 
  ls -1 DAT001234-* | grep t -v | grep n -v > readcskralplot.i001234
# 
# names of sub paths csk00????;
# gfortran -fbounds-check readcskralplot.f -o readcskralplot
# ifort -C -check bounds readcskralplot.f -o readcskralplot
#
  ./readcskralplot < readcskralplot.i001234 > readcskralplot.out001234
  mv fort.19 readcskralplot.fort001234
