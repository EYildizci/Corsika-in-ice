#!/bin/bash
#
# job_submit -p1 -cp -t350 -m1000 readparticall.sh000082
#
# create list of particle data files and run `readparticall` program:
# ------------------------------------------------------------------------
# usage: ./readparticall.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# cd csk000082/
# 
  ls -1 DAT000082-* | grep t -v | grep n -v > readparticall.i000082
# 
# names of sub paths csk00????;
# gfortran -fbounds-check readparticall.f -o readparticall
# f77 -fbounds-check readparticall.f -o readparticall
# ifort -C -check bounds readparticall.f -o readparticall
#
  ./readparticall < readparticall.i000082 > readparticall.out000082
