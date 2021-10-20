#!/bin/bash
#
# job_submit -p1 -cp -t660 -m1000 showanalyhc3.sh000070
#
# create file list and run `showanalyhc3` program:
# ----------------------------------------------------------------
# cd csk000070/
#
  ls -1 DAT000070-* | grep t -v | grep n -v > showanalyhc3.i000070
# 
# names of sub paths csk00????;
# gfortran -fbounds-check showanalyhc3.f -o showanalyhc3
# f77 -fbounds-check -m32 showanalyhc3.f -o showanalyhc3
# ifort -C -check bounds showanalyhc3.f -o showanalyhc3
#
  ./showanalyhc3 < showanalyhc3.i000070 > showanalyhc3.out000070
  mv fort.9 showanalyhc3.fort000070
