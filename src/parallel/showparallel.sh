#!/bin/bash
# 
# create a short tabular of parallel corsika simulations:
# ------------------------------------------------------------------------
# Primary   lg(E)  theta    phi  runtsk  sizeGBy  procs
#    dectcut  dectmax  t(min)  files  RATIO  obslev  Xmagn  Zmagn
#       ecutha  ecutmu  ecutel  ecutga  thilev  wmax  lg(thirad)
# ------------------------------------------------------------------------
# usage: ./showparallel.sh
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# 
  ls -1 csk00*/Job00*.out > showparallel.jobinfos
# 
# names of sub paths csk00????;
# hc3 job protocols Job00????_%jobid.err, Job00????_%jobid.out;
#
# gfortran -fbounds-check showparallel.f -o showparallel
# f77 -fbounds-check -m32 showparallel.f -o showparallel
# ifort -C -check bounds showparallel.f -o showparallel
#
# ----- in case of csk002141/Job00*.err does not contain current number of files:
# ls -l csk002141/DAT00*.lst | wc | awk '{printf("%7d\n",$1)}' > csk002141/Job00*.err
#
  ./showparallel < showparallel.jobinfos > showparallel.l16-work-jobinfos
#
