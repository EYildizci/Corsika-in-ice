#!/bin/bash
#
# showsimulist.sh:
# ================
#           create tabular of air shower simulations of corsika
#           structure by reading first record of all particle
#           data files DATiiiiii and additional reading of the
#           corresponding protocol file DATiiiiii.lst .
# --------------------------------------------------------------------
#                                   juergen.oehlschlaeger@kit.edu
# --------------------------------------------------------------------
#
# prim, lg(E), theta, phi, nsh, task, size, obslev, h1stme, thilev,
#       thiwmax, lg(thirad), verspgm, models, rundate, Xmagn, Zmagn.
#
  ls -l DAT?????? | ./showsimulist > showsimulist.cd3-e20m100thi
#
#       f77 -fbounds-check -m32 showsimulist.f -o showsimulist
#       gfortran -fbounds-check showsimulist.f -o showsimulist
#       ifort -C -check bounds showsimulist.f -o showsimulist
#
