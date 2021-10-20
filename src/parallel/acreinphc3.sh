#!/bin/bash
#
# acreinphc3.sh:
# ==============
#     Automatic creation of a single or successive steering files named
#     `parallel-00iiii` and shell script files named `jobhc3-00iiii`
#     to run parallel corsika simulations on the KIT-CS hc3 processors 
#     (i.e. hc3.scc.kit.edu) with MPI parallelization system;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# compilation:
#     f77 -fbounds-check -m32 acreinphc3.f -o acreinphc3
#     gfortran -fbounds-check acreinphc3.f -o acreinphc3
#     ifort -C -check bounds acreinphc3.f -o acreinphc3
# execution:
      ./acreinphc3
      chmod +x jobhc3-*
      # now jobhc3-00iiii, parallel-00iiii exist. 
# submit on hc3:
#     execution of script jobhc3-00iiii:
# ./jobhc3-000100    
#     or direct job_submit (long command): 
# job_submit -cp -p48 -t90 -m2000 -oJob000100_%j.out -eJob000100_%j.err \
#        mpirun mpi_corsika73525_stnd_SIBYLL_urqmd_runner parallel-000100
#
