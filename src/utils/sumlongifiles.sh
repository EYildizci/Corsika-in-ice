#!/bin/bash
#
# job_submit -p1 -cp -t30 -m1000 sumlongifiles.sh001567
#
# sum up all `.long` files:
# ------------------------------------------------------------------------
#                                       juergen.oehlschlaeger@kit.edu
# ------------------------------------------------------------------------
# write names of long files to text file:
  ls -1 DAT001567-*.long > sumlongifiles.i001567
# sum up contents of files:
  ./sumlongifiles < sumlongifiles.i001567 > sumlongifiles.out001567
