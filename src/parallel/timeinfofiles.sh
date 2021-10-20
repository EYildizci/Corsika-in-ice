#!/bin/bash
# 
# timeinfofiles.sh:
# =================
#     rename all still identical names `time.txt` of all csk00????
#     subdirectories to corresponding names like `time.txt00????`,
#     only if not yet done; version for hc3 (i.e. hc3.scc.kit.edu),
#     but also valid on /data/corsdat*/joe/ paths;
#     parallel run numbers taken up to 009999.
# -----------------------------------------------------------------
#                                juergen.oehlschlaeger@kit.edu
# -----------------------------------------------------------------
# - - - - write original `csk00????/time.txt` name as 7th line: 
  for timename in $( ls -1 csk00????/time.txt ); do 
    nlines=`wc $timename | awk '{ printf("%7d",$1) }'`
    if [ $nlines -eq 6 ]; then
      ls -1 $timename >> $timename 
    fi
  done
# - - - - write list of original csk00????/time.txt names:
  filist=`echo "timeinfofiles.inptxt" | awk '{ printf("%s",$1) }'`
  /bin/ls -1 csk00????/time.txt > $filist
  nfiles=`wc $filist | awk '{ printf("%7d",$1) }'`
  if [ $nfiles -gt 0 ]; then
# - - - - rename all files not yet renamed by run number: 
    for timename in $( cat $filist ); do   
      echo "corsika" $timename > timefile-nnnnnn.tmp
      jtsize=`wc timefile-nnnnnn.tmp |  awk '{ printf("%d",$3) }'`
      if [ $jtsize -lt 28 ]; then
        sed -i "/corsika/ s/csk/   /" timefile-nnnnnn.tmp 
        sed -i "/corsika/ s/\// /" timefile-nnnnnn.tmp 
        runnr=`grep "corsika" timefile-nnnnnn.tmp | awk '{ printf("%d",$2) }'`
        /bin/rm timefile-nnnnnn.tmp
        newname=`echo $timename $runnr | awk '{ printf("%s%06d",$1,$2) }'`
        /bin/mv $timename $newname
        echo $newname
      fi
    done
  else
# - - - - no files have to be renamed:
    echo "          ..." $filist "is empty ... "
  fi
#
