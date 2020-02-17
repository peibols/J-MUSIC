#!/usr/bin/env bash

# combine multiple surface files
for ii in `ls surface_eps_*_*.dat`
do
    jj=`echo $ii | cut -f 1-3 -d _ `
    cat $ii >> $jj.dat
done
rm surface_eps_*_*.dat 2> /dev/null
