#!/bin/bash

cd /lustre/kmooley/projects/JAGWAR/GW151226/anna

export PATH=$PATH:/lustre/kmooley/programs/aipslite/parseltongue-2.0-dev/bin/:/lustre/kmooley/programs/aipslite/parseltongue-2.0-dev/bin/:/lustre/kmooley/programs/aipslite/obit-1.1.486-1-64b

export PYTHONPATH=$PYTHONPATH:/lustre/kmooley/programs/aipslite/core:/lustre/kmooley/programs/python/pyfits-3.1/lib/python2.7/site-packages:/lustre/kmooley/programs/python/pywcs-1.11-4.8.2/lib/python2.7/site-packages:/lustre/kmooley/programs/python/pyds9-1.7/installation/lib/python2.6/site-packages:/lustre/kmooley/pipeline/CASA/scripts:/lustre/kmooley/programs/python/pyephem/ephem-3.7.5.3/instl/lib64/python2.6/site-packages

# --wprojplanes 32

#for NITER in `seq 1000 1000 5000`; do
for NITER in `seq 5000 5000 5000`; do
    COMMAND="/lustre/kmooley/pipeline/CASA/scripts/cleanfld.py \
        --niter $NITER --imsize 1500 --pixsize 2 --clipmax 0 \
        --cyclefactor 4.5 --beam 8.0 \
        --wprojplanes 32
        --multiscale [0,3,10]
        --parallel /lustre/kmooley/projects/JAGWAR/GW151226/anna"
    echo $COMMAND
    $COMMAND
    COMMAND="/lustre/kmooley/pipeline/AIPS/flatn.py \
        --imsize 2000,2000 --coord 62.00,43.00"
    echo $COMMAND
    $COMMAND
    rm -f 2*.fits
    rm -f 3*.fits
    DIR="Run_wpp_Multiscale_NITER_$NITER"
    mkdir $DIR
    mv *.fits $DIR
    mv *.log $DIR
    rm -rf data
done 


