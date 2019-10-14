#!/bin/bash
#for i in $( ls ./finalcats ); do

#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f g  
#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f r 
#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f i 
#python ./galsim_no_pixel_run.py -i ./$1/$2 -o ./$1/$2 -f z

python ./galsim_no_pixel_run.py -i ./example_catalog_for_galsim.cat -o test -f g
python ./galsim_no_pixel_run.py -i ./example_catalog_for_galsim.cat -o test -f r
python ./galsim_no_pixel_run.py -i ./example_catalog_for_galsim.cat -o test -f z

