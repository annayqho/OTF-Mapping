""" Another version of the script from Steve Myers,
this time modified for CASA Pipeline MS """

# set up for rectangular sky survey blocks

# STM 2016-05-31 CASA 4.6.0 run imaging for all subtiles
# Versions:
# STM 2016-06-01 full version, updated scheme with timing
# STM 2016-06-03 use_spws and use_pointing added
# STM 2016-06-07 automatic column, use_target_intents & fields
# STM 2016-06-13 version for MS produced by CASA Integrated Pipeline (CIPL)

import time

startRunTime=time.time()
prevRunTime = startRunTime

#================================================================================
# Setup information for the tile and subtiles to be imaged
#================================================================================
# Dataset inputs: mydataset, data_dir
# Tile inputs: tile_center_ra,tile_center_dec,Num_subtile_ra,Num_subtile_dec,L_subtile_arcsec
# Subtile inputs: i_subtile, j_subtile

# Name of the dataset (for MS naming)
mydataset = 'TSKY0001.sb32154065.eb32157201.57530.47058900463'

# 
# data_dir = '/lustre/aoc/projects/vlass/smyers/Run_TSKY0001.sb32154065.eb32157201.57530/Calibrate_Pipeline4.6.0/'
# This default might work for you if I happened to name it this way
# data_dir = '/lustre/aoc/projects/vlass/smyers/Run_'+mydataset+'/CalScript/'
data_dir = '/lustre/aoc/projects/vlass/smyers/Run_'+mydataset+'/Calibrate_Pipeline_4.6.0/'

# Name of calibrated MS and which datacolumn the calibrated data is in
# calibrated_ms = data_dir + mydataset + '_calibrated_target.ms' # e.g. from script
calibrated_ms = data_dir + mydataset + '.ms' # e.g. from CIPL
# this is the name that comes out of the pipeline

# calibrated_ms_datacolumn = 'all'
# Now there is 'auto' option to detect whether there is CORRECTED_DATA and split only this
# this needs 2016-06-06 version of imaging script to support this
calibrated_ms_datacolumn = 'auto'

# Check intents and field names, only in script v20160607 or later
use_target_intent = '*TARGET*' # picks out scans with this intent only

# Use matching with regex in getfieldcone.py v20160607 or later if this is set
use_target_fields = ['^0','^1','^2'] # picks out only OTFM fields which start with 0,1,2
# this is set up for OTF: by convention, all the target data fields
# start with those

# Location and Name of imaging script(s) to use
# use_script_dir = './'
# here we set to SMyers testing script area
use_script_dir = '/lustre/aoc/projects/vlass/smyers/Scripts/'
# would point this to mine
scriptname = 'run_QuickLook_SubtilePreset_submosaic_tclean_mfs2048MHz_Pipe.py'
scriptfile = use_script_dir + scriptname

# Either set these or use preset values 
# these are for the M31
tile_center_epo = 'J2000'
tile_center_ra = '00:42:44.30'
tile_center_dec = '41.16.09.0'
#
# Choose Num_subtile_ra, Num_subtile_dec, subtile_delta_arcsec to cover the tile area
# L_subtile_arcsec should be equal to or a bit bigger than subtile_delta_arcsec
#
L_subtile_arcsec = 3600.0 # this is size of subtile to image
subtile_delta_arcsec = 3600.0 # this is distance between subtile centers
# in principle you might make these overlapping

# for M31 was 8 deg by 2 deg
Num_subtile_ra = 8
Num_subtile_dec = 2

# L_subtile_arcsec = 2320.0 # this is size of subtile to image
# subtile_delta_arcsec = 2320.0 # this is distance between subtile centers
# Num_subtile_ra = 13
# Num_subtile_dec = 3
# L_subtile_arcsec = 2000.0
# subtile_delta_arcsec = 1800.0
# Num_subtile_ra = 17
# Num_subtile_dec = 3
# L_subtile_arcsec = 1200.0
# subtile_delta_arcsec = 1200.0
# Num_subtile_ra = 25
# Num_subtile_dec = 5

# for a 10" beam, you'd probably want to use 2-3"
subtile_pixelsize = 1.0 # arcsec
# should be 2000"
# tells you how much bigger the actual image is than the subtile size
subtile_padding = 2000 # in pixels, added to L_subtile_arcsec for the submosaic imaging

subtile_dirname = 'Imaging_TSKY0001_sb32154065_57530'
subtile_prefix = 'QuickLook_L3600'

# Selection string list for spws for imaging
# assumes all are good or will be used
use_spws = ['']

# Use the POINTING table in SDM->MS?
# we aren't yet making a pointing table, so leave False
use_pointing = False

# Enable autoboxing (or not)
# Claire Chandler's autoboxing...could play with this, but
# makes it take much longer
# but this image does clearly need clean boxes
# tclean now has its own autoboxing, but haven't played with it
use_ccbox = False

#
# For this example pick a single subtile
# for j_subtile in [0]:
#     for i_subtile in [5]:
#
# Loop over subtiles, here we do all of them (serially)
for j_subtile in range(Num_subtile_dec):
    for i_subtile in range(Num_subtile_ra):
        # This subtile (i,j) indexed to i=0,Num_ra-1 j=0,Num_dec-1
        execfile(scriptfile)
#
# Done
#
print 'Subtile imaging complete'
#
# 
endRunTime=time.time()
TotalRunTime = (endRunTime - startRunTime)
print 'Total run time was: %10.3f seconds ' % (TotalRunTime)
