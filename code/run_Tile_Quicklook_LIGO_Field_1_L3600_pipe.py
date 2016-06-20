# Original script by STM, revised by AYQH for VLA-GW151226 data

# Takes a rectangular tile and breaks it into rectangular subtiles
# Tiles have boundaries along lines of RA & Dec
# Tile inputs: 
### tile_center_ra,tile_center_dec,
### Num_subtile_ra,Num_subtile_dec,
### L_subtile_arcsec
# Generates locations of subtiles subtile_center_dir



import time
import copy

startRunTime=time.time()
prevRunTime = startRunTime

#================================================================================
# Setup information for the tile and subtiles to be imaged
#================================================================================

# Choose calibrated ms dataset file
# 11 Feb dataset:
# mydataset = '16A-237.sb31782759.eb31845879.57429.90817564815'
# 14 Feb dataset:
mydataset = '16A-237.sb31782757.eb31851884.57432.83502763889'
data_dir = '/lustre/aoc/observers/aho/' 
calibrated_ms = data_dir + mydataset + '.ms' 

# 'auto' option: detect whether there is CORRECTED_DATA and split only this
calibrated_ms_datacolumn = 'auto'

# Use matching with regex in getfieldcone.py 
# Picks out only OTFM fields which start with 0,1,w
use_target_fields = ['^0','^1','^2'] 

# Choose imaging script(s)
use_script_dir = '/users/aho/VLA_GW_Followup/code/'
scriptname = 'run_QuickLook_submosaic_tclean_mfs2048MHz_Pipe.py'
scriptfile = use_script_dir + scriptname

# If you have only generated split/flag/calibration script, 
# set pointer to that script here
# Note: you can generate this using make_restorefile.py
# use_restorescript = data_dir + 'restore_'+ mydataset + '.py'

# Either set these or use preset values 
# Center of Stripe82 tile A at (315.0deg,0deg)
tile_center_epo = 'J2000'
# generated using my script make_tiles.py
# (RA, dec) = (54.70525568 deg, 37. deg)
tile_center_ra = '04:08:00.00' # in hh:mm:ss I assume?
tile_center_dec = '43.00.00.0' # in dd:mm:ss I assume?

#tile_center_ra = '21:00:00.00'
#tile_center_dec = '00.00.00.0'

# Set up the tiling and subtiles
# Choose Num_subtile_ra, Num_subtile_dec, subtile_delta_arcsec to cover the tile area
# L_subtile_arcsec should be equal to or a bit bigger than subtile_delta_arcsec
#
# For 1deg x 1deg subtiles:
# L_subtile_arcsec = 3600.0 # this is size of final subimage
# subtile_delta_arcsec = 3600.0 # this is distance between subtile centers
# Just for testing?
L_subtile_arcsec = 60.0
subtile_delta_arcsec = 60.0

# For all 30deg x 1.23deg of tile
# Num_subtile_ra = 29 
# Num_subtile_dec = 2
# For single central subtile
Num_subtile_ra = 1 
Num_subtile_dec = 1
# Other subtiles:
# L_subtile_arcsec = 2320.0 
# subtile_delta_arcsec = 2320.0 
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

# Choose the pixel size in arcsec. Set this to be at maximum 1/2.5 of the PSF width (FWHM), 
# ideally 1/4 or 1/5 even. For QuickLook PSF~2.5" so 1" here is pushing it.
subtile_pixelsize = 1.0 # arcsec

# Choose "padding" to add to L_subtile_arcsec to make the full image.
# This should be 2 x 0.8 x PB(FWHM). For S-band at 2GHz (low-end of band)
# PB(FWHM)~21' so 2 x 0.8 x 21' ~ 2000" so thats what I use here. Adjust for your band.
# Can specify in pixels or arcsec, prefer arcsec:
# subtile_padding = 2000 # in pixels, added to L_subtile_arcsec for the submosaic imaging
subtile_padding_arcsec = 2000.0 # added to L_subtile_arcsec for the submosaic imaging

subtile_dirname = 'Imaging_'+mydataset
subtile_prefix = 'QuickLook_L'+str(int(L_subtile_arcsec))

# Selection string list for spws for imaging
# This in the dataset after split but before imaging (spw renumbered to 0~15)
# Stripe82 A upper 3 windows bad but not flagged, reverse baseband order
# use_spws = ['0~4,8~15']
use_spws = ['1~14']
# We're commenting this out for the GW follow up, for now.
# Looks like although spw 0 is bad, it's getting a lower weight.

# Use the POINTING table in SDM->MS?
use_pointing = False
# There's no useful pointing table yet

# Enable autoboxing (or not)
# use_ccbox = False
use_ccbox = True
use_maxboxcycles = 3

# Imaging parameters (depend on your observation)
use_threshold = '0.000180Jy' # VLASS depth x1.5, same as GW follow up
# use_restore = '2.5arcsec'    # circular restoring beam size if desired
use_restore = '8.0arcsec' # for the LIGO follow up

# These you tune for quality of imaging
use_fld_cycleniter = 10000
# use_fld_niter = 5000         # max number of clean iter
#use_fld_cycleniter = 750
use_fld_cycleniter = 1000    # max number of iter per major cycle

#------------------------------------------------------------------------------------
# Done setting stuff manually, now run the thing.
#------------------------------------------------------------------------------------
# Derived tile and subtile info
tile_center_dir = me.direction(tile_center_epo,tile_center_ra,tile_center_dec)
tile_center_ra_rad = tile_center_dir['m0']['value']
tile_center_dec_rad = tile_center_dir['m1']['value']
subtile_delta_rad = subtile_delta_arcsec/206264.806
subtile_delta_ra_rad = subtile_delta_rad/pl.cos(tile_center_dec_rad)
#
# For this example pick subtiles manually:
# for j_subtile in [0]:
#     for i_subtile in [5]:
#
# Loop over subtiles, here we do all of them (serially)
for j_subtile in range(Num_subtile_dec):
    for i_subtile in range(Num_subtile_ra):
        # This subtile (i,j) indexed to i=0,Num_ra-1 j=0,Num_dec-1
        # The stuff that figures out where subtiles are:
        d_subtile_ra_rad = (i_subtile - 0.5*(Num_subtile_ra-1.0))*subtile_delta_ra_rad
        d_subtile_dec_rad = (j_subtile - 0.5*(Num_subtile_dec-1.0))*subtile_delta_rad
        #
        ra_subtile_rad = tile_center_ra_rad + d_subtile_ra_rad
        if ra_subtile_rad>(pl.pi):
            ra_subtile_rad = 2.*pl.pi - ra_subtile_rad
        elif ra_subtile_rad<=(-pl.pi):
            ra_subtile_rad = -2.*pl.pi - ra_subtile_rad
        #
        dec_subtile_rad = tile_center_dec_rad + d_subtile_dec_rad
        if dec_subtile_rad>=(0.5*pl.pi):
            dec_subtile_rad = pl.pi - dec_subtile_rad
        #
        subtile_center_dir = copy.deepcopy(tile_center_dir)
        subtile_center_dir['m0']['value'] = ra_subtile_rad
        subtile_center_dir['m1']['value'] = dec_subtile_rad
        #
        # Run the script to image this subtile at location subtile_center_dir:
        # Could also pass as string e.g. by setting
        # subtile_center = 'J2000 00:42:44.30 41.16.09.0'
        execfile(scriptfile)
#
# Done
#
print('Subtile imaging complete')
#
# 
endRunTime=time.time()
TotalRunTime = (endRunTime - startRunTime)
print('Total run time was: %10.3f seconds ' % (TotalRunTime))
