# Original script by STM, revised by AYQH for VLA-GW151226 data
# 20 June 2016: Condensing

# Takes a rectangular tile and breaks it into rectangular subtiles,
# with boundaries along lines of RA & Dec.
# Then runs the imaging script on these subtiles.

import time
import copy

startRunTime=time.time()
prevRunTime = startRunTime

#================================================================================
# Functions
#================================================================================

def subtile_padding(band):
    """ Padding to add to the subtile size for the full image
    This should be 2 x 0.8 x PB (FWHM)
    
    Parameters
    ----------
    band: string, the frequency band, e.g. S or L

    Returns
    -------
    padding: float, the padding in arcsec
    """

    if band == "S":
        PB = 21 # arcminutes
    else:
        print("you sure you don't want S-band?")
    padding = 2 * 0.8 * 21 * 60 # convert to arcsec
    return padding


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

# Set tile center
tile_center_epo = 'J2000'
# generated using my script make_tiles.py
# (RA, dec) = (54.70525568 deg, 37. deg)
tile_center_ra = '04:08:00.00' # in hh:mm:ss I assume?
tile_center_dec = '43.00.00.0' # in dd:mm:ss I assume?

# Set up the subtiles: number and separation
Num_subtile_ra = 1 # single subtile
Num_subtile_dec = 1
L_subtile_arcsec = 60.0 # in arcseconds (so 3600.0 for a degree)
# size of final subimage
subtile_delta_arcsec = 60.0 # distance btwn subtile centers
if L_subtile_arcsec < subtile_delta_arcsec:
    print("Error: L_subtile_arcsec should be equal to \
            or a bit bigger than subtile_delta_arcsec")

# Choose the pixel size in arcsec. 
# Set this to be at maximum 1/2.5 of the PSF width (FWHM), 
# ideally 1/4 or 1/5 even. 
# For QuickLook PSF~2.5" so 1" here is pushing it.
subtile_pixelsize = 1.0 # arcsec

subtile_padding_arcsec = subtile_padding("S") 

imaging_dirname = 'Imaging_'+mydataset
subtile_prefix = 'QuickLook_L'+str(int(L_subtile_arcsec))

# Selection string list for spws for imaging
# This in the dataset after split but before imaging (spw renumbered to 0~15)
# By eye, looked to be RFI in spw 0 and 15
imaging_spwstr = ['1~14']

# We don't want any specific intents
myintentstr = ''

# For selecting fields in the subtile
dobox = True

# There's no useful pointing table in SDM yet 
clear_pointing = True

# Enable autoboxing (or not)
doccbox = True

# Autoboxing parameters if doccbox = True
if doccbox:
    maxboxcycles = 3
    box_niter = 10000
    box_cyclefactor = 4.5
    box_cycleniter = 500
    peaksnrlimit = 5.0

# UV taper
douvtaper = False
myuvtaper = ['7.0arcsec']

# Number of Taylor terms for MFS
fld_nterms = 1

# Multiscale
fld_multiscale=[0]

# Imaging parameters (specific to this observation)
fld_threshold = '0.000180Jy' # survey depth x1.5
fld_threshold_nobox = fld_threshold
use_restore = '8.0arcsec' # circular restoring beam size

# These you tune for quality of imaging
# clean iterations w/o box
fld_niter = 5000 # max number of clean iter
fld_cyclefactor = 4.5 # max # minor cycle iterations
fld_cycleniter = 1000 # max number of iter per major cycle

# Parameters for processing
doimaging = True
docleanup = False
dousescratch = True
dostats = True

#------------------------------------------------------------------------------------
# Derive tile and subtile locations, and run imaging script
#------------------------------------------------------------------------------------

# Derive tile and subtile info

# Convert RA, Dec, subtile spacing to radians
# (Requires CASA)
tile_center_dir = me.direction(tile_center_epo,tile_center_ra,tile_center_dec)
tile_center_ra_rad = tile_center_dir['m0']['value']
tile_center_dec_rad = tile_center_dir['m1']['value']
subtile_delta_rad = subtile_delta_arcsec/206264.806
# Making the small angle approximation
subtile_delta_ra_rad = subtile_delta_rad/pl.cos(tile_center_dec_rad)

# Loop over subtiles
for j_subtile in range(Num_subtile_dec):
    for i_subtile in range(Num_subtile_ra):
        # Shift to subtile center:
        d_subtile_ra_rad = (i_subtile - 0.5*(Num_subtile_ra-1.0))*subtile_delta_ra_rad
        d_subtile_dec_rad = (j_subtile - 0.5*(Num_subtile_dec-1.0))*subtile_delta_rad
        
        ra_subtile_rad = tile_center_ra_rad + d_subtile_ra_rad
        if ra_subtile_rad>(pl.pi):
            ra_subtile_rad = 2.*pl.pi - ra_subtile_rad
        elif ra_subtile_rad<=(-pl.pi):
            ra_subtile_rad = -2.*pl.pi - ra_subtile_rad
        
        dec_subtile_rad = tile_center_dec_rad + d_subtile_dec_rad
        if dec_subtile_rad>=(0.5*pl.pi):
            dec_subtile_rad = pl.pi - dec_subtile_rad
        
        subtile_center_dir = copy.deepcopy(tile_center_dir)
        subtile_center_dir['m0']['value'] = ra_subtile_rad
        subtile_center_dir['m1']['value'] = dec_subtile_rad
        
        # Run the script to image this subtile at location subtile_center_dir:
        execfile(scriptfile)

print('Subtile imaging complete')

endRunTime=time.time()
TotalRunTime = (endRunTime - startRunTime)
print('Total run time was: %10.3f seconds ' % (TotalRunTime))
