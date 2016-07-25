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
        padding = 2 * 0.8 * 21 * 60 # convert to arcsec
        # this gives 2016, need to round for efficient imagesize
        return 2000.0
    else:
        print("you sure you don't want S-band?")


#================================================================================
# Setup information for the tile and subtiles to be imaged
#================================================================================

# Choose calibrated ms dataset file
# 11 Feb dataset:
# mydataset = '16A-237.sb31782759.eb31845879.57429.90817564815'
# 14 Feb dataset:
mydataset = '16A-237.sb31782757.eb31851884.57432.83502763889'
# data_dir = '/lustre/aoc/observers/aho/VLASS_Field/' 
data_dir = '/lustre/aoc/observers/aho/Run_16A-237.sb31782757.eb31851884.57432.83502763889'

# sample VLASS field, for testing PyBDSM mask
# mydataset = 'TSKY0001.sb32295801.eb32296475.57549.31722762731'

calibrated_ms = data_dir + mydataset + '.ms' 

# 'auto' option: detect whether there is CORRECTED_DATA and split only this
calibrated_ms_datacolumn = 'auto'

# for VLASS only:
# use_target_intent = '*TARGET*' # picks out scans with this intent only

# Use matching with regex in getfieldcone.py 
# Picks out only OTFM fields which start with 0,1,w
use_target_fields = ['^0','^1','^2'] 

# Choose imaging script(s)
use_script_dir = '/users/aho/VLA_GW_Followup/code/'
#scriptname = 'image_tile.py'
scriptname = 'run_QuickLook_submosaic_tclean_mfs2048MHz_Pipe.py'
scriptfile = use_script_dir + scriptname

# Set tile center
tile_center_epo = 'J2000'
# generated using my script make_tiles.py
# (RA, dec) from Kunal = (54.70525568 deg, 37. deg)
# corresponds to:
# tile_center_ra = '03:38:49.26'
# tile_center_Dec = '37:00:00.00'
tile_center_ra = '04:08:00.00' # in hh:mm:ss I assume?
tile_center_dec = '43.00.00.0' # in dd:mm:ss I assume?

# from Steve's VLASS script:
# tile_center_ra = '21:00:00.00'
# tile_center_dec = '00.00.00.0'

# Set up the subtiles: number and separation
L_subtile_arcsec = 3600.0 # in arcseconds (so 3600.0 for a degree)
# size of final subimage
subtile_delta_arcsec = 3600.0 # distance btwn subtile centers
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
#subtile_prefix = 'QuickLook_L'+str(int(L_subtile_arcsec))

# for Steve's VLASS box:
subtile_prefix = 'QuickLook_CCBox_L'+str(int(L_subtile_arcsec))

# Selection string list for spws for imaging
# This in the dataset after split but before imaging (spw renumbered to 0~15)
# By eye, looked to be RFI in spw 0 and 15
imaging_spwstr = ['1~14']

# for Steve's VLASS box:
# imaging_spwstr = ['0~4,8~15']
#imaging_spwstr = ['11~12']

# We don't want any specific intents
myintentstr = ''

# For selecting fields in the subtile
dobox = True

# There's no useful pointing table in SDM yet 
clear_pointing = True

# Enable autoboxing (or not)
# doccbox = False
doccbox = True

# Enable input mask (or not)
# mask = '/lustre/aoc/observers/aho/VLASS_Field/img.TSKY0001.sb32295801.eb32296475.57549.31722762731_QuickLook_CCBox_L3600.clean.mask'
# mask = 'test.pybdsm_gaus_model.image'
mask = ''

# Autoboxing parameters if doccbox = True
if doccbox:
    maxboxcycles = 3
    box_niter = 10000
    box_cyclefactor = 4.5
    box_cycleniter = 500
    peaksnrlimit = 5.0

# Masking parameters if mask != ''
if mask != '':
    maxboxcycles = 3
    box_niter = 10000
    box_cyclefactor = 4.5
    box_cycleniter = 500
    peaksnrlimit = 5.0


# Widefield parameters
fld_wprojplanes=1
fld_facets=1
fld_gridder='mosaic'
fld_pblimit=0.2
fld_normtype='flatnoise'

# More imaging stuff
myweight = 'briggs'
myrmode='norm'
myrobust=1.0

# Save model
dosavemodel = 'modelcolumn'

# Deconvolver and multiscale
# fld_deconvolver='hogbom'
# for Multi-Taylor:
use_multiscale = True
fld_deconvolver='mtmfs'
fld_specmode = 'mfs'
fld_reffreq = '3.0GHz'

# Number of Taylor terms for MFS
fld_nterms = 2
use_nterms = 2

# UV taper
douvtaper = False
if douvtaper == True:
    myuvtaper = ['7.0arcsec']
elif douvtaper == False:
    myuvtaper = []
else:
    print("uvtaper set incorrectly?")

# Multiscale
# fld_multiscale = [0]
fld_multiscale=[0,3,10]
use_multiscale = [0,3,10]
# what I did for Kunal...multiscale [0,3,10]

# Imaging parameters (specific to this observation)
fld_threshold = '0.000180Jy' # survey depth x1.5, GW same as VLASS
fld_threshold_nobox = fld_threshold

# for LIGO
use_restore = '8.0arcsec' # circular restoring beam size
# for Steve's VLASS box:
# use_restore = '2.5arcsec'    # circular restoring beam size if desired

# These you tune for quality of imaging
# clean iterations w/o box
fld_niter = 10000 # max number of clean iter
fld_cyclefactor = 4.5 # max # minor cycle iterations
fld_cycleniter = 1000 # max number of iter per major cycle

# Parameters for processing
doimaging = True
docleanup = False
dousescratch = True
dostats = True
# parallel =  True
parallel = False

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

if tile_center_ra_rad>(pl.pi):
    tile_center_ra_rad = 2.*pl.pi - tile_center_ra_rad
elif tile_center_ra_rad <= (-pl.pi):
    tile_center_ra_rad = -2.*pl.pi - tile_center_ra_rad 

if tile_center_dec_rad >= (0.5*pl.pi):
    tile_center_dec_rad = pl.pi - tile_center_dec_rad 

subtile_center_dir = copy.deepcopy(tile_center_dir)
subtile_center_dir['m0']['value'] = tile_center_ra_rad
subtile_center_dir['m1']['value'] = tile_center_dec_rad 

# Run the script to image this subtile at location subtile_center_dir:
execfile(scriptfile)

print('Subtile imaging complete')

endRunTime=time.time()
TotalRunTime = (endRunTime - startRunTime)
print('Total run time was: %10.3f seconds ' % (TotalRunTime))
