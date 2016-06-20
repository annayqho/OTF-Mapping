# Original script from STM, modified by AYQH

import copy
import time
import os
from os import F_OK

#================================================================================
logbuffer = []
mytext = 'VLASS QuickLook submosaic TClean MFS TT0-only JointMosaic'
myscriptvers = '2016-06-14 STM (4.6.0)'

logstring = mytext
print(logstring)
logbuffer.append(logstring)
logstring = myscriptvers
print(logstring)
logbuffer.append(logstring)

#================================================================================
# Set submosaic center
#================================================================================

# subtile_center_dir #
# was created in the run_Tile script
mycenter_dir = copy.deepcopy(subtile_center_dir)
epo = mycenter_dir['refer']
rapos = mycenter_dir['m0']
decpos = mycenter_dir['m1']
ral = qa.angle(rapos,form=["tim"],prec=9)
decl = qa.angle(decpos,prec=10)
mycenter = epo + ' ' + ral[0] + ' ' + decl[0]
logstring = 'SubMosaic center = '+mycenter
print(logstring)
logbuffer.append(logstring)

#================================================================================
# Other setup stuff specific to this dataset
#================================================================================
# calibrated_ms was created in run_Tile
splitfile = calibrated_ms
if os.access(splitfile,F_OK):
    logstring = 'Found calibrated ms '+splitfile
else:
    logstring = 'ERROR: could not find '+splitfile
print(logstring)
logbuffer.append(logstring)

# calibrated_ms_datacolumn was created in run_Tile
splitdatacolumn = calibrated_ms_datacolumn
logstring = 'Will use datacolumn '+splitdatacolumn
print(logstring)
logbuffer.append(logstring)


# NOTE: this splitfile is assumed to contain only the target fields that you
# want to image from, though in principle the getfieldirbox will pull them from
# even the full pipeline calibrated MS (though you may pull in the setup fields
# near those corners. What I do is run:
# msfile = mydataset+'.ms'
# splitfile = mydataset+'_calibrated_target.ms'
# targetfields = '0*' # picks out only OTFM fields which for me happen to start with 0*
# targetintent = '*TARGET*'
# mstransform(vis=msfile,outputvis=splitfile,datacolumn='corrected',field=targetfields,intent=targetintent,correlation='RR,LL') # not doing polarization yet

# Option to clear POINTING table before imaging
if clear_pointing:
    logstring = 'Will clear MS POINTING table(s)'
else:
    logstring = 'MS POINTING table(s) will be used if they contain data'
print(logstring)
casalog.post(logstring)
logbuffer.append(logstring)

#================================================================================
# Setup stuff not dataset specific
#================================================================================

# location for writing out images
# subtile_prefix & imaging_dirname were set in run_Tile 
postfix ='_'+ subtile_prefix
tiling='_subtile_%i_%i' % (i_subtile,j_subtile)
scriptprefix=mydataset+postfix+tiling
imaging_dir = imaging_dirname+postfix+tiling

# use_script_dir is set 
execfile(use_script_dir+'getfieldcone.py')

# Set up some parameters for processing
sdmfile = calibrated_ms.split('.ms')[0]
workfile = mydataset + '_calibrated_target_working.ms'
doimaging = True
docleanup = False
dousescratch = True
dosavemodel = 'modelcolumn'
dostats = True

# imaging parms
visname = imaging_dir+'/'+workfile
clnname = imaging_dir+'/img.'+outname+postfix+'.clean'
dirtyname = imaging_dir+'/img.'+outname+postfix+'.dirty'

dostartmodel = True

# Widefield parameters
fld_wprojplanes=1
fld_facets=1
fld_gridder='mosaic'
fld_pblimit=0.2
fld_normtype='flatnoise'
fld_wbawp = True

# Deconvolver and multiscale
fld_deconvolver='hogbom'

# If multi-Taylor term
# fld_deconvolver='mtmfs'
# Number of Taylor terms for MFS
if 'use_nterms' in locals() or 'use_nterms' in globals():
    fld_nterms = use_nterms
else:
    fld_nterms = 1
# Multiscale
if 'use_multiscale' in locals() or 'use_multiscale' in globals():
    fld_multiscale=copy.deepcopy(use_multiscale)
else:
    fld_multiscale=[0]

fld_specmode = 'mfs'
fld_reffreq = '3.0GHz'

if 'use_threshold' in locals() or 'use_threshold' in globals():
    fld_threshold = use_threshold
    fld_threshold_q = qa.quantity(fld_threshold)
    fld_threshold_qJy = qa.convert(fld_threshold_q,'Jy')
    fld_thresholdJy = fld_threshold_qJy['value']
    #
    fld_threshold_nobox = use_threshold
    fld_threshold_nobox_q = qa.quantity(fld_threshold_nobox)
    fld_threshold_nobox_qJy = qa.convert(fld_threshold_nobox_q,'Jy')
    fld_thresholdJy_nobox = fld_threshold_nobox_qJy['value']
else:
    # Based on MFS sensitivity 120uJy 1500MHz bw
    # 1.5sigma
    # fld_thresholdJy = 0.000450 # for 120MHz bw
    fld_thresholdJy = 0.000180 # for 1500MHz bw
    fld_threshold = str(fld_thresholdJy)+'Jy'
    # Without box
    # fld_thresholdJy_nobox = 0.000450 # for 120MHz bw
    fld_thresholdJy_nobox = 0.000180 # for 1500MHz bw
    fld_threshold_nobox = str(fld_thresholdJy_nobox)+'Jy'

# clean iterations w/o box
if 'use_fld_niter' in locals() or 'use_fld_niter' in globals():
    fld_niter = use_fld_niter
else:
    fld_niter = 5000
if 'use_fld_cyclefactor' in locals() or 'use_fld_cyclefactor' in globals():
    fld_cyclefactor = use_fld_cyclefactor
else:
    fld_cyclefactor = 4.5
# Max number minor cycle iterations
# Set this lower to compensate for changing PSF over submosaic
if 'use_fld_cycleniter' in locals() or 'use_fld_cycleniter' in globals():
    fld_cycleniter = use_fld_cycleniter
else:
    # fld_cycleniter = 750
    fld_cycleniter = 1000

# B-config 3 GHz (beam approx 2.5' robust 1)
# for fieldcone = 0.25deg subim size 1deg = 3600"
# for fieldcone = 0.5deg subim size 1.5deg = 5400"
# for fld_subim_size = 2000" want fieldcone = 2000" and fld_size = 4000
# for fld_subim_size = 2400" want fieldcone = 2200" and fld_size = 4400
# fld_size = 4000
# cellsize = 1.0
cellsize = subtile_pixelsize
fld_cell = str(cellsize)+'arcsec'
# fld_size = 4320
L_subtile_pixels = int(L_subtile_arcsec/subtile_pixelsize)
if 'subtile_padding' in locals() or 'subtile_padding' in globals():
    padding = subtile_padding
    padding_arcsec = cellsize*padding
else:
    if 'subtile_padding_arcsec' in locals() or 'subtile_padding_arcsec' in globals():
        padding_arcsec = subtile_padding_arcsec
    else:
        padding_arcsec = 2000.0
    #
    padding = int(padding_arcsec/cellsize)
#
fld_size = L_subtile_pixels + padding
logstring = 'Using field image size %i with cell size %s ' % (fld_size,fld_cell)
print(logstring)
logbuffer.append(logstring)

dosubim = True
# for fieldcone = 0.25deg subim size 0.25deg = 900"
#fld_subim_size=1800
#fld_subim = '900,900,2699,2699'
# for fieldcone = 0.5deg subim size 0.5deg = 1800"
# for fld_subim_size = 2000" want fieldcone = 2000" and fld_size = 4000
# for fld_subim_size = 2400" want fieldcone = 2200" and fld_size = 4400
# fld_subim_size=2000
# fld_subim_size=2320
fld_subim_size=L_subtile_pixels
ilow = fld_size/2 - (fld_subim_size/2)
iup = fld_size/2 + (fld_subim_size/2) - 1
fld_subim = str(ilow)+','+str(ilow)+','+str(iup)+','+str(iup)
logstring = 'Using field subimage blc,trc of (%i,%i) ' % (ilow,iup)
print(logstring)
logbuffer.append(logstring)

# Use common restoring beam?
if 'use_restore' in locals() or 'use_restore' in globals():
    fld_bmaj = use_restore
    fld_bmin = use_restore
    dorestore = not use_restore==''
    fld_bpa = '0deg'
else:
    # S-band B-config
    dorestore=True
    fld_bmaj = '2.5arcsec'
    fld_bmin = '2.5arcsec'
    fld_bpa = '0deg'

# B-config 3 GHz (tapered to 7.5')
# douvtaper=True
# myuvtaper=['7.0arcsec']
# fld_size = 1200
# fld_cell = '2.5arcsec'
# fld_niter = 400
# fld_threshold = '0.0006Jy'
# fld_cyclefactor = 3.5
# fld_subim = '150,150,1050,1050'
# fld_bmaj = '7.5arcsec'
# fld_bmin = '7.5arcsec'
# fld_bpa = '0deg'

#myweight = 'natural'
#myweight = 'uniform'
myweight = 'briggs'
myrmode='norm'
#myrobust=0.5
myrobust=1.0

print('')
logbuffer.append(' ')

#====================================================================
# list of fields to be cleaned
#mynumfields = 1280
#fldnos = range(mynumfields)
#fldnos = [758]
fldnos = []
#
#docone=False
# Use field (box/cone)
beamsearchradius_arcsec = 1000.0
docone=True
# If box Actual dist is linmos_subim_size + 1000"
# mydistance='2000arcsec'
mydistance_arcsec = 0.5*L_subtile_arcsec + beamsearchradius_arcsec
mydistance=str(mydistance_arcsec)+'arcsec'
print(mydistance)
# 
if docone:
    logstring = 'Selecting fields within cone radius of '+mydistance
    print(logstring)
    logbuffer.append(logstring)

#====================================================================
# Save the parameters used 
params = {}
# Version of script
params['version'] = myscriptvers
# User set params
params['user'] = {}
params['user']['sdmfile'] = sdmfile
params['user']['outname'] = outname
# Tiling parameters
params['user']['subtile_center'] = mycenter
params['user']['L_subtile_arcsec'] = L_subtile_arcsec
params['user']['subtile_delta_arcsec'] = subtile_delta_arcsec
params['user']['padding_arcsec'] = padding_arcsec
#
if doimaging:
    params['user']['imaging_spwstr'] = imaging_spwstr
    params['user']['imaging_dir'] = imaging_dir
    if docone:
        params['user']['mydistance'] = mydistance
    params['user']['fldnos'] = fldnos
    params['user']['fld_size'] = fld_size
    params['user']['fld_cell'] = fld_cell
    params['user']['fld_niter'] = fld_niter
    params['user']['fld_cycleniter'] = fld_cycleniter
    params['user']['fld_threshold_nobox'] = fld_threshold_nobox
    params['user']['fld_cyclefactor'] = fld_cyclefactor
    params['user']['dosubim'] = dosubim
    if dosubim:
        params['user']['fld_subim'] = fld_subim
    params['user']['submos_dir'] = mycenter_dir
    params['user']['fld_deconvolver'] = fld_deconvolver
    if fld_deconvolver=='mtmfs':
        params['user']['fld_nterms'] = fld_nterms
        params['user']['fld_multiscale'] = fld_multiscale
    params['user']['fld_specmode'] = fld_specmode
    if fld_specmode=='mfs':
        params['user']['fld_reffreq'] = fld_reffreq
    else:
        params['user']['fld_nchan'] = fld_nchan
        params['user']['fld_start'] = fld_start
        params['user']['fld_width'] = fld_width
        params['user']['fld_interp'] = fld_interp
        params['user']['fld_reffreq'] = fld_reffreq
    params['user']['fld_gridder'] = fld_gridder
    if fld_gridder=='widefield' or fld_gridder=='awproject':
        params['user']['fld_nterms'] = fld_nterms
        params['user']['fld_wprojplanes'] = fld_wprojplanes
        params['user']['fld_facets'] = fld_facets
    if dorestore:
        params['user']['fld_bmaj'] = fld_bmaj
        params['user']['fld_bmin'] = fld_bmin
        params['user']['fld_bpa'] = fld_bpa
    params['user']['douvtaper'] = douvtaper
    if douvtaper:
        params['user']['myuvtaper'] = myuvtaper
    params['user']['myweight'] = myweight
    if myweight=='briggs':
        params['user']['myrmode'] = myrmode
        params['user']['myrobust'] = myrobust
    # Autoboxing parameters
    params['user']['doccbox'] = doccbox
    if doccbox:
        params['user']['box_niter'] = box_niter
        params['user']['box_cycleniter'] = box_cycleniter
        params['user']['box_cyclefactor'] = box_cyclefactor
        params['user']['peaksnrlimit'] = peaksnrlimit
        params['user']['maxboxcycles'] = maxboxcycles
        params['user']['fld_threshold'] = fld_threshold
      

stagename = []
stagetime = []
steplist = []
steptimes = {}

logbuffer = []
if dostats:
    statbuffer = []

#====================================================================
# Start actual processing
#====================================================================
startTime=time.time()
startProc=time.clock()
prevTime = startTime

logstring = 'Starting imaging of '+splitfile+' using script '+myscriptvers
print(logstring)
casalog.post(logstring)
logbuffer.append(logstring)

# get dataset sizes in MB (1024x1024 bytes!)
f = os.popen('du -ms '+splitfile)
fstr = f.readline()
f.close()
datasize_ms = float( fstr.split("\t")[0] )
logstring = 'MS is '+str(datasize_ms)+' MB'
print(logstring)
casalog.post(logstring)
logbuffer.append(logstring)

#
if os.path.exists(imaging_dir):
    print('Using existing image output directory '+imaging_dir)
else:
    print('Creating directory '+imaging_dir)
    os.makedirs(imaging_dir)

if docone:
   
    # fldnos = getfieldcone(splitfile,distance=mydistance,center=mycenter)
    # Use direction measure to subtile center
    # fldnos = getfieldircone(splitfile,distance=mydistance,center_dir=subtile_center_dir)
    # Use bounding box distance rather than cone
    # Check if there are regex to match with
    if 'use_target_fields' in locals() or 'use_target_fields' in globals():
        mymatchregex = use_target_fields
    else:
        mymatchregex = ''
    if mymatchregex=='' or mymatchregex==[]:
        # no name match, backward compatible
        fldnos = getfieldirbox(splitfile,distance=mydistance,center_dir=mycenter_dir)
    else:
        # use regex match, needs getfieldcone.py v20160607 or later
        logstring = 'Matching field names using regex string(s) : '+str(mymatchregex)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        fldnos = getfieldirbox(splitfile,distance=mydistance,center_dir=mycenter_dir,matchregex=mymatchregex)
logstring = 'Will image a total of '+str(len(fldnos))+' fields'+' using specmode='+fld_specmode
print(logstring)
casalog.post(logstring)
logbuffer.append(logstring)
logstring = 'Will image fields = '+str(fldnos)
print(logstring)
casalog.post(logstring)
logbuffer.append(logstring)
#
# Construct field selection
fldstrs = ''
for field in fldnos:
    if fldstrs=='':
        fldstrs = str(field)
    else:
        fldstrs = fldstrs + ',' + str(field)

if doimaging:
    #
    # Make a working copy of ms
    print('Splitting '+splitfile+' to '+visname)
    os.system('rm -rf '+visname+'*')
    #os.system('cp -rf '+splitfile+' '+visname)
    if splitdatacolumn=='auto':
        # check if CORRECTED_DATA is there use that if so
        tb.open(splitfile)
        split_colnames = tb.colnames()
        tb.close()
        if 'CORRECTED_DATA' in split_colnames:
            mydatacolumn = 'corrected'
        else:
            mydatacolumn = 'data'
    else:
        mydatacolumn = splitdatacolumn
    # 
    # Check for specific intents
    if 'use_target_intent' in locals() or 'use_target_intent' in globals():
        myintentstr = use_target_intent
    else:
        myintentstr = ''
    #
    if 'use_restorescript' in locals() or 'use_restorescript' in globals():
        # If requested, run instead the restore/split script generated beforehand
        stepname = 'split/restore'
        targetfields = fldstrs
        targetintent = myintentstr
        usemsfile = splitfile
        targetfile = visname
        execfile(use_restorescript)
    else:
        stepname = 'split'
        mstransform(splitfile,visname,field=fldstrs,intent=myintentstr,datacolumn=mydatacolumn)
    #
    # Make a listing
    listobs(visname,listfile=visname+'.listobs')
    #
    # If requested clear POINTING table in ms
    if clear_pointing:
        # Remove the POINTING tables
        # for name in visnamelist:
        tb.open(visname+'/POINTING', nomodify = False)
        a = tb.rownumbers()
        tb.removerows(a)
        tb.close()
        #
        logstring = 'Cleared MS POINTING table(s) for '+visname
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
    #
    # Extract info from SPECTRAL_WINDOW subtable
    tb.open(visname+"/SPECTRAL_WINDOW")
    nchanarr=tb.getcol("NUM_CHAN")
    spwnamearr=tb.getcol("NAME")
    reffreqarr=tb.getcol("REF_FREQUENCY")
    bwarr=tb.getcol("TOTAL_BANDWIDTH")
    chanfreqarr=tb.getcol("CHAN_FREQ")
    chanwidtharr=tb.getcol("CHAN_WIDTH")
    tb.close()
    nspw = len(nchanarr)
    print('Found '+str(nspw)+' spw in SPECTRAL_WINDOW table')
    spwlookup = {}
    for isp in range(nspw):
        spwlookup[isp] = {}
        spwlookup[isp]['nchan'] = nchanarr[isp]
        spwlookup[isp]['name'] = str( spwnamearr[isp] )
        spwlookup[isp]['reffreq'] = reffreqarr[isp]
        spwlookup[isp]['bandwidth'] = bwarr[isp]
        nch=nchanarr[isp]
        bw=bwarr[isp]
        rfq=reffreqarr[isp]
        cfq=rfq+(0.5*float(nch-1)/float(nch))*bw
        spwlookup[isp]['centerfreq'] = cfq
        #
        bw_str = '%.3fMHz' % (bw/1.e6)
        cfq_str = '%.3fMHz' % (cfq/1.e6)
        spwlookup[isp]['centerfreq_str'] = cfq_str
        spwlookup[isp]['bandwidth_str'] = bw_str
    print('Extracted information for '+str(nspw)+' SpectralWindows')
    #
    currTime=time.time()
    stagedur = currTime-prevTime
    # stepname = 'split'
    if steptimes.has_key(stepname):
        steptimes[stepname]+=stagedur
    else:
        steplist.append(stepname)
        steptimes[stepname]=stagedur
    stagestr = stepname
    stagetime.append(stagedur)
    stagename.append(stagestr)
    print(stagestr+' took '+str(stagedur)+' sec')
    prevTime = currTime

# Make standard (untapered) image set
startimagingTime = currTime
if doimaging:
    #
    spwstr = imaging_spwstr
    foundbox=False
    #
    # Joint MFS mosaic of all fields in split MS
    fldlist = range(len(fldnos))
    fieldstr = ''
    clnim = clnname
    dirtyim = dirtyname
    #
    # Frequency info
    # spw_center = spwlookup[ispw]['centerfreq_str']
    # spw_width = spwlookup[ispw]['bandwidth_str']
    #
    # clean up previous images
    os.system('rm -rf '+clnim+'.*')
    os.system('rm -rf '+dirtyim+'.*')
    #
    if fld_nterms>1:
        clnrestored=[]
        clnresidual = []
        clnmodel=[]
        clnpsf = []
        clnpb = [clnim+'.pb']
        dirtyresidual = []
        for tt in range(fld_nterms):
            myim=clnim+'.image.tt'+str(tt)
            clnrestored.append(myim)
            resim=clnim+'.residual.tt'+str(tt)
            clnresidual.append(resim)
            modim=clnim+'.model.tt'+str(tt)
            clnmodel.append(modim)
            psfim=clnim+'.psf.tt'+str(tt)
            clnpsf.append(psfim)
            dirtyresim=dirtyim+'.residual.tt'+str(tt)
            dirtyresidual.append(dirtyresim)
    else:
        clnrestored=[clnim+'.image']
        clnresidual = [clnim+'.residual']
        clnmodel=[clnim+'.model']
        clnpsf = [clnim+'.psf']
        clnpb = [clnim+'.pb']
        dirtyresidual = [dirtyim+'.residual']
    #
    if dorestore:
        myrestore=[fld_bmaj,fld_bmin,fld_bpa]
    else:
        myrestore=[]
    iterdone=0
    itercycle=0
    # ==========================================
    # The part that (optionally) does autoboxing
    # ==========================================
    if doccbox:
        #
        # make a dirty image to begin with
        #
        try:
            # this is where we create the primary beam .pb image
            tclean(visname,
                   imagename=clnim,
                   field=fieldstr,
                   spw=spwstr,
                   imsize=[fld_size,fld_size],
                   cell=[fld_cell,fld_cell],
                   phasecenter=mycenter,
                   stokes='I',
                   startmodel='',
                   specmode=fld_specmode,
                   reffreq=fld_reffreq,
                   # nchan=fld_nchan,
                   # start=spw_center,
                   # width=spw_width,
                   # interpolation=fld_interp,
                   # resttfreq=[fld_reffreq],
                   gridder=fld_gridder,
                   pblimit=fld_pblimit,
                   normtype=fld_normtype,
                   # wbawp=fld_wbawp,
                   # wprojplanes=fld_wprojplanes,
                   # facets=fld_facets,
                   deconvolver=fld_deconvolver,
                   # scales=fld_multiscale,
                   restoringbeam=myrestore,
                   niter=0,
                   threshold=fld_threshold,
                   cycleniter=fld_cycleniter,
                   cyclefactor=fld_cyclefactor,
                   usemask='user',
                   mask='',
                   interactive=False,
                   weighting=myweight,
                   robust=myrobust,
                   uvtaper=myuvtaper,
                   makeimages='choose',
                   calcres=True,
                   calcpsf=True,
                   writepb=True,
                   savemodel=dosavemodel)
        except:
            logstring = 'WARNING: Failed creating dirty submosaic'
            print(logstring)
            casalog.post(logstring)
            logbuffer.append(logstring)
        #
        # Save parameters from this run
        os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
        #
        currTime=time.time()
        stagedur = currTime-prevTime
        stepname = 'tclean'
        if steptimes.has_key(stepname):
            steptimes[stepname]+=stagedur
        else:
            steplist.append(stepname)
            steptimes[stepname]=stagedur
        stagestr = stepname+' dirty cycle '+str(itercycle)
        stagetime.append(stagedur)
        stagename.append(stagestr)
        print(stagestr+' took '+str(stagedur)+' sec')
        prevTime = currTime
        #
        # ===== Construct the mask for CCBox
        os.system('cp -r '+clnresidual[0]+' '+dirtyresidual[0])
        #
        # construct a PSF with the Gaussian core subtracted (new CJC 10May2016)
        psfresid = clnim+'.psf.subtracted'
        os.system('rm -rf '+psfresid)
        immath(imagename=clnpsf[0],mode='evalexpr',expr='1.0*IM0',outfile=psfresid)
        ia.open(psfresid)
        psfimstat1=ia.statistics()
        blcx=psfimstat1['maxpos'][0]-20
        trcx=psfimstat1['maxpos'][0]+20
        blcy=psfimstat1['maxpos'][1]-20
        trcy=psfimstat1['maxpos'][1]+20
        blctrc=str(blcx)+','+str(blcy)+','+str(trcx)+','+str(trcy)
        clrec=ia.fitcomponents(box=blctrc)
        ia.modify(clrec['results'],subtract=True)
        ia.close()
        psfimstat2=imstat(imagename=psfresid)
        psfmin=max(abs(psfimstat2['min'][0]),psfimstat2['max'][0])
        #
        logstring = 'Using PSF sidelobe level for masking = '+str(psfmin)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        imageimstat=imstat(imagename=dirtyresidual[0])
        immax=imageimstat['max'][0]
        imRMS=imageimstat['rms'][0]
        logstring = "Dirty image RMS = "+str(imRMS)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        pksnr=immax/imRMS
        logstring = 'Dirty image Peak/rms = '+str(pksnr)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        # threshold for initial mask is defined by immax*n*abs(psfmin) (n>=2),
        # unless abs(psfmin)>0.5, in which case do something different...
        #
        if abs(psfmin)<0.5:
            threshfraction=psfmin*int(0.5/psfmin)
        else:
            threshfraction=1.05*psfmin
        #
        thresh1=immax*threshfraction
        logstring = 'Cycle '+str(itercycle)+' initial threshold = '+str(thresh1)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        #
        maskname=clnim+'.cycle'+str(itercycle)
        immath(imagename=dirtyresidual[0],mode='evalexpr',
               expr='iif(IM0>'+str(thresh1)+',1.0,0.0)',
               outfile=maskname+'_mask',stokes='I')
        #
        fwhm1=cellsize*sqrt((psfimstat1['maxpos'][0]-psfimstat1['minpos'][0])**2+(psfimstat1['maxpos'][1]-psfimstat1['minpos'][1])**2)
        fwhm1str=str(fwhm1)+'arcsec'
        logstring = 'Smoothing mask with FWHM='+fwhm1str
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        imsmooth(imagename=maskname+'_mask',kernel='gauss',major=fwhm1str,minor=fwhm1str,
                 pa='0deg',outfile=maskname+'_sm_mask')
        #
        tmpimstat=imstat(imagename=maskname+'_sm_mask')
        maskpk=tmpimstat['max'][0]
        #
        threshmask = maskname+'_sm_thresh_mask'
        immath(imagename=maskname+'_sm_mask',mode='evalexpr',
                  expr='iif(IM0>'+str(maskpk/2.0)+',1.0,0.0)',
                  outfile=threshmask,stokes='I')
        logstring = 'Created smoothed thresholded mask image '+threshmask
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        maskstat=imstat(threshmask)
        logstring = 'Mask image contains '+str(int(maskstat['sum'][0]))+' active pixels'
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        currTime=time.time()
        stagedur = currTime-prevTime
        stepname = 'masking'
        if steptimes.has_key(stepname):
            steptimes[stepname]+=stagedur
        else:
            steplist.append(stepname)
            steptimes[stepname]=stagedur
        stagestr = stepname+' initial'
        stagetime.append(stagedur)
        stagename.append(stagestr)
        print(stagestr+' took '+str(stagedur)+' sec')
        prevTime = currTime
        #
        # =============================
        # Do first real clean iteration
        #
        box_threshold = 3.0*imRMS
        if box_threshold<fld_thresholdJy:
            box_threshold = fld_thresholdJy
        logstring = 'Cleaning submosaic  with mask image '+threshmask+' to '+str(box_threshold)+'Jy'
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        # call tclean with interactive=0 to return iterations
        imresult=tclean(visname,
                        imagename=clnim,
                        field=fieldstr,
                        spw=spwstr,
                        imsize=[fld_size,fld_size],
                        cell=[fld_cell,fld_cell],
                        phasecenter=mycenter,
                        stokes='I',
                        startmodel='',
                        specmode=fld_specmode,
                        reffreq=fld_reffreq,
                        # nchan=fld_nchan,
                        # start=spw_center,
                        # width=spw_width,
                        # interpolation=fld_interp,
                        # restfreq=[fld_reffreq],
                        gridder=fld_gridder,
                        pblimit=fld_pblimit,
                        normtype=fld_normtype,
                        # wbawp=fld_wbawp,
                        # wprojplanes=fld_wprojplanes,
                        # facets=fld_facets,
                        deconvolver=fld_deconvolver,
                        # scales=fld_multiscale,
                        restoringbeam=myrestore,
                        niter=box_niter,
                        threshold=box_threshold,
                        cycleniter=box_cycleniter,
                        cyclefactor=box_cyclefactor,
                        usemask='user',
                        mask=threshmask,
                        interactive=0,
                        weighting=myweight,
                        robust=myrobust,
                        uvtaper=myuvtaper,
                        makeimages='choose',
                        calcres=False,
                        calcpsf=False,
                        savemodel=dosavemodel)
        itercycle+=1
        os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
        if imresult.has_key('iterdone'):
            iterdone=imresult['iterdone']
            logstring = "Imaging for cycle "+str(itercycle)+" completed with "+str(iterdone)+" iterations"
        else:
            iterdone+=1
            logstring = "Imaging for cycle "+str(itercycle)+" completed"
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        oldRMS=imRMS
        imageimstat=imstat(imagename=clnresidual[0])
        immax=imageimstat['max'][0]
        imRMS=imageimstat['rms'][0]
        logstring = "Residual RMS = "+str(imRMS)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        RMSratio=((oldRMS-imRMS)/(oldRMS))
        logstring = "(OLD RMS - NEW RMS) / OLD RMS = "+str(RMSratio)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        pksnr=immax/imRMS
        logstring = 'Peak/rms = '+str(pksnr)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        #
        currTime=time.time()
        stagedur = currTime-prevTime
        stepname = 'tclean'
        if steptimes.has_key(stepname):
            steptimes[stepname]+=stagedur
        else:
            steplist.append(stepname)
            steptimes[stepname]=stagedur
        stagestr = stepname+' boxed cycle '+str(itercycle)
        stagetime.append(stagedur)
        stagename.append(stagestr)
        print(stagestr+' took '+str(stagedur)+' sec')
        prevTime = currTime
        #
        if pksnr>peaksnrlimit and iterdone>0:
            doboxed=True
        else:
            doboxed=False
        #
        # ==========
        # now iterate
        #
        while doboxed:
            if pksnr>peaksnrlimit and iterdone>0 and imRMS>fld_thresholdJy:
                logstring = "RMS ratio too high: more cleaning required ..."
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                maskname=clnim+'.cycle'+str(itercycle)
                os.system('cp -rf '+clnim+'.mask '+maskname+'_oldmask')
                thresh2=immax*threshfraction
                logstring = ' Cycle '+str(itercycle)+' new initial threshold = '+str(thresh2)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                immath(imagename=clnresidual[0],mode='evalexpr',
                       expr='iif(IM0>'+str(thresh2)+',1.0,0.0)',
                       outfile=maskname+'_mask',stokes='I')
                imsmooth(imagename=maskname+'_mask',kernel='gauss',major=fwhm1str,minor=fwhm1str,
                         pa='0deg',outfile=maskname+'_sm_mask')
                tmpimstat=imstat(imagename=maskname+'_sm_mask')
                maskpk=tmpimstat['max'][0]
                immath(imagename=maskname+'_sm_mask',mode='evalexpr',
                       expr='iif(IM0>'+str(maskpk/2.0)+',1.0,0.0)',
                       outfile=maskname+'_sm_thresh_mask',stokes='I')
                threshmask = maskname+'_sm_sum_mask'
                immath(imagename=[maskname+'_sm_thresh_mask',
                                  maskname+'_oldmask'],mode='evalexpr',
                       expr='max(IM0,IM1)',outfile=threshmask,
                       stokes='I')
                #
                logstring = 'Created smoothed thresholded summed mask image '+threshmask
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                #
                maskstat=imstat(threshmask)
                logstring = 'Mask image contains '+str(int(maskstat['sum'][0]))+' active pixels'
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                #
                currTime=time.time()
                stagedur = currTime-prevTime
                stepname = 'masking'
                if steptimes.has_key(stepname):
                    steptimes[stepname]+=stagedur
                else:
                    steplist.append(stepname)
                    steptimes[stepname]=stagedur
                stagestr = stepname+' iter '+str(iterdone)
                stagetime.append(stagedur)
                stagename.append(stagestr)
                print(stagestr+' took '+str(stagedur)+' sec')
                prevTime = currTime
                #
                box_threshold=3.0*imRMS
                if box_threshold<fld_thresholdJy:
                       box_threshold = fld_thresholdJy
                logstring = 'Cleaning submosaic with mask image '+threshmask+' to '+str(box_threshold)+'Jy'
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                # for 4.5.1 have to copy this mask to .mask
                os.system('cp -rf '+threshmask+' '+clnim+'.mask ')
                #
                # call tclean with interactive=0 to return iterations
                imresult=tclean(visname,
                                imagename=clnim,
                                field=fieldstr,
                                spw=spwstr,
                                imsize=[fld_size,fld_size],
                                cell=[fld_cell,fld_cell],
                                phasecenter=mycenter,
                                stokes='I',
                                startmodel='',
                                specmode=fld_specmode,
                                reffreq=fld_reffreq,
                                # nchan=fld_nchan,
                                # start=spw_center,
                                # width=spw_width,
                                # interpolation=fld_interp,
                                # restfreq=[fld_reffreq],
                                gridder=fld_gridder,
                                pblimit=fld_pblimit,
                                normtype=fld_normtype,
                                # wprojplanes=fld_wprojplanes,
                                # facets=fld_facets,
                                deconvolver=fld_deconvolver,
                                # scales=fld_multiscale,
                                restoringbeam=myrestore,
                                niter=box_niter,
                                threshold=box_threshold,
                                cycleniter=box_cycleniter,
                                cyclefactor=box_cyclefactor,
                                usemask='user',
                                # mask=threshmask,
                                mask='',
                                interactive=0,
                                weighting=myweight,
                                robust=myrobust,
                                uvtaper=myuvtaper,
                                makeimages='choose',
                                calcres=False,
                                calcpsf=False,
                                # overwrite=True,
                                savemodel=dosavemodel)
                itercycle+=1
                os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
                if imresult.has_key('iterdone'):
                    iterdone=imresult['iterdone']
                    logstring = 'Imaging for cycle '+str(itercycle)+' completed with '+str(iterdone)+' iterations'
                else:
                    iterdone+=1
                    logstring = 'Imaging for cycle '+str(itercycle)+' completed'
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                #
                oldRMS=imRMS
                imageimstat=imstat(imagename=clnresidual[0])
                immax=imageimstat['max'][0]
                imRMS=imageimstat['rms'][0]
                logstring = "Residual RMS = "+str(imRMS)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                #
                RMSratio=((oldRMS-imRMS)/(oldRMS))
                logstring = "(OLD RMS - NEW RMS) / OLD RMS = "+str(RMSratio)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                #
                pksnr=immax/imRMS
                logstring = 'Peak/rms = '+str(pksnr)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                #
                currTime=time.time()
                stagedur = currTime-prevTime
                stepname = 'tclean'
                if steptimes.has_key(stepname):
                    steptimes[stepname]+=stagedur
                else:
                    steplist.append(stepname)
                    steptimes[stepname]=stagedur
                stagestr = stepname+' boxed cycle '+str(itercycle)
                stagetime.append(stagedur)
                stagename.append(stagestr)
                print(stagestr+' took '+str(stagedur)+' sec')
                prevTime = currTime
                #
                # Possibly terminate boxed cleaning
                if itercycle>=maxboxcycles:
                    doboxed=False
                    logstring = "Reached max number of boxed cycles, terminating boxing"
                    print(logstring)
                    casalog.post(logstring)
                    logbuffer.append(logstring)
            else:
                logstring = "No more boxed cleaning needed, finishing up"
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                doboxed=False
        #
        # Save last mask
        os.system('cp -r '+clnim+'.mask '+clnim+'_lastmask_mask')
    else:
        logstring = 'Will NOT do any autoboxing'
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
    # ==========================================
    # END of (optional) autoboxing code
    # ==========================================
    #
    # ================================
    # Now final clean without box mask
    # ================================
    #
    logstring = 'Final Cleaning submosaic without mask '
    print(logstring)
    casalog.post(logstring)
    logbuffer.append(logstring)
    os.system('rm -rf '+clnim+'.mask ')
    # call tclean with interactive=0 to return iterations
    try:
        imresult=tclean(visname,
                        imagename=clnim,
                        field=fieldstr,
                        spw=spwstr,
                        imsize=[fld_size,fld_size],
                        cell=[fld_cell,fld_cell],
                        phasecenter=mycenter,
                        stokes='I',
                        startmodel='',
                        specmode=fld_specmode,
                        reffreq=fld_reffreq,
                        # nchan=fld_nchan,
                        # start=spw_center,
                        # width=spw_width,
                        # interpolation=fld_interp,
                        # restfreq=[fld_reffreq],
                        gridder=fld_gridder,
                        pblimit=fld_pblimit,
                        normtype=fld_normtype,
                        # wbawp=fld_wbawp,
                        # wprojplanes=fld_wprojplanes,
                        # facets=fld_facets,
                        deconvolver=fld_deconvolver,
                        # scales=fld_multiscale,
                        restoringbeam=myrestore,
                        niter=fld_niter,
                        threshold=fld_threshold_nobox,
                        cycleniter=fld_cycleniter,
                        cyclefactor=fld_cyclefactor,
                        usemask='user',
                        mask='',
                        interactive=0,
                        weighting=myweight,
                        robust=myrobust,
                        uvtaper=myuvtaper,
                        makeimages='auto',
                        # makeimages='choose',
                        calcres=True,
                        calcpsf=True,
                        savemodel=dosavemodel)
        #
        itercycle+=1
        os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
        if imresult.has_key('iterdone'):
            iterdone=imresult['iterdone']
            logstring = 'Final imaging for cycle '+str(itercycle)+' completed with '+str(iterdone)+' iterations'
        else:
            iterdone+=1
            logstring = 'Final imaging for cycle '+str(itercycle)+' completed'
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
    except:
        logstring = 'WARNING: Failed final cleaning submosaic without mask '
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
    #
    currTime=time.time()
    stagedur = currTime-prevTime
    stepname = 'tclean'
    if steptimes.has_key(stepname):
        steptimes[stepname]+=stagedur
    else:
        steplist.append(stepname)
        steptimes[stepname]=stagedur
    stagestr = stepname+' final unboxed cycle '+str(itercycle)
    stagetime.append(stagedur)
    stagename.append(stagestr)
    print(stagestr+' took '+str(stagedur)+' sec')
    prevTime = currTime

#
# Make list of all clean images to mosaic together
# subim if needed
if doimaging or dostats or dosubim:
    # Construct field list and list of subimages
    clndict={}
    num_images = 0
    num_good = 0
    #
    spwstr='TT0'
    clndict[spwstr]={}
    num_images += 1
    clnim = clnname
    if fld_nterms>1:
        # We will use only tt0
        clnimage=clnim+'.image.tt0'
        clnpb = clnim+'.pb'
        clnres = clnim+'.residual.tt0'
    else:
        clnimage=clnim+'.image'
        clnpb = clnim+'.pb'
        clnres = clnim+'.residual'
    if dosubim:
        pbsubim=clnpb+'.subim'
        ressubim = clnres+'.subim'
        (blc_i,blc_j,trc_i,trc_j)=fld_subim.split(',')
        mysize_ra = int(trc_i)-int(blc_i)+1
        mysize_dec = int(trc_j)-int(blc_j)+1
        logstring = 'Creating subimages of size '+str(mysize_ra)+' x '+str(mysize_dec)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        if os.access(clnimage,F_OK):
            # Now subim this image
            clnsubim=clnimage+'.subim'
            imsubimage(clnimage,clnsubim,box=fld_subim)
        else:
            logstring = 'WARNING: could not find '+clnimage
            print(logstring)
            casalog.post(logstring)
            logbuffer.append(logstring)
        if os.access(clnres,F_OK):
            residsubim=clnres+'.subim'
            imsubimage(clnres,residsubim,box=fld_subim)
        else:
            logstring = 'WARNING: could not find '+clnres
            print(logstring)
            casalog.post(logstring)
            logbuffer.append(logstring)
        if os.access(clnpb,F_OK):
            pbsubim=clnpb+'.subim'
            imsubimage(clnpb,pbsubim,box=fld_subim)
        else:
            logstring = 'WARNING: could not find '+clnpb
            print(logstring)
            casalog.post(logstring)
            logbuffer.append(logstring)
    else:
        clnsubim=clnimage
        pbsubim=clnpb
        ressubim = clnres
    if os.access(clnsubim,F_OK):
        num_good += 1
    else:
        logstring = 'WARNING: could not find clean mosaic '
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
    #
    clndict[spwstr]['clnlist'] = [clnsubim]
    clndict[spwstr]['pblist'] = [pbsubim]
    clndict[spwstr]['reslist'] = [ressubim]
    #
    logstring = 'Found '+str(num_good)+' cleaned field planes out of a total '+str(num_images)+' expected'
    print(logstring)
    casalog.post(logstring)
    logbuffer.append(logstring)
    #
    currTime=time.time()
    stagedur = currTime-prevTime
    stepname = 'image verification and subim'
    if steptimes.has_key(stepname):
        steptimes[stepname]+=stagedur
    else:
        steplist.append(stepname)
        steptimes[stepname]=stagedur
    stagestr = stepname
    stagetime.append(stagedur)
    stagename.append(stagestr)
    print(stagestr+' took '+str(stagedur)+' sec')
    prevTime = currTime

# Get stats
statdict={}
statdict['Fields']={}
statdict['Channels']={}
statdict['All']={}
fldsigma={}
rmslist_all=[]
rmslist_region_all=[]
if dostats and num_good>0:
    spwstr='TT0'
    statdict['Channels'][spwstr]={}
    nimg = len(clndict[spwstr]['clnlist'])
    rmslist=[]
    print('Stats for '+spwstr+' :')
    # clean image stats
    ifld=0
    imname = clndict[spwstr]['clnlist'][ifld]
    mystat = imstat(imname,chans='')
    statdict['Channels'][spwstr]['Restored']={}
    statdict['Channels'][spwstr]['Restored']['Stats']=mystat
    cln_max = mystat['max'][0]
    cln_min = mystat['min'][0]
    #
    statstring = 'Restored max = '+str(cln_max)+' min = '+str(cln_min)
    print(statstring)
    statbuffer.append(statstring)
    #
    imname = clndict[spwstr]['reslist'][ifld]
    mystat = imstat(imname,chans='')
    statdict['Channels'][spwstr]['Residual']={}
    statdict['Channels'][spwstr]['Residual']['Stats']=mystat
    res_sigma = mystat['sigma'][0]
    res_max = mystat['max'][0]
    rmslist_all.append(res_sigma)
    #
    statstring = 'Residual sigma = '+str(res_sigma)+' max = '+str(res_max)
    print(statstring)
    statbuffer.append(statstring)
    #
    # PB stats
    imname = clndict[spwstr]['pblist'][ifld]
    mystat = imstat(imname,chans='')
    statdict['Channels'][spwstr]['PB']=mystat
    pb_min = mystat['min'][0]
    pb_max = mystat['max'][0]
    print('PB max = '+str(pb_max)+' min = '+str(pb_min))
    statstring = 'Spw '+spwstr+': PB max = '+str(pb_max)+' min = '+str(pb_min)
    statbuffer.append(statstring)
    #
    # rms statistics in sub-areas of mosaic
    rmslist_region = []
    statdict['Channels'][spwstr]['Residual']['Regions'] = {}
    statdict['Channels'][spwstr]['Residual']['Regions']['Stats'] = {}
    nreg=0
    if dosubim:
        (blc_i,blc_j,trc_i,trc_j)=fld_subim.split(',')
        mysize_ra = int(trc_i)-int(blc_i)+1
        mysize_dec = int(trc_j)-int(blc_j)+1
    else:
        mysize_ra = fld_size
        mysize_dec = fld_size
    #
    for j in range(10):
        jlow = int(mysize_dec*float(j)/10.0)
        jup = int(mysize_dec*float(j+1)/10.0) - 1
        for i in range(10):
            ilow = int(mysize_ra*float(i)/10.0)
            iup = int(mysize_ra*float(i+1)/10.0) - 1
            boxstr = str(ilow)+','+str(jlow)+','+str(iup)+','+str(jup)
            imname = clndict[spwstr]['reslist'][ifld]
            mystat = imstat(imname,box=boxstr,chans='')
            nreg += 1
            statdict['Channels'][spwstr]['Residual']['Regions']['Stats'][nreg] = mystat
            resid_sigma = mystat['sigma'][0]
            rmslist_region_all.append(resid_sigma)
    #
    #
    # Median rms per channel over all spw
    print(' ')
    statbuffer.append(' ')
    print('Statistics for all fields: ')
    statbuffer.append('Statistics for all fields: ')
    rmsarr_all = pl.array(rmslist_all)
    median_sigma_all = pl.median(rmsarr_all)
    max_sigma_all = rmsarr_all.max()
    min_sigma_all = rmsarr_all.min()
    num_sigma_all = len(rmslist_all)
    statdict['All']['Median']={}
    statdict['All']['Median']['Residual']={}
    statdict['All']['Median']['Residual']['sigma']=median_sigma_all
    statdict['All']['Median']['Residual']['max']=max_sigma_all
    statdict['All']['Median']['Residual']['min']=min_sigma_all
    statdict['All']['Median']['Residual']['number']=num_sigma_all
    statstring = 'Residual Mosiac All channels: Median sigma num = '+str(num_sigma_all)+' median = '+str(median_sigma_all)
    print(statstring)
    statbuffer.append(statstring)
    statstring = 'Residual Mosiac All channels: Median sigma max = '+str(max_sigma_all)+' min = '+str(min_sigma_all)
    print(statstring)
    statbuffer.append(statstring)
    #
    rmsarr_region_all = pl.array(rmslist_region_all)
    median_sigma_region_all = pl.median(rmsarr_region_all)
    max_sigma_region_all = rmsarr_region_all.max()
    min_sigma_region_all = rmsarr_region_all.min()
    num_sigma_region_all = len(rmslist_region_all)
    statdict['All']['Median']['Regions'] = {}
    statdict['All']['Median']['Regions']['Residual'] = {}
    statdict['All']['Median']['Regions']['Residual']['number'] = num_sigma_region_all
    statdict['All']['Median']['Regions']['Residual']['sigma'] = median_sigma_region_all
    statdict['All']['Median']['Regions']['Residual']['max']=max_sigma_region_all
    statdict['All']['Median']['Regions']['Residual']['min']=min_sigma_region_all
    statstring = 'Residual Mosaic All Channels: subregions num = '+str(num_sigma_region_all)+' median sigma = '+str(median_sigma_region_all)
    print(statstring)
    statbuffer.append(statstring)
    statstring = 'Residual Mosaic All Channels: subregions max = '+str(max_sigma_region_all)+' min sigma = '+str(min_sigma_region_all)
    print(statstring)
    statbuffer.append(statstring)
    #
    currTime=time.time()
    stagedur = currTime-prevTime
    stepname = 'statistics'
    if steptimes.has_key(stepname):
        steptimes[stepname]+=stagedur
    else:
        steplist.append(stepname)
        steptimes[stepname]=stagedur
    stagestr = stepname
    stagetime.append(stagedur)
    stagename.append(stagestr)
    print(stagestr+' took '+str(stagedur)+' sec')
    prevTime = currTime
    
clnmosfile = clnrestored[0]

#====================================================================
casalog.post('Completed imaging of MS '+splitfile)
#====================================================================

endProc=time.clock()
endTime=time.time()

walltime = (endTime - startTime)
cputime = (endProc - startProc)

import datetime
datestring=datetime.datetime.isoformat(datetime.datetime.today())

myvers = casalog.version()
myuser = os.getenv('USER')
myhost = str( os.getenv('HOST') )
mycwd = os.getcwd()
myos = os.uname()
mypath = os.environ.get('CASAPATH')

#
if myvers.__contains__('#'):
    mybuild = string.split(list.pop(myvers.split('#')),')')[0]
    buildstring = '.r'+str(mybuild)
elif myvers.__contains__(' r'):
    mybuild = string.split(list.pop(myvers.split(' r')),')')[0]
    buildstring = '.r'+str(mybuild)
elif myvers.__contains__('(r'):
    mybuild = string.split(list.pop(myvers.split('(r')),')')[0]
    buildstring = '.r'+str(mybuild)
else:
    mybuild = 'unknown'
    buildstring = ''

# Print to terminal, and also save most things to a logfile
outfile='out.'+scriptprefix+'.'+datestring+buildstring+'.log'
mylogfile = open(outfile,'w')

def lprint(msg, lfile):
    """
    Prints msg to both stdout and lfile.
    """
    print(msg)
    print >>mylogfile, msg
    
lprint(mytext, mylogfile)
lprint('Running '+myvers+' on host '+myhost, mylogfile)
lprint('  at '+datestring, mylogfile)
lprint('  using '+mypath, mylogfile)

lprint('', mylogfile)
lprint('---', mylogfile)
lprint('Script version: '+params['version'], mylogfile)
lprint('User set parameters used in execution:', mylogfile)
lprint('---', mylogfile)

for keys in params['user'].keys():
    lprint('  %s  =  %s ' % ( keys, params['user'][keys] ), mylogfile)

lprint('', mylogfile)

new_regression = {}
# Dataset size info
new_regression['datasize'] = {}
new_regression['datasize']['ms'] = datasize_ms

total = {}
total['wall'] = (endTime - startTime)
total['cpu'] = (endProc - startProc)
total['rate_ms'] = (datasize_ms/(endTime - startTime))

nstages = stagetime.__len__()
timing = {}
timing['total'] = total
timing['nstages'] = nstages
timing['stagename'] = stagename
timing['stagetime'] = stagetime

# Save timing to regression dictionary
new_regression['timing'] = timing

# Logging messages
if len(logbuffer)>0:
    lprint('', mylogfile)
    lprint('* Logging:                                      *', mylogfile)
    print >>mylogfile,''
    for logstring in logbuffer:
        print >>mylogfile,logstring

# Final stats and timing
if dostats:
    lprint('', mylogfile)
    lprint('* Results:                                      *', mylogfile)
    print >>mylogfile,'' 
    for statstring in statbuffer:
        print >>mylogfile,statstring

lprint('', mylogfile)
lprint('********* Benchmarking *************************', mylogfile)
lprint('*                                              *', mylogfile)
lprint('Total wall clock time was: %10.3f ' % total['wall'], mylogfile)
lprint('Total CPU        time was: %10.3f ' % total['cpu'], mylogfile)
lprint('MS  processing rate MB/s was: %8.1f ' % total['rate_ms'], mylogfile)
lprint('', mylogfile)
lprint('* Breakdown:                                   *', mylogfile)

lprint('', mylogfile)
lprint('* Breakdown by stage:                           *', mylogfile)
for i in range(nstages):
    lprint('* %40s * time was: %10.3f ' % (stagename[i],stagetime[i]), mylogfile)

lprint('', mylogfile)
lprint('* Breakdown by steps:                           *', mylogfile)
for stepname in steplist:
    lprint('* %40s * time was: %10.3f ' % (stepname,steptimes[stepname]), mylogfile)

lprint('************************************************', mylogfile)

lprint("", mylogfile)
lprint("Done with " + mydataset, mylogfile)

lprint("", mylogfile)
lprint("Final Joint Mosaic image is " + clnmosfile, mylogfile)

mylogfile.close()
print("Results are in "+outfile)

#====================================================================
# Done
#====================================================================
