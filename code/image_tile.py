import copy
import time
import os
import subprocess
from os import F_OK
from getfieldcone import *


class Image(object):
    """ Image a tile of VLA data """


    def __init__(self, ms_file, ms_datacolumn, target_intent, target_fields,
            script_file, tile_center_epo, tile_center_ra, tile_center_dec, 
            L_tile_arcsec, tile_pixelsize, tile_padding_arcsec, imaging_dir,
            spwstr, dobox, clear_pointing, doccbox, mask, maxboxcycles,
            box_niter, box_cyclefactor, box_cycleniter, peaksnrlimit, 
            fld_wprojplanes, fld_facets, fld_gridder, fld_pblimit,
            fld_normtype, myweight, myrmode, myrobust, dosavemodel,
            use_multiscale, fld_deconvolver, fld_specmode, fld_reffreq,
            fld_nterms, uv_taper, fld_multiscale, fld_threshold, use_restore,
            fld_niter, fld_cyclefaactor, fld_cycleniter, doimaging, docleanup,
            dousescratch, dostats, parallel):
        self.wall_time = None
        self.cpu_time = None
        self.current_stage = None
        self.logbuffer = []
        self.stagename = []
        self.stagetime = []
        self.statbuffer = []
        intro = 'VLA OTFM Imaging Script, 2016-07-26 AYQH'
        add_logstring(intro)
        self.ms_file = ms_file
        set_ms_file(self.ms_file)
        self.ms_datacolumn = ms_datacolumn
        set_ms_datacolumn(self.ms_datacolumn)
        self.clear_pointing = clear_pointing
        set_pointing_table(clear_pointing)
        self.target_intent = target_intent
        self.script_file = script_file
        self.fld_size = set_fld_size(
                tile_pixelsize, L_tile_arcsec, tile_padding_arcsec)
        tile_center = set_center(
                tile_center_epo, tile_center_ra, tile_center_dec)
        self.fldnos = get_fields(L_tile_arcsec, ms_file, tile_center, 
                                 target_fields, fld_specmode)
        self.imaging_dir = set_imaging_dir(imaging_dir)
        self.spwstr = spwstr
        self.dobox = dobox 
        self.doccbox = doccbox
        self.mask = mask
        self.maxboxcycles = maxboxcycles
        self.box_niter = box_niter
        self.box_cyclefactor = box_cyclefactor
        self.box_cycleniter = box_cycleniter
        self.peaksnrlimit = peaksnrlimit, 
        self.fld_wprojplanes = fld_wprojplanes
        self.fld_facets = fld_facets
        self.fld_gridder = fld_gridder
        self.fld_pblimit = fld_pblimit
        self.fld_normtype = fld_normtype
        self.myweight = myweight
        self.myrmode = myrmode
        self.myrobust = myrobust
        self.dosavemodel = dosavemodel
        self.use_multiscale = use_multiscale
        self.fld_deconvolver = fld_deconvolver
        self.fld_specmode = fld_specmode
        self.fld_reffreq = fld_reffreq
        self.fld_nterms = fld_nterms
        self.uv_taper = uv_taper
        self.fld_multiscale = fld_multiscale
        self.fld_threshold = fld_threshold
        self.myrestore = set_myrestore(userestore)
        self.fld_niter = fld_niter
        self.fld_cyclefaactor = fld_cyclefactor
        self.fld_cycleniter = fld_cycleniter
        self.doimaging = doimaging
        self.docleanup = docleanup
        self.dousescratch = dousescratch
        self.dostats = dostats
        self.parallel = parallel


    def add_logstring(logstring):
        """ Print, add to casalog, append to the logbuffer 
        
        Parameters
        ----------
        logstring (str): string to append
        """
        print(logstring)
        casalog.post(logstring)
        self.logbuffer.append(logstring)


    def lprint(msg, lfile):
        """ Prints msg to both stdout and lfile. """
        print(msg)
        print >>mylogfile, msg


    def find_sources():
        """ Use PyBDSM to identify point sources.
        Ideally, would feed it a fits image, have it spit out
        an image in CASA or .fits format.
        bash script could call CASA that does some bit of cleaning.
        Then run PyBDSM then run the rest of CASA. """
        print("Coming soon...")


    def set_center(epo, ra, dec):
        """ Set submosaic center & record it
        
        Parameters
        ----------
        epo: (str) the epoch, e.g. J2000
        ra: (str) the tile center RA, format 'XXhXXmXXs'
        dec: (str) the tile center Dec, format 'XXdXXmXXs'
        """
        tile_center = me.direction(epo, ra, dec)
        logstring = 'Tile Center: % (RA, Dec) = (%s, %s)' %(epo, ra, dec)
        add_logtring(logstring)
        return tile_center


    def set_ms_file(filename):
        """ Check that file exists, add to log """
        if os.access(filename, F_OK):
            logstring = 'Found calibrated ms '+splitfile
            # get dataset sizes in MB (1024x1024 bytes!)
            fstr = subprocess.check_output(["du", "-ms", filename])
            datasize_ms = fstr.split("\t")[0]
            logstring = 'The size of the ms file is %s MB' %datasize_ms
            add_logstring(logstring)
        else:
            logstring = 'ERROR: could not find '+splitfile
        add_logstring(logstring)


    def set_ms_datacolumn(datacolumn):
        """ Data column in logstring """
        logstring = 'Using datacolumn ' + datacolumn
        add_logstring(logstring)


    def set_pointing_table(clear_pointing):
        """ Pointing table in logstring """
        # Option to clear POINTING table before imaging
        if clear_pointing:
            logstring = 'Will clear MS POINTING table(s)'
        else:
            logstring = 'MS POINTING table(s) will be used if contain data'
        add_logstring(logstring)


    def set_imaging_dir(imaging_dir):
        """ Set up the name of the output image directory """
        logstring = "Output images go in dir %s" %imaging_dir
        add_logstring(logstring)
        if os.path.exists(imaging_dir):
            logstring = 'Using existing image output directory: ' + imaging_dir
        else:
            logstring = 'Creating directory ' + imaging_dir
            os.makedirs(imaging_dir)
        add_logstring(logstring)
        return imaging_dir


    def set_thresholds(self):
        """ Set thresholds for processing """
        fld_threshold_q = qa.quantity(self.fld_threshold)
        fld_threshold_qJy = qa.convert(fld_threshold_q,'Jy')
        fld_thresholdJy = fld_threshold_qJy['value']

        fld_threshold_nobox_q = qa.quantity(self.fld_threshold_nobox)
        fld_threshold_nobox_qJy = qa.convert(fld_threshold_nobox_q,'Jy')
        fld_thresholdJy_nobox = fld_threshold_nobox_qJy['value']


    def set_fld_size(tile_pixelsize, L_tile_arcsec, padding_arcsec):
        """ Set the size of the field, including padding

        Parameters
        ----------
        tile_pixelsize: size of a pixel in arcsecs
        L_tile_arcsec: length of tile in arcsec
        tile_padding_arcsec: amount of padding in arcsec

        Returns
        -------
        fld_size: length of the field, tile + padding
        """
        L_subtile_pixels = int(L_tile_arcsec/tile_pixelsize)
        padding = int(padding_arcsec/tile_pixelsize)
        fld_size = L_subtile_pixels + padding
        logstring = 'Using field image size %i with cell size %s arcsec' % (
                fld_size,tile_pixelsize)
        add_logstring(logstring)
        logbuffer.append(logstring)
        return fld_size


    def set_myrestore(userestore):
        """ Set parameters for restoring beam 

        Parameters
        ----------
        fld_bmaj: 
        fld_bmin:
        fld_bpa:

        Returns
        -------
        myrestore (list): beam parameters
        """
        if use_restore == '':
            myrestore = []
        else:
            fld_bmaj = use_restore
            fld_bmin = use_restore
            fld_bpa = '0deg'
            myrestore = [fld_bmaj, fld_bmin, fld_bpa]
        return myrestore


    def get_fields(L_tile_arcsec, ms_file, tile_center, tgt_fields, specmode):
        """ Generate a list of fields to be cleaned 

        Parameters
        ----------
        L_tile_arcsec
        ms_file
        tile_center
        tgt_fields: target fields
        specmode: fld_specmode
        
        Returns
        -------
        fldnos (list): fields to be cleaned
        """
        fldnos = []

        if dobox:
            beamsearchradius_arcsec = 1000.0
            mydist = 0.5 * L_tile_arcsec + beamsearchradius_arcsec
            logstring = 'Selecting fields w/in box of length %s arcsec' %mydist
            add_logstring(logstring)

            if target_fields==[]:
                logstring = 'No target field preference provided'
                fldnos = getfieldirbox(
                        ms_file, distance=mydist, center_dir=tile_center)
            else:
                logstring = 'Matching field names using regex string(s): %s' \
                        %(str(target_fields))
                fldnos = getfieldirbox(
                        ms_file, distance=mydist, center_dir=tile_center, 
                        matchregex=target_fields)
            add_logstring(logstring)
            logstring = 'Will image a total of %s fields using specmode %s' \
                    %(str(len(fldnos)), str(fld_specmode))
            add_logstring(logstring)
            logstring = 'Will image fields = %s' %str(fldnos)
            add_logstring(logstring)
        return fldnos


    def initialize_stage(self, stage_name):
        """ 
        Initialize a stage of the process 
        Start the clock, print a message
        """
        if self.current_stage != None:
            print("There is another stage already running!")
        else:
            self.current_stage = stage_name
            logstring = "Initializing stage %s" %stage_name
            add_logstring(logstring)
            self.start_clock()


    def end_stage(self, stage_name):
        """
        End a stage of the process
        End the clock, measure length of time
        """
        if self.current_stage == None:
            print("There is no stage being run!")
        else:
            wall_time = self.wall_time
            cpu_time = self.cpu_time
            self.check_clock
            total_wall_time = self.wall_time - wall_time
            total_cpu_time = self.cpu_time - cpu_time
            self.stagename.append(stagename)
            self.stagetime.append(total_wall_time)
            self.current_stage = None
            logstring = "Ending stage %s. Total wall time was %s, Total CPU \
                    time was %s." %(stage_name, total_wall_time, total_cpu_time)
            add_logstring(logstring)


    def check_clock(self):
        """ Check the clock """
        self.wall_time = time.time()
        self.cpu_time = time.clock() 


    def initialize_imaging(self):
        """ Initialize imaging """
        logstring = 'Starting imaging of %s using %s' %(
                self.ms_file, self.script_file)
        add_logstring(logstring)


    def create_ms_working_copy(self):
        """
        Create a working copy of the ms & make a listing
        """
        print("Creating a working copy of the MS...")
        sdmfile = self.ms_file.split('.ms')[0]
        workfile = sdmfile + '_working_copy.ms'
        visname = self.imaging_dir + '/' + workfile
        print('Splitting ' + self.ms_file + ' to ' + visname)
        os.system('rm -rf ' + visname + '*')
        if splitdatacolumn=='auto':
            # detect whether there is CORRECTED_DATA & split only this
            tb.open(self.ms_file)
            split_colnames = tb.colnames()
            tb.close()
            if 'CORRECTED_DATA' in split_colnames:
                mydatacolumn = 'corrected'
            else:
                mydatacolumn = 'data'
        else:
            mydatacolumn = splitdatacolumn
        stepname = 'split'
        fldstrs = ', '.join(self.fldnos)
        mstransform(
                self.ms_file, visname, field=fldstrs,
                intent=self.target_intent, datacolumn=mydatacolumn)
        # make a listing
        listobs(visname,listfile=visname+'.listobs')
        print("Working copy of MS created")


    def clear_pointing(ms_file):
        """ Clear POINTING table in ms file """
        tb.open(ms_file + '/POINTING', nomodify = False)
        a = tb.rownumbers()
        tb.removerows(a)
        tb.close()
        logstring = 'Cleared MS POINTING table(s) for '+ ms_file
        add_logstring(logstring)
        

    def make_dirty_image(self):
        """ Make a dirty image, including the primary beam .pb image """
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
                   gridder=fld_gridder,
                   pblimit=fld_pblimit,
                   normtype=fld_normtype,
                   deconvolver=fld_deconvolver,
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
                   savemodel=dosavemodel,
                   parallel=parallel)
            logstring = "Created dirty mosaic"
        except:
            logstring = 'WARNING: Failed creating dirty mosaic'
        add_logstring(logstring)
        
        # Save parameters from this run
        os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
        
# Make standard (untapered) image set
startimagingTime = currTime
if doimaging:
    spwstr = imaging_spwstr
    foundbox=False
    
    # =====================================================
    # Make a joint MFS mosaic of all fields in the split MS
    # =====================================================
    fldlist = range(len(fldnos))
    fieldstr = ''
    clnname = imaging_dir+'/img.'+mydataset+postfix+'.clean'
    clnim = clnname
    dirtyname = imaging_dir+'/img.'+mydataset+postfix+'.dirty'
    dirtyim = dirtyname

    # clean up previous images
    os.system('rm -rf '+clnim+'.*')
    os.system('rm -rf '+dirtyim+'.*')
    
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
    
    iterdone=0
    itercycle=0


    # ===============================================
    # (Optionally) use an input mask, say from PyBDSM
    # ===============================================
    if mask != '':
        # ===== Have a look at the provided mask
        threshmask = mask
        maskstat=imstat(threshmask)
        npix = int(maskstat['sum'][0])
        logstring = 'Mask image contains %s active pixels' %str(npix)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
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
        
        # =============================
        # Do first real clean iteration
        # =============================
        box_threshold = fld_threshold_nobox
        logstring = 'Cleaning submosaic with mask image %s to %s Jy' \
                %(threshmask, str(box_threshold))
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        # call tclean with interactive=0 to return iterations
        imresult=tclean(visname,
                        imagename=clnim,
                        field=fieldstr,
                        spw=spwstr,
                        imsize=fld_size,
                        cell=fld_cell,
                        phasecenter=mycenter,
                        stokes='I',
                        startmodel='',
                        specmode=fld_specmode,
                        reffreq=fld_reffreq,
                        gridder=fld_gridder,
                        pblimit=fld_pblimit,
                        normtype=fld_normtype,
                        deconvolver=fld_deconvolver,
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
                        #makeimages='choose',
                        #makeimages='auto',
                        calcres=True,
                        calcpsf=True,
                        savemodel=dosavemodel,
                        #savemodel='none',
                        parallel=parallel)
        itercycle+=1
        os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
        if imresult.has_key('iterdone'):
            iterdone=imresult['iterdone']
            logstring = "Imaging for cycle %s completed with %s iterations" \
                    %(str(itercycle), str(iterdone))
        else:
            iterdone+=1
            logstring = "Imaging for cycle %s completed" %str(itercycle)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        imageimstat=imstat(imagename=clnresidual[0])
        immax=imageimstat['max'][0]
        imRMS=imageimstat['rms'][0]
        logstring = "Residual RMS = "+str(imRMS)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        pksnr=immax/imRMS
        logstring = 'Peak/rms = '+str(pksnr)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
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
        
        if pksnr > peaksnrlimit and iterdone > 0:
            doboxed = True
        else:
            doboxed = False

        # Calculate threshfraction
        # construct a PSF with the Gaussian core subtracted 
        psfresid = clnim+'.psf.subtracted'
        os.system('rm -rf '+psfresid)
        immath(imagename=clnpsf[0],mode='evalexpr',expr='1.0*IM0',outfile=psfresid,stokes='I')
        ia.open(psfresid)
        psfimstat1=ia.statistics()
        # maxpos and minpos are the coordinates of the
        # max and min pixel values respectively
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
        
        logstring = 'Using PSF sidelobe level for masking = '+str(psfmin)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        # threshold for initial mask is defined by immax*n*abs(psfmin) (n>=2),
        # unless abs(psfmin)>0.5, in which case do something different...
        if abs(psfmin)<0.5:
            threshfraction=psfmin*int(0.5/psfmin)
        else:
            threshfraction=1.05*psfmin
 
        
        # ==========
        # now iterate
        # ==========
        while doboxed:
            if pksnr>peaksnrlimit and iterdone>0 and imRMS>fld_thresholdJy:
                logstring = "RMS ratio too high: more cleaning required ..."
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                maskname=clnim+'.cycle'+str(itercycle)
                os.system('cp -rf '+clnim+'.mask '+maskname+'_oldmask')
                thresh2=immax*threshfraction
                logstring = 'Cycle %s new initial threshold = %s' \
                        %(str(itercycle), str(thresh2))
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                immath(
                        imagename=clnresidual[0],mode='evalexpr',
                        expr='iif(IM0>'+str(thresh2)+',1.0,0.0)',
                        outfile=maskname+'_mask',stokes='I')
                imsmooth(
                        imagename=maskname+'_mask',kernel='gauss',
                        major=fwhm1str,minor=fwhm1str,
                        pa='0deg',outfile=maskname+'_sm_mask')
                tmpimstat=imstat(imagename=maskname+'_sm_mask')
                maskpk=tmpimstat['max'][0]
                immath(
                        imagename=maskname+'_sm_mask',mode='evalexpr',
                        expr='iif(IM0>'+str(maskpk/2.0)+',1.0,0.0)',
                        outfile=maskname+'_sm_thresh_mask',stokes='I')
                threshmask = maskname+'_sm_sum_mask'
                immath(imagename=[maskname+'_sm_thresh_mask',
                                  maskname+'_oldmask'],mode='evalexpr',
                       expr='max(IM0,IM1)',outfile=threshmask,
                       stokes='I')
                
                logstring = 'Created smoothed thresholded summed mask image '+threshmask
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                maskstat=imstat(threshmask)
                logstring = 'Mask image contains %s active pixels' \
                        %str(int(maskstat['sum'][0]))
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
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
                
                box_threshold=3.0*imRMS
                if box_threshold<fld_thresholdJy:
                       box_threshold = fld_thresholdJy
                logstring = 'Cleaning submosaic with mask image %s to %s Jy' \
                        %(threshmask, str(box_threshold))
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)

                # for 4.5.1 have to copy this mask to .mask
                os.system('cp -rf '+threshmask+' '+clnim+'.mask ')
                
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
                                gridder=fld_gridder,
                                pblimit=fld_pblimit,
                                normtype=fld_normtype,
                                deconvolver=fld_deconvolver,
                                restoringbeam=myrestore,
                                niter=box_niter,
                                threshold=box_threshold,
                                cycleniter=box_cycleniter,
                                cyclefactor=box_cyclefactor,
                                usemask='user',
                                mask='',
                                interactive=0,
                                weighting=myweight,
                                robust=myrobust,
                                uvtaper=myuvtaper,
                                makeimages='choose',
                                calcres=False,
                                calcpsf=False,
                                savemodel=dosavemodel,
                                parallel=parallel)
                itercycle+=1
                os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
                if imresult.has_key('iterdone'):
                    iterdone=imresult['iterdone']
                    logstring = 'Imaging for cycle %s completed with %s iterations' \
                            %(str(itercycle), str(iterdone))
                else:
                    iterdone+=1
                    logstring = 'Imaging for cycle '+str(itercycle)+' completed'
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                oldRMS=imRMS
                imageimstat=imstat(imagename=clnresidual[0])
                immax=imageimstat['max'][0]
                imRMS=imageimstat['rms'][0]
                logstring = "Residual RMS = "+str(imRMS)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                RMSratio=((oldRMS-imRMS)/(oldRMS))
                logstring = "(OLD RMS - NEW RMS) / OLD RMS = "+str(RMSratio)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                pksnr=immax/imRMS
                logstring = 'Peak/rms = '+str(pksnr)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
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
        
        # Save last mask
        os.system('cp -r '+clnim+'.mask '+clnim+'_lastmask_mask')

    else:
        logstring = 'Will NOT use an input mask'
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
    # ==========================================
    # END of (optional) input mask code
    # ==========================================
    
    # ==========================================
    # (Optionally) do autoboxing
    # ==========================================
    if doccbox:
        # ===== Construct the mask for CCBox
        os.system('cp -r '+clnresidual[0]+' '+dirtyresidual[0])
        
        # construct a PSF with the Gaussian core subtracted 
        psfresid = clnim+'.psf.subtracted'
        os.system('rm -rf '+psfresid)
        immath(imagename=clnpsf[0],mode='evalexpr',expr='1.0*IM0',outfile=psfresid,stokes='I')
        ia.open(psfresid)
        psfimstat1=ia.statistics()
        # maxpos and minpos are the coordinates of the
        # max and min pixel values respectively
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
        
        logstring = 'Using PSF sidelobe level for masking = '+str(psfmin)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        imageimstat=imstat(imagename=dirtyresidual[0])
        immax=imageimstat['max'][0]
        imRMS=imageimstat['rms'][0]
        logstring = "Dirty image RMS = "+str(imRMS)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        pksnr=immax/imRMS
        logstring = 'Dirty image Peak/rms = '+str(pksnr)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        # threshold for initial mask is defined by immax*n*abs(psfmin) (n>=2),
        # unless abs(psfmin)>0.5, in which case do something different...
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
        
        maskname=clnim+'.cycle'+str(itercycle)
        immath(imagename=dirtyresidual[0],mode='evalexpr',
               expr='iif(IM0>'+str(thresh1)+',1.0,0.0)',
               outfile=maskname+'_mask',stokes='I')
       
        dist1 = psfimstat1['maxpos'][0]-psfimstat1['minpos'][0]
        dist2 = psfimstat1['maxpos'][1]-psfimstat1['minpos'][1]
        fwhm1 = cellsize * sqrt(dist1**2+dist2**2)
        fwhm1str=str(fwhm1)+'arcsec'
        logstring = 'Smoothing mask with FWHM='+fwhm1str
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        imsmooth(imagename=maskname+'_mask',kernel='gauss',major=fwhm1str,minor=fwhm1str,
                 pa='0deg',outfile=maskname+'_sm_mask')

        tmpimstat=imstat(imagename=maskname+'_sm_mask')
        maskpk=tmpimstat['max'][0]
        threshmask = maskname+'_sm_thresh_mask'
        immath(imagename=maskname+'_sm_mask',mode='evalexpr',
                  expr='iif(IM0>'+str(maskpk/2.0)+',1.0,0.0)',
                  outfile=threshmask,stokes='I')
        logstring = 'Created smoothed thresholded mask image '+threshmask
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        maskstat=imstat(threshmask)
        npix = int(maskstat['sum'][0])
        logstring = 'Mask image contains %s active pixels' %str(npix)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
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
        
        # =============================
        # Do first real clean iteration
        # =============================
        box_threshold = 3.0*imRMS
        if box_threshold<fld_thresholdJy:
            box_threshold = fld_thresholdJy
        logstring = 'Cleaning submosaic with mask image %s to %s Jy' \
                %(threshmask, str(box_threshold))
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
                        gridder=fld_gridder,
                        pblimit=fld_pblimit,
                        normtype=fld_normtype,
                        deconvolver=fld_deconvolver,
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
                        savemodel=dosavemodel,
                        parallel=parallel)
        itercycle+=1
        os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
        if imresult.has_key('iterdone'):
            iterdone=imresult['iterdone']
            logstring = "Imaging for cycle %s completed with %s iterations" \
                    %(str(itercycle), str(iterdone))
        else:
            iterdone+=1
            logstring = "Imaging for cycle %s completed" %str(itercycle)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        oldRMS=imRMS
        imageimstat=imstat(imagename=clnresidual[0])
        immax=imageimstat['max'][0]
        imRMS=imageimstat['rms'][0]
        logstring = "Residual RMS = "+str(imRMS)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        RMSratio=((oldRMS-imRMS)/(oldRMS))
        logstring = "(OLD RMS - NEW RMS) / OLD RMS = "+str(RMSratio)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
        pksnr=immax/imRMS
        logstring = 'Peak/rms = '+str(pksnr)
        print(logstring)
        casalog.post(logstring)
        logbuffer.append(logstring)
        
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
        
        if pksnr>peaksnrlimit and iterdone>0:
            doboxed=True
        else:
            doboxed=False
        
        # ==========
        # now iterate
        # ==========
        while doboxed:
            if pksnr>peaksnrlimit and iterdone>0 and imRMS>fld_thresholdJy:
                logstring = "RMS ratio too high: more cleaning required ..."
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                maskname=clnim+'.cycle'+str(itercycle)
                os.system('cp -rf '+clnim+'.mask '+maskname+'_oldmask')
                thresh2=immax*threshfraction
                logstring = 'Cycle %s new initial threshold = %s' \
                        %(str(itercycle), str(thresh2))
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                immath(
                        imagename=clnresidual[0],mode='evalexpr',
                        expr='iif(IM0>'+str(thresh2)+',1.0,0.0)',
                        outfile=maskname+'_mask',stokes='I')
                imsmooth(
                        imagename=maskname+'_mask',kernel='gauss',
                        major=fwhm1str,minor=fwhm1str,
                        pa='0deg',outfile=maskname+'_sm_mask')
                tmpimstat=imstat(imagename=maskname+'_sm_mask')
                maskpk=tmpimstat['max'][0]
                immath(
                        imagename=maskname+'_sm_mask',mode='evalexpr',
                        expr='iif(IM0>'+str(maskpk/2.0)+',1.0,0.0)',
                        outfile=maskname+'_sm_thresh_mask',stokes='I')
                threshmask = maskname+'_sm_sum_mask'
                immath(imagename=[maskname+'_sm_thresh_mask',
                                  maskname+'_oldmask'],mode='evalexpr',
                       expr='max(IM0,IM1)',outfile=threshmask,
                       stokes='I')
                
                logstring = 'Created smoothed thresholded summed mask image '+threshmask
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                maskstat=imstat(threshmask)
                logstring = 'Mask image contains %s active pixels' \
                        %str(int(maskstat['sum'][0]))
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
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
                
                box_threshold=3.0*imRMS
                if box_threshold<fld_thresholdJy:
                       box_threshold = fld_thresholdJy
                logstring = 'Cleaning submosaic with mask image %s to %s Jy' \
                        %(threshmask, str(box_threshold))
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)

                # for 4.5.1 have to copy this mask to .mask
                os.system('cp -rf '+threshmask+' '+clnim+'.mask ')
                
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
                                gridder=fld_gridder,
                                pblimit=fld_pblimit,
                                normtype=fld_normtype,
                                deconvolver=fld_deconvolver,
                                restoringbeam=myrestore,
                                niter=box_niter,
                                threshold=box_threshold,
                                cycleniter=box_cycleniter,
                                cyclefactor=box_cyclefactor,
                                usemask='user',
                                mask='',
                                interactive=0,
                                weighting=myweight,
                                robust=myrobust,
                                uvtaper=myuvtaper,
                                makeimages='choose',
                                calcres=False,
                                calcpsf=False,
                                savemodel=dosavemodel,
                                parallel=parallel)
                itercycle+=1
                os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
                if imresult.has_key('iterdone'):
                    iterdone=imresult['iterdone']
                    logstring = 'Imaging for cycle %s completed with %s iterations' \
                            %(str(itercycle), str(iterdone))
                else:
                    iterdone+=1
                    logstring = 'Imaging for cycle '+str(itercycle)+' completed'
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                oldRMS=imRMS
                imageimstat=imstat(imagename=clnresidual[0])
                immax=imageimstat['max'][0]
                imRMS=imageimstat['rms'][0]
                logstring = "Residual RMS = "+str(imRMS)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                RMSratio=((oldRMS-imRMS)/(oldRMS))
                logstring = "(OLD RMS - NEW RMS) / OLD RMS = "+str(RMSratio)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
                pksnr=immax/imRMS
                logstring = 'Peak/rms = '+str(pksnr)
                print(logstring)
                casalog.post(logstring)
                logbuffer.append(logstring)
                
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
    
    # ================================
    # Now final clean without box mask
    # ================================
    logstring = 'Final Cleaning submosaic without mask '
    # logstring = 'Cleaning submosaic with mask %s' %mask
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
                        gridder=fld_gridder,
                        pblimit=fld_pblimit,
                        normtype=fld_normtype,
                        deconvolver=fld_deconvolver,
                        restoringbeam=myrestore,
                        niter=fld_niter,
                        threshold=fld_threshold_nobox,
                        cycleniter=fld_cycleniter,
                        cyclefactor=fld_cyclefactor,
                        usemask='user',
                        mask='',
                        #mask = mask,
                        interactive=0,
                        weighting=myweight,
                        robust=myrobust,
                        uvtaper=myuvtaper,
                        makeimages='auto',
                        calcres=True,
                        calcpsf=True,
                        savemodel=dosavemodel,
                        parallel=parallel)
        
        itercycle+=1
        os.system('cp tclean.last '+clnim+'_tclean_'+str(itercycle)+'.last')
        if imresult.has_key('iterdone'):
            iterdone=imresult['iterdone']
            logstring = 'Final imaging for cycle %s completed with %s iterations' \
                    %(str(itercycle), str(iterdone))
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
    
    clndict[spwstr]['clnlist'] = [clnsubim]
    clndict[spwstr]['pblist'] = [pbsubim]
    clndict[spwstr]['reslist'] = [ressubim]
    
    logstring = 'Found %s cleaned field planes out of a total %s expected' \
            %(str(num_good), str(num_images))
    print(logstring)
    casalog.post(logstring)
    logbuffer.append(logstring)
    
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
    
    statstring = 'Restored max = '+str(cln_max)+' min = '+str(cln_min)
    print(statstring)
    statbuffer.append(statstring)
    
    imname = clndict[spwstr]['reslist'][ifld]
    mystat = imstat(imname,chans='')
    statdict['Channels'][spwstr]['Residual']={}
    statdict['Channels'][spwstr]['Residual']['Stats']=mystat
    res_sigma = mystat['sigma'][0]
    res_max = mystat['max'][0]
    rmslist_all.append(res_sigma)
    
    statstring = 'Residual sigma = '+str(res_sigma)+' max = '+str(res_max)
    print(statstring)
    statbuffer.append(statstring)
    
    # PB stats
    imname = clndict[spwstr]['pblist'][ifld]
    mystat = imstat(imname,chans='')
    statdict['Channels'][spwstr]['PB']=mystat
    pb_min = mystat['min'][0]
    pb_max = mystat['max'][0]
    print('PB max = '+str(pb_max)+' min = '+str(pb_min))
    statstring = 'Spw '+spwstr+': PB max = '+str(pb_max)+' min = '+str(pb_min)
    statbuffer.append(statstring)
    
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
    statstring = 'Residual Mosiac All channels: \
            Median sigma num = '+str(num_sigma_all)+\
            ' median = '+str(median_sigma_all)
    print(statstring)
    statbuffer.append(statstring)
    statstring = 'Residual Mosiac All channels: Median sigma max = '+str(max_sigma_all)+\
            ' min = '+str(min_sigma_all)
    print(statstring)
    statbuffer.append(statstring)
    
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
    statstring = 'Residual Mosaic All Channels: subregions num = '\
            +str(num_sigma_region_all)+\
            ' median sigma = '+str(median_sigma_region_all)
    print(statstring)
    statbuffer.append(statstring)
    statstring = 'Residual Mosaic All Channels: subregions max = '\
            +str(max_sigma_region_all)\
            +' min sigma = '+str(min_sigma_region_all)
    print(statstring)
    statbuffer.append(statstring)
    
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
