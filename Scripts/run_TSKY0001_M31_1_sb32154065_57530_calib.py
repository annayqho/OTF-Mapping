""" This is an example calibration script from Steve Myers. 
This calibrates a particular M31 VLASS test block. 
You would change stuff in the header to make it work on another dataset."""

#
# Basic Calibration Script for CASA 4.2 based on old TRSR0015 script
# STM 2014-11-13 version for 13B-370 S300_1 sb28581653_56669
# STM 2015-05-13 coverted for TSKY0001 M31
# STM 2015-06-08 better autoflagging and calibration
# STM 2015-06-23 save final flags
# STM 2015-11-10 test in 4.5.0
# STM 2016-02-23 test in 4.5.1, new directory location
# STM 2016-03-03 test in 4.6.0, better benchmarking
# STM 2016-03-07 use PB 2013 flux scale
# STM 2016-05-24 for TSKY0001.sb32154065.eb32157201.57530
#================================================================================
myscriptvers = '2016-03-07 STM (4.6.x)'
mydataset = 'TSKY0001.sb32154065.eb32157201.57530.47058900463'
scriptprefix='run_TSKY0001_M31_1_sb32154065_57530_calib_pb2013'
prefix='TSKY0001_M31_1_sb32154065_57530'

import time
import os
from os import F_OK

# Set up some parameters for processing
sdmdir = '/lustre/aoc/projects/vlass/smyers/'
# sdmdir = '../Data/' # it's too big to go on my computer
sdmname = mydataset
sdmfile = sdmdir + sdmname
rawmsdir = '/lustre/aoc/projects/vlass/smyers/Run_TSKY0001.sb32154065.eb32157201.57530/'
#rawmsdir = './'
rawmsname = mydataset + '.ms'
rawmsfile = rawmsdir + rawmsname
spwinput = ''
outname = prefix

msfile = outname + '_working.ms'
splitfile = outname + '_calibrated_target.ms'

doimport = True
mytbuff = 0.225
dohanning = True
doinpflagzero = False
doinpshadow = False

doshadow=False
dotfcrop1=True
tfcrop_timecutoff=4.0
tfcrop_freqcutoff=4.0
dorflag1=True
rflag1_timedevscale=4.0
rflag1_freqdevscale=4.0
myrflag1corr='ABS_ALL'
myrflag1ntime=True
dorflag2=True
rflag2_timedevscale=4.0
rflag2_freqdevscale=4.0
myrflag2corr='ABS_ALL'
myrflag2ntime=True
dorflag3=True
rflag3_timedevscale=4.0
rflag3_freqdevscale=4.0
myrflag3corr='ABS_ALL'
myrflag3ntime=True
do_rflag_onepass = True

docleanup=True

#==========================================
# Dataset specific flagging
# Generic flagging from examination of BPASS cal
# flagstr1=["spw='0~15:0~2;61~63'",
#           "spw='0:0~7'",
#           "spw='1:26~34;50~53'",
#           "spw='5:24~26;37~39;51~52'",
#           "spw='6:4~8;16~21;37~45'",
#           "spw='8:0~7;50~61'",
#           "spw='9:3;9~10;16~20;31~46'",
#           # "spw='10:13~16;32~56'",
#           "spw='10'",
#           "spw='11:17~23;29~35;56~61'",
#           "spw='13:53~58'",
#           "spw='14:25~28;52'",
#           "spw='15:3'",
#           "antenna='ea08' spw='8~15'", # something wrong
# 	  "scan='1'"]

flagstr1=["spw='0~15:0~2;61~63'",
          "spw='0:0~7'",
          "spw='8:0~7'",
          "spw='9'",
          "spw='10'",
          "antenna='ea05",
	  "scan='1'"]

flagstr2=[]

# Dataset specific channel range for G0 (avoid flagged regions)
#g0spw = '0~15:30~33'
g0spw = '0:30~33,1:35~38,2~8:30~33,9~10:26~29,11:36~39,12~15:30~33'

# Dataset specific calibrators
fluxfield = 'J0137+3309'
bpassfield = 'J0137+3309'
phasefield = 'J0038+4137'
ampfield = 'J0038+4137'

fluxmodel='3C48_S.im'
#fluxstandard = 'Perley-Butler 2010'
fluxstandard = 'Perley-Butler 2013'

calfields = 'J0137+3309,J0038+4137'
calfieldlist = ['J0137+3309','J0038+4137']

# Dataset specific targets
targetfields = '0*'
targetintent = '*TARGET*'

#==========================================
# 
fluxvalue=1.0
fluxfreq='3.0GHz'
fluxspix=0.0

# calrefant = 'ea27,ea09,ea14'
calrefant = 'ea24,ea25,ea20'

allcalspw = '0~15:3~60'
bpassspw = '0~15:3~60'
dobasebanddelays=True
basebandA = '0~7'
basebandB = '8~15'
basebandchan = '3~60'
calspws = '0~15:8~50'

dofluxscale=True
fluxscalescans = ''

# averaging parms
#average_time = '1s'
#average_chan = 1
average_time = ''
average_chan = 1

#====================================================================
# Save the parameters used 
params = {}
# Version of script
params['version'] = myscriptvers
# User set params
params['user'] = {}
params['user']['sdmfile'] = sdmfile
params['user']['msfile'] = msfile
params['user']['outname'] = outname
params['user']['doimport'] = doimport
params['user']['mytbuff'] = mytbuff
params['user']['dohanning'] = dohanning
params['user']['doinpflagzero'] = doinpflagzero
params['user']['doinpshadow'] = doinpshadow
params['user']['dotfcrop1'] = dotfcrop1
if dotfcrop1:
    params['user']['tfcrop_timecutoff']=tfcrop_timecutoff
    params['user']['tfcrop_freqcutoff']=tfcrop_freqcutoff
params['user']['dorflag1'] = dorflag1
if dorflag1:
    params['user']['rflag1_timedevscale']=rflag1_timedevscale
    params['user']['rflag1_freqdevscale']=rflag1_freqdevscale
    params['user']['myrflag1corr']=myrflag1corr
    params['user']['myrflag1ntime']=myrflag1ntime
params['user']['dorflag2'] = dorflag2
if dorflag2:
    params['user']['rflag2_timedevscale']=rflag2_timedevscale
    params['user']['rflag2_freqdevscale']=rflag2_freqdevscale
    params['user']['myrflag2corr']=myrflag2corr
    params['user']['myrflag2ntime']=myrflag2ntime
if dorflag3:
    params['user']['rflag3_timedevscale']=rflag3_timedevscale
    params['user']['rflag3_freqdevscale']=rflag3_freqdevscale
    params['user']['myrflag3corr']=myrflag3corr
    params['user']['myrflag3ntime']=myrflag3ntime
params['user']['flagstr1'] = flagstr1
params['user']['flagstr2'] = flagstr2
params['user']['fluxfield'] = fluxfield
params['user']['bpassfield'] = bpassfield
params['user']['phasefield'] = phasefield
params['user']['ampfield'] = ampfield
params['user']['targetfields'] = targetfields
params['user']['calrefant'] = calrefant
params['user']['g0spw'] = g0spw
params['user']['bpassspw'] = bpassspw
params['user']['calspws'] = calspws
params['user']['average_time'] = average_time
params['user']['average_chan'] = average_chan
params['user']['fluxmodel'] = fluxmodel
params['user']['fluxstandard'] = fluxstandard
#

# Function to track timing for benchmarks
def stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes):
    currTime=time.time()
    stagedur = currTime-prevTime
    if steptimes.has_key(stepname):
        steptimes[stepname]+=stagedur
    else:
        steplist.append(stepname)
        steptimes[stepname]=stagedur
    stagetag = stepname+': '+stagestr
    stagetime.append(stagedur)
    stagename.append(stagetag)
    print stagetag+' took '+str(stagedur)+' sec'
    casalog.post(stagetag+' took '+str(stagedur)+' sec')
    prevTime = currTime
    return prevTime

stagename = []
stagetime = []
steplist = []
steptimes = {}

#====================================================================
# Start actual processing
#====================================================================
startTime=time.time()
startProc=time.clock()
prevTime = startTime

print 'Starting calibration of MS '+msfile+' using script '+myscriptvers
casalog.post('Starting calibration of MS '+msfile+' using script '+myscriptvers)

if doimport:
    print 'Importing data'
    os.system('rm -rf '+msfile+'*')
    #
    stagestr = 'rm msfile'
    stepname = 'misc'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
    #
    importevla(asdm=sdmfile,vis=msfile,flagbackup=False,online=True,tbuff=mytbuff,flagzero=doinpflagzero,shadow=doinpshadow,applyflags=False)
    #
    stagestr = 'import from sdm'
    stepname = 'importevla'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
    # 
    # get dataset sizes in MB (1024x1024 bytes!)
    f = os.popen('du -ms '+sdmfile)
    fstr = f.readline()
    f.close()
    datasize_sdm = float( fstr.split("\t")[0] )
    print 'SDM is '+str(datasize_sdm)+' MB'
    # 
else:
    if msfile!=rawmsfile:
        if os.access(msfile,F_OK):
            os.system('rm -rf '+msfile+'*')
            stagestr = 'rm msfile'
            stepname = 'misc'
            prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
        #
        # copy from original ms
        print 'Copying '+rawmsfile+' to '+msfile
        os.system('cp -rf '+rawmsfile+' '+msfile)
        #
        stagestr = 'cp msfile'
        stepname = 'misc'
        prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
    else:
        os.system('rm -rf '+msfile+'.flagversions')
        print 'Clearing flags'
        flagdata(vis=msfile,mode='unflag',flagbackup=False)
        #
        stagestr = 'unflag'
        stepname = 'flagdata'
        prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
    
# get dataset sizes in MB (1024x1024 bytes!)
f = os.popen('du -ms '+msfile)
fstr = f.readline()
f.close()
datasize_ms = float( fstr.split("\t")[0] )
print 'MS is '+str(datasize_ms)+' MB'

# Online flags
print 'Applying Online Flags'
if doimport:
    flagcmd(vis=msfile,inpmode='table',action='apply',flagbackup=False)
else:
    flagcmd(vis=msfile,inpmode='xml',tbuff=mytbuff,action='apply',flagbackup=False)
#
stagestr ='apply online flags'
stepname = 'flagcmd'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# PLOT Online flags
print 'Plotting Online Flags'
flagcmd(vis=msfile,inpmode='table',useapplied=True,action='plot',plotfile='plot.'+outname+'.onlineflags.png')
#
stagestr ='plot online flags'
stepname = 'flagcmd'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# Flag for zeroes
if not doinpflagzero:
    print 'Flagging Zeros'
    flagdata(vis=msfile,mode='clip',clipzeros=True,flagbackup=False)
    #
    currTime=time.time()
    stagestr ='flag zeros'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

print 'Listing data'
listobs(msfile,listfile=msfile+'.listobs')
#
stagestr ='list working ms'
stepname = 'listobs'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# Flag for shadowing
if doshadow:
    print 'Flagging Shadowed Antennas'
    flagdata(vis=msfile,mode='shadow',tolerance=2.0,flagbackup=False)
    #
    stagestr ='flag data to shadowed antennas'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# split off calibrators
calfile = outname + '_cal.ms'
print 'Splitting Calibrators to '+calfile
os.system('rm -rf '+calfile+'*')
mstransform(vis=msfile,outputvis=calfile,field=calfields,datacolumn='data',hanning=dohanning)
if dohanning:
    stagestr =' split cals w/hanning'
else:
    stagestr =' split cals'
stepname = 'mstransform'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
                
# initial Manual flags
if flagstr1!='':
    print 'Applying manual flagging'
    flagdata(vis=calfile,mode='list',inpfile=flagstr1,action='apply',flagbackup=True)
    #
    stagestr ='manual flags initial'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# Initial TFCROP for RFI
if dotfcrop1:
    print 'Applying TFCROP'
    flagdata(vis=calfile,mode='tfcrop',correlation='',timecutoff=tfcrop_timecutoff,freqcutoff=tfcrop_freqcutoff,timefit='line',freqfit='poly',maxnpieces=9,flagdimension='freqtime',extendflags=False,action='apply',savepars=True,outfile=calfile+'.tfcrop1parms')
    #
    stagestr ='tfcrop1 cals'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
    #
    print 'Running FLAGDATA extend'
    flagdata(vis=calfile,correlation='',mode='extend',extendpols=True,growtime=50.0,growaround=True,flagnearfreq=False,action='apply')
    #
    stagestr ='tfcrop1 cals extend'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# clean up any previous calibration
os.system('rm -rf cal.'+outname+'.*')
#
stagestr = 'rm caltables'
stepname = 'misc'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# ***SetJy

print 'Running SetJy'
if fluxmodel!='':
    print 'Will setjy using standard '+fluxstandard+' for '+fluxfield
    setjy(vis=calfile,field=fluxfield,standard=fluxstandard,model=fluxmodel,scalebychan=True,usescratch=False)
else:
    setjy(vis=calfile,field=fluxfield,standard='manual',
          fluxdensity=[fluxvalue,0,0,0],spix=fluxspix,reffreq=fluxfreq,
          scalebychan=True,usescratch=False)
#
stagestr ='setjy flux cal'
stepname = 'setjy'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# ===============================================================================
# ***Start Calibration

# prior calibration : antpos
print 'Generating Antenna Position Corrections'
antpos_name='cal.'+outname+'.antpos'
gencal(calfile,caltable=antpos_name,caltype='antpos')
if os.access(antpos_name,F_OK):
    # antpos table exists
    print 'Successfully Generated antpos correction table '+antpos_name
else:
    # no antpos corrections
    gencal(calfile,caltable=antpos_name,caltype='ph',parameter=[0.0])
    print 'No antenna position corrections, using zero phases '
stagestr ='gencal antpos'
stepname = 'calibration'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# prior calibration : gaincurves
#print 'Generating Gaincurves'
#gc_name='cal.'+outname+'.gc'
#gencal(calfile,caltable=gc_name,caltype='gc')

# intial phasecal on the calibrators
print 'Initial Phasecal before Delay/Bandpass'
g0_name='cal.'+outname+'.G0'
gaincal(vis=calfile,caltable=g0_name,gaintable=[antpos_name],
        field=calfields,spw=g0spw,gaintype='G',
        refant=calrefant,calmode='p',solint='int',minsnr=3)

#plotcal(caltable=g0_name,xaxis='time',yaxis='phase',iteration='antenna', \
#         plotrange=[-1,-1,-180,180])
#
stagestr ='cal init phase'
stepname = 'calibration'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# delay
print 'Generating Antenna Delays'
k0_name='cal.'+outname+'.K0'
if dobasebanddelays:
    # NEW FEATURE: multi-band delays using combine='spw' per baseband
    #delayspw = basebandA + ':' + basebandchan
    delayspw = basebandA
    gaincal(vis=calfile,caltable=k0_name,gaintable=[antpos_name,g0_name], \
            field=bpassfield,spw=delayspw,gaintype='K', 
            gainfield=['',bpassfield],\
            refant=calrefant,combine='scan,spw',solint='inf',minsnr=3)
    #delayspw = basebandB + ':' + basebandchan
    delayspw = basebandB
    gaincal(vis=calfile,caltable=k0_name,gaintable=[antpos_name,g0_name], \
            field=bpassfield,spw=delayspw,gaintype='K', \
            gainfield=['',bpassfield],\
            refant=calrefant,combine='scan,spw',solint='inf',minsnr=3,append=True)
    delayspwmap = [0,0,0,0,0,0,0,0,8,8,8,8,8,8,8,8]
else:
    delayspw = bpassspw
    gaincal(vis=calfile,caltable=k0_name,gaintable=[antpos_name,g0_name], \
            field=bpassfield,spw=delayspw,gaintype='K', \
            gainfield=['',bpassfield],\
            refant=calrefant,combine='scan',solint='inf',minsnr=3)
    delayspwmap = []

#plotcal(k0_name,xaxis='antenna')
#
stagestr ='cal delay'
stepname = 'calibration'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# bandpass
print 'Generating Bandpasses'
b0_name='cal.'+outname+'.B0'
bandpass(vis=calfile,caltable=b0_name,gaintable=[antpos_name,g0_name,k0_name], \
         field=bpassfield, spw=bpassspw, \
         gainfield=['',bpassfield,''],interp=['nearest','nearest','nearest'], \
         spwmap=[[],[],delayspwmap], \
         bandtype='B',combine='scan',solint='inf',refant=calrefant,solnorm=False)

#plotcal(caltable=b0_name,xaxis='freq',yaxis='amp',iteration='antenna')
#
#plotcal(caltable=b0_name,xaxis='freq',yaxis='phase',iteration='antenna', \
#        plotrange=[-1,-1,-180,180])
#
stagestr ='bandpass cal'
stepname = 'calibration'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# RFLAG for RFI and bandpass calibrated data
if dorflag1:
    #
    for fld in calfieldlist:
        applycal(vis=calfile,gaintable=[antpos_name,g0_name,k0_name,b0_name], \
                     spwmap=[[],[],delayspwmap,[]], \
                     field=fld,gainfield=['',fld,'',''], \
                     interp=['nearest','nearest','nearest','nearest'],\
                     calwt=False,parang=False,applymode='calflag',flagbackup=False)
    #
    if do_rflag_onepass:
        #
        print 'Calculating and Applying RFLAG'
        flagdata(vis=calfile,mode='rflag',datacolumn='corrected',correlation=myrflag1corr,timedevscale=rflag1_timedevscale,freqdevscale=rflag1_freqdevscale,savepars=True,outfile=calfile+'.rflag1parms')
    else:
        tdevfile1 = 'cal.'+outname+'.rflag1.timedev.txt'
        fdevfile1 = 'cal.'+outname+'.rflag1.freqdev.txt'
        os.system('rm -rf '+tdevfile1)
        os.system('rm -rf '+fdevfile1)
        #
        print 'Calculating RFLAG'
        flagdata(vis=calfile,mode='rflag',datacolumn='corrected',correlation=myrflag1corr,timedev=tdevfile1,freqdev=fdevfile1,timedevscale=rflag1_timedevscale,freqdevscale=rflag1_freqdevscale,action='calculate')
        #
        print 'Applying RFLAG'
        flagdata(vis=calfile,mode='rflag',datacolumn='corrected',correlation=myrflag1corr,timedev=tdevfile1,freqdev=fdevfile1,timedevscale=rflag1_timedevscale,freqdevscale=rflag1_freqdevscale,action='apply',savepars=True,outfile=calfile+'.rflag1parms')
    #
    print 'Extending RFLAG'
    flagdata(vis=calfile,correlation='',mode='extend',extendpols=True,growtime=50.0,growaround=True,flagneartime=myrflag1ntime,flagnearfreq=False,action='apply')
    #
    stagestr ='rflag1'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
    #
    # Need to redo bandpass calibration using this flagging
    #
    b0_name='cal.'+outname+'.B1'
    bandpass(vis=calfile,caltable=b0_name,gaintable=[antpos_name,g0_name,k0_name], \
                 field=bpassfield, spw=bpassspw, \
                 gainfield=['',bpassfield,''],interp=['nearest','nearest','nearest'], \
                 spwmap=[[],[],delayspwmap], \
                 bandtype='B',combine='scan',solint='inf',refant=calrefant,solnorm=False)
    #
    stagestr ='bandpass redo'
    stepname = 'calibration'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# gain calibrators
# phase per scan
print 'Generating Phase Gains per scan'
g1scan_name='cal.'+outname+'.G1scan'
gaincal(vis=calfile,caltable=g1scan_name,gaintable=[antpos_name,k0_name,b0_name], \
            spwmap=[[],delayspwmap,[]], \
            field=calfields,spw=calspws,refant=calrefant,solnorm=F, \
            solint='inf',gaintype='G',calmode='p')
#
stagestr ='gaincal scan phase'
stepname = 'calibration'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# phase per integration
print 'Generating Phase Gains per int'
g1int_name='cal.'+outname+'.G1int'
gaincal(vis=calfile,caltable=g1int_name,gaintable=[antpos_name,k0_name,b0_name], \
            spwmap=[[],delayspwmap,[]], \
            field=calfields,spw=calspws,refant=calrefant,solnorm=F, \
            solint='int',gaintype='G',calmode='p')
#
stagestr ='gaincal int phase'
stepname = 'calibration'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# amp per scan
print 'Generating Amplitude Gains'
g2scan_name='cal.'+outname+'.G2scan'
for fld in calfieldlist:
    gaincal(vis=calfile, caltable=g2scan_name,gaintable=[antpos_name,k0_name,b0_name,g1int_name], \
                spwmap=[[],delayspwmap,[],[]], \
                field=fld, gainfield=['','','',fld], refant=calrefant, spw=calspws, \
                solnorm=F,solint='inf',gaintype='G',calmode='a',append=True)
#
stagestr ='gaincal scan amp'
stepname = 'calibration'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# Fluxscale if needed
if dofluxscale:
    print 'Generating Flux-Scaled Amplitude Gains'
    # Solve for incremental fluxscale table
    f3scan_name='cal.'+outname+'.F3scan'
    fluxscale(calfile,caltable=g2scan_name,fluxtable=f3scan_name,
              reference=fluxfield,transfer=calfields,
              incremental=True,listfile=f3scan_name+'.listfile',fitorder=1)
    #
    # plotcal(caltable=g1scan_name,xaxis='time',yaxis='phase',iteration='antenna',plotrange=[-1,-1,-180,180])
    # plotcal(caltable=f3scan_name,xaxis='time',yaxis='amp',iteration='antenna')
    #
    stagestr ='fluxscale'
    stepname = 'calibration'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# apply to calibrators
print 'Applying Calibration to Calibrators'
if dofluxscale:
    for fld in calfieldlist:
        if fld==fluxfield:
            applycal(vis=calfile,gaintable=[antpos_name,k0_name,b0_name,g1int_name,g2scan_name], \
                         spwmap=[[],delayspwmap,[],[],[]], \
                         field=fluxfield,gainfield=['','','',fld,fld], \
#                         interp=['nearest','nearest','nearest','nearest','nearest'], \
                         interp=['nearest','nearest','nearest','linear','nearest'], \
                         calwt=False,parang=False,applymode='calflag',flagbackup=True)
        else:
            applycal(vis=calfile,gaintable=[antpos_name,k0_name,b0_name,g1int_name,g2scan_name,f3scan_name], \
                         spwmap=[[],delayspwmap,[],[],[],[]], \
                         field=fld,gainfield=['','','',fld,fld,fld], \
#                         interp=['nearest','nearest','nearest','nearest','nearest','nearest'],\
                         interp=['nearest','nearest','nearest','linear','nearest','nearest'],\
                         calwt=False,parang=False,applymode='calflag',flagbackup=False)
else:
    for fld in calfieldlist:
        applycal(vis=calfile,gaintable=[antpos_name,k0_name,b0_name,g1int_name,g2scan_name], \
                     spwmap=[[],delayspwmap,[],[],[]], \
                     field=fld,gainfield=['','','',fld,fld], \
#                     interp=['nearest','nearest','nearest','nearest','nearest'],\
                     interp=['nearest','nearest','nearest','linear','nearest'],\
                     calwt=False,parang=False,applymode='calflag',flagbackup=False)
#
stagestr ='applycal cals'
stepname = 'applycal'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
#
# Final RFLAG for RFI
if dorflag2:
    #
    if do_rflag_onepass:
        #
        print 'Calculating and Applying RFLAG2 to Cal Data'
        flagdata(vis=calfile,field=calfields,mode='rflag',datacolumn='corrected',correlation=myrflag2corr,timedevscale=rflag2_timedevscale,freqdevscale=rflag2_freqdevscale,savepars=True,outfile=calfile+'.rflag2calparms')
        #
    else:
        tdevfile2 = 'cal.'+outname+'.rflag2.timedev.txt'
        fdevfile2 = 'cal.'+outname+'.rflag2.freqdev.txt'
        os.system('rm -rf '+tdevfile2)
        os.system('rm -rf '+fdevfile2)
        #
        print 'Calculating RFLAG2 cals'
        flagdata(vis=calfile,field=calfields,mode='rflag',datacolumn='corrected',correlation=myrflag2corr,timedev=tdevfile2,freqdev=fdevfile2,timedevscale=rflag2_timedevscale,freqdevscale=rflag2_freqdevscale,action='calculate')
        #
        print 'Applying RFLAG2 to Cal Data'
        flagdata(vis=calfile,field=calfields,mode='rflag',datacolumn='corrected',correlation=myrflag2corr,timedev=tdevfile2,freqdev=fdevfile2,timedevscale=rflag2_timedevscale,freqdevscale=rflag2_freqdevscale,action='apply',savepars=True,outfile=calfile+'.rflag2calparms')
    #
    print 'Extending RFLAG2 cals'
    flagdata(vis=calfile,field=calfields,correlation='',mode='extend',extendpols=True,growtime=50.0,growaround=True,flagneartime=myrflag2ntime,flagnearfreq=False,action='apply')
    #
    stagestr ='rflag2 cals'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
#
# Save final flags
flagmanager(vis=calfile,mode='save',versionname='final',comment='final flags')
stagestr ='save final cal flags'
stepname = 'flagmanager'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# split off target data, parallel hands only
print 'Splitting target data RR,LL only'
os.system('rm -rf '+splitfile+'*')
mstransform(vis=msfile,outputvis=splitfile,datacolumn='data',field=targetfields,intent=targetintent,correlation='RR,LL',hanning=dohanning)
# 
if dohanning:
    stagestr ='split target w/hanning'
else:
    stagestr ='split target'
stepname = 'mstransform'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
#
# get dataset sizes in MB (1024x1024 bytes!)
f = os.popen('du -ms '+splitfile)
fstr = f.readline()
f.close()
datasize_final = float( fstr.split("\t")[0] )
print 'Split Target MS is '+str(datasize_final)+' MB'
#
print 'Listing Split Target MS'
listobs(splitfile,listfile=splitfile+'.listobs')
# 
stagestr ='listobs final target split'
stepname = 'listobs'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
#
print 'Applying Calibration to Target'
if dofluxscale:
    applycal(vis=splitfile,
             gaintable=[antpos_name,k0_name,b0_name,g1scan_name,g2scan_name,f3scan_name], \
                 spwmap=[[],delayspwmap,[],[],[],[]], \
                 field='',gainfield=['','','',phasefield,ampfield,ampfield], \
                 interp=['nearest','nearest','nearest','linear','linear','linear'], \
                 calwt=False,parang=False,applymode='calflagstrict',flagbackup=True)
else:
    applycal(vis=splitfile,
             gaintable=[antpos_name,k0_name,b0_name,g1scan_name,g2scan_name], \
                 spwmap=[[],delayspwmap,[],[],[],[]], \
                 field='',gainfield=['','','',phasefield,ampfield], \
                 interp=['nearest','nearest','nearest','linear','linear'], \
                 calwt=False,parang=False,applymode='calflagstrict',flagbackup=True)
#
stagestr ='applycal target'
stepname = 'applycal'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# Final Manual flags on target
flagstr = flagstr1 + flagstr2
if len(flagstr)>0:
    print 'Final Manual Flagging of Target Data'
    flagdata(vis=splitfile,mode='list',inpfile=flagstr,action='apply',flagbackup=True)
    #
    stagestr ='flagdata manual final'
    stepname = 'flagdata'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
    #
    # Final RFLAG for RFI
    if dorflag3:
        print 'Final RFLAG of Target Data'
        rflag3datacolumn='corrected'
        #
        if do_rflag_onepass:
            print 'Calculating and Applying RFLAG3'
            flagdata(vis=splitfile,field=targetfields,intent=targetintent,mode='rflag',datacolumn=rflag3datacolumn,correlation=myrflag3corr,timedevscale=rflag3_timedevscale,freqdevscale=rflag3_freqdevscale)
        else:
            tdevfile3 = 'cal.'+outname+'.rflag3.timedev.txt'
            fdevfile3 = 'cal.'+outname+'.rflag3.freqdev.txt'
            os.system('rm -rf '+tdevfile3)
            os.system('rm -rf '+fdevfile3)
            #
            print 'Calculating RFLAG3'
            flagdata(vis=splitfile,field=targetfields,intent=targetintent,mode='rflag',datacolumn=rflag3datacolumn,correlation=myrflag3corr,timedev=tdevfile3,freqdev=fdevfile3,timedevscale=rflag3_timedevscale,freqdevscale=rflag3_freqdevscale,action='calculate')
            #
            print 'Applying RFLAG3 to Target Data'
            flagdata(vis=splitfile,field=targetfields,intent=targetintent,mode='rflag',datacolumn=rflag3datacolumn,correlation=myrflag3corr,timedev=tdevfile3,freqdev=fdevfile3,timedevscale=rflag3_timedevscale,freqdevscale=rflag3_freqdevscale,action='apply',savepars=True,outfile=calfile+'.rflag2parms')
        #
        print 'Extending RFLAG3'
        flagdata(vis=splitfile,field=targetfields,intent=targetintent,correlation='',mode='extend',extendpols=True,growtime=50.0,growaround=True,flagneartime=myrflag3ntime,flagnearfreq=False,action='apply')
        #
        stagestr ='rflag3 final target'
        stepname = 'flagdata'
        prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)
#
# Save final flags
flagmanager(vis=splitfile,mode='save',versionname='final',comment='final flags')
stagestr ='save final target flags'
stepname = 'flagmanager'
prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

# Final clean up
if docleanup:
    os.system('rm -rf '+msfile+'*')
    stagestr ='rm working msfiles'
    stepname = 'misc'
    prevTime = stagemark(prevTime,stepname,stagestr,stagetime,stagename,steplist,steptimes)

#====================================================================
casalog.post('Completed calibration of MS '+msfile)
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

# Get version string
versstring = ''
myvershash = myvers.split(' ')
if myvershash.count('Version')>0:
    indx = myvershash.index('Version')
    if indx<(len(myvershash)-1):
        versstring = '_Ver'+myvershash[indx+1]

# Get build number - Deal with all the previous possible formats
if myvers.__contains__('#'):
    mybuild = string.split(list.pop(myvers.split('#')),')')[0]
    buildstring = versstring+'.r'+str(mybuild)
elif myvers.__contains__('(r'):
    mybuild = string.split(list.pop(myvers.split('(r')),')')[0]
    buildstring = versstring+'.r'+str(mybuild)
elif myvers.__contains__(' r'):
    mybuild = string.split(list.pop(myvers.split(' r')),')')[0]
    buildstring = versstring+'.r'+str(mybuild)
else:
    mybuild = 'unknown'
    buildstring = versstring

# Print to terminal, and also save most things to a logfile
outfile='out.'+prefix+'.'+datestring+buildstring+'.log'
mylogfile = open(outfile,'w')

def lprint(msg, lfile):
    """
    Prints msg to both stdout and lfile.
    """
    print msg
    print >>mylogfile, msg
    
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
new_regression['datasize']['raw'] = datasize_sdm
new_regression['datasize']['ms'] = datasize_ms

total = {}
total['wall'] = (endTime - startTime)
total['cpu'] = (endProc - startProc)
total['rate_raw'] = (datasize_sdm/(endTime - startTime))
total['rate_ms'] = (datasize_ms/(endTime - startTime))

nstages = stagetime.__len__()
timing = {}
timing['total'] = total
timing['nstages'] = nstages
timing['stagename'] = stagename
timing['stagetime'] = stagetime

# Save timing to regression dictionary
new_regression['timing'] = timing

# Final stats and timing

lprint('', mylogfile)
lprint('********* Benchmarking *************************', mylogfile)
lprint('*', mylogfile)
lprint('Total wall clock time was: %10.3f ' % total['wall'], mylogfile)
lprint('Total CPU        time was: %10.3f ' % total['cpu'], mylogfile)
lprint('Raw processing rate MB/s was: %8.1f ' % total['rate_raw'], mylogfile)
lprint('MS  processing rate MB/s was: %8.1f ' % total['rate_ms'], mylogfile)

lprint('*', mylogfile)
lprint('* Breakdown by stage: ', mylogfile)
for i in range(nstages):
    lprint('* %40s * time was: %10.3f ' % (stagename[i],stagetime[i]), mylogfile)

lprint('*', mylogfile)
lprint('* Breakdown by steps: ', mylogfile)
for stepname in steplist:
    lprint('* %40s * time was: %10.3f ' % (stepname,steptimes[stepname]), mylogfile)

lprint('************************************************', mylogfile)

lprint("", mylogfile)
lprint("Done with " + mydataset, mylogfile)

mylogfile.close()
print "Results are in "+outfile

#====================================================================
# Done
#====================================================================
