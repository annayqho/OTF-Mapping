""" Get the RMS values from clean_log files """

import numpy as np
import matplotlib.pyplot as plt
import glob
from matplotlib import rc
rc('font', family='serif')
rc('text', usetex=True)

def get_vers(prefix):
    vers = ['Run_%sNITER_1000' %prefix,
            'Run_%sNITER_2000' %prefix,
            'Run_%sNITER_3000' %prefix,
            'Run_%sNITER_4000' %prefix,
            'Run_%sNITER_5000' %prefix]
    return vers

def get_rms(rms_blank, vers):
    rms = np.zeros(rms_blank.shape)
    for ii,ver in enumerate(vers):
        fname = direc + ver + "/clean.log"
        a = np.loadtxt(fname, dtype=float)
        rms[:,ii] = a[:,9]
    return rms

direc = "/lustre/kmooley/projects/JAGWAR/GW151226/anna/"
vers_ms = get_vers("Multiscale_")
vers = get_vers("")
vers_wp_ms = get_vers("wpp_Multiscale_")
nvers = len(vers)
nfields = 80
rms_blank = np.zeros((nfields, nvers))
rms = get_rms(rms_blank, vers)
rms_ms = get_rms(rms_blank, vers_ms)
rms_wp_ms = get_rms(rms_blank, vers_wp_ms)
diff = rms - rms_ms

fig = plt.figure()
plt.tick_params(axis="both", labelsize=14)
x = [1000, 2000, 3000, 4000, 5000]

# cols = ['k', 'cyan', 'magenta', 'g', 'orange']
# for ii,val in enumerate(x):
#     plt.hist(diff[:,ii], bins=10, range=(-300,200), 
#             color=cols[ii], label="niter %s" %val,
#             histtype="step")
# 
 
#choose = rms[:,0]!=rms[:,-1]
#choose = rms_ms[:,0] < 280
#choose = diff[:,0] < 0 
choose = rms[:,0] > 0 # all of them
for ii,line in enumerate(rms_ms[choose]): 
    plt.plot(x, line, c='cyan', alpha=0.2)
    plt.plot(x, rms_wp_ms[choose][ii], c='magenta', alpha=0.2)
plt.plot(x, np.median(rms_ms[choose], axis=0), c='cyan', lw=5, label="No Multiscale")
plt.plot(x, np.median(rms_wp_ms[choose], axis=0), c='magenta', lw=5, label="Multiscale")
plt.xlabel("Number of Iterations", fontsize=16)
plt.ylabel(r"RMS ($\mu$Jy)", fontsize=16)
plt.title(r"Multiscale, no wproj, RMS \textless RMS MS (for niter=1000)")
# #ss \, 280 $\mu$Jy", fontsize=16)
plt.legend()
#plt.savefig("multiscale_no_wproj")
plt.show()
# 
