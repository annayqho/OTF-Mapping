""" Pull out the RA and Dec of all fields """

import numpy as np
import matplotlib.pyplot as plt
from math import floor


def deg_to_ra(deg):
    """ Convert deg to RA in hh:mm:ss.ss """
    hours = deg * 24 / 360
    hh = floor(hours)
    minutes = (hours-hh)*60
    mm = floor(minutes)
    seconds = (minutes-mm)*60
    return "%s:%s:%s" %(hh,mm,seconds)


def ra_to_deg(ra):
    """ Convert RA to degrees.
    Input format: XX:XX:XX.XXXX """
    hh,mm,ss = ra.split(':')
    hours = float(hh) + float(mm)/60 + float(ss)/3600
    deg = hours*360/24
    return deg

def decl_to_deg(decl):
    """ Convert declination to degrees.
    Input format: +XX.XX.XX.XXXXX """
    dd,mm,ss,frac = decl.split('.')
    deg = float(dd) + float(mm)/60 + (float(ss+'.'+frac))/3600
    return deg

def get_field_positions():
    direc = "/Users/annaho/Dropbox/Research/VLA-GW"
    inputf = open("%s/casapy-20160615-210844.log" %direc, "r")
    log = inputf.readlines()
    inputf.close()

    sizes = np.array([len(val) for val in log])

    # When I plot the sizes, I see a jump at sizes[5510] and sizes[10988].
    # I assume that everything in between is the fields...

    # OK, the field descriptions start at log[5513] 
    # and there are 5474 scans, so it makes sense that the last one
    # would be at 10987

    field_data_raw = np.array([val.split(' ') for val in log[5513:10988]])
    field_data = np.array([list(filter(None, val)) for val in field_data_raw])

    # the values are 
    # DATE, TIME+stuff, ID, Code, Name, RA, Dec, Epoch, SrcID, nRows

    ra = field_data[:,5]
    dec = field_data[:,6]
    ra_deg = np.array([ra_to_deg(val) for val in ra])
    dec_deg = np.array([decl_to_deg(val) for val in dec])
    # plot(ra_deg, dec_deg)

    return ra_deg, dec_deg

def plot(ra_deg, dec_deg):
    plt.scatter(ra_deg, dec_deg, s=2, c='k')
    plt.gca().invert_xaxis()
    plt.xlabel("RA (deg)", fontsize=16)
    plt.ylabel("Dec (deg)", fontsize=16)
    plt.title("11 Feb: Region 2, Epoch 1", fontsize=16)
    #plt.show()
    plt.savefig("11feb_field_map.png")

if __name__=="__main__":
    get_field_positions()

