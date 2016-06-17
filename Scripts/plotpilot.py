""" VLA Pilot Survey fields, modified by AYQH from S. Myers"""

import scipy as spy
import numpy as np
import matplotlib.pyplot as plt

def get_vlass_fields():
    # tiles = {}
    tiles = {'COSMOS':{'center':[150.0,4.0,],'width':[10.0,8.0],'color':'b'},
             'Cygnus':{'center':[307.5,40.0,],'width':[12.0,8.0],'color':'g'},
             'Stripe82':{'center':[0.0,0.0,],'width':[120.0,2.6667],'color':'r'},
             'Cepheus':{'center':[345.0,62.0,],'width':[20.0,8.0],'color':'c'},
             'CDFS':{'center':[52.5,-27.0,],'width':[11.25,8.0],'color':'m'},
             'GCen':{'center':[267.0,-29.0,],'width':[11.25,8.0],'color':'y'},
             'SGCap':{'center':[0.0,6.0,],'width':[85.0,16.0],'color':'orange'},
             'NGCap':{'center':[202.5,55.0,],'width':[105.0,10.0],'color':'k'},
        }

    tilelist = tiles.keys()

    markers = {}
    # markers = {'E-N1':{'center':[240.0,55.0],'marker':'x','color':'k'},
    #            }

    markerlist = markers.keys()

    #
    # linestyle = ['_', '-', '--', ':']
    # color = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    #
    fig = plt.figure()
    plt.xlim(0.0,360.0)
    plt.ylim(-40.0,100.0)

    for tile in tilelist:
        blcx = tiles[tile]['center'][0] - 0.5*tiles[tile]['width'][0]
        trcx = tiles[tile]['center'][0] + 0.5*tiles[tile]['width'][0]
        blcy = tiles[tile]['center'][1] - 0.5*tiles[tile]['width'][1]
        trcy = tiles[tile]['center'][1] + 0.5*tiles[tile]['width'][1]
        x = [blcx,trcx,trcx,blcx,blcx]
        y = [blcy,blcy,trcy,trcy,blcy]
        #
        plt.plot(x,y,linestyle='-',color=tiles[tile]['color'],label=tile)
        #
        # check for RA wrap
        if blcx<0.0:
            blcx1 = blcx + 360.0
            trcx1 = trcx + 360.0
            x1 = [blcx1,trcx1,trcx1,blcx1,blcx1]
            plt.plot(x1,y,linestyle='-',color=tiles[tile]['color'])
        elif trcx>360.0:
            blcx1 = blcx - 360.0
            trcx1 = trcx - 360.0
            x1 = [blcx1,trcx1,trcx1,blcx1,blcx1]
            plt.plot(x1,y,linestyle='-',color=tiles[tile]['color'])

    # Now any markers
    for marker in markerlist:
        xm = [markers[marker]['center'][0]]
        ym = [markers[marker]['center'][1]]
        # plt.plot(xm,ym,marker=markers[marker]['marker'],color=markers[marker]['color'],label=marker)
        plt.scatter(xm,ym,marker=markers[marker]['marker'],c=markers[marker]['color'],label=marker)

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    mytitle = 'VLASS Pilot Tile Location'
    plt.title(mytitle)
    # plt.legend(loc='upper left')
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #        ncol=2, mode="expand", borderaxespad=0.)
    plt.legend(loc='upper center', ncol=4, mode="expand")

    return fig
