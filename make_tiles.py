""" Given field locations in (ra, dec) generate 1deg x 1deg tiles """

import numpy as np
import matplotlib.pyplot as plt
from itertools import permutations, repeat
import itertools
from read_log import get_field_positions


def get_tiles(ra, dec):
    # identify the bounds
    padding = 0
    min_ra = min(ra) - padding
    max_ra = max(ra) + padding
    min_dec = min(dec) - padding
    max_dec = max(dec) + padding
    
    # divide up the bounds into a grid 1 deg in size
    ra_grid = np.arange(min_ra, max_ra+1, step=1)
    dec_grid = np.arange(min_dec, max_dec+1, step=1)

    # define each box in the grid by its center
    ra_centers = []
    dec_centers = []
    for ii in range(len(ra_grid)-1):
        ra_centers.append((ra_grid[ii+1]+ra_grid[ii])/2)
    for ii in range(len(dec_grid)-1):
        dec_centers.append((dec_grid[ii+1]+dec_grid[ii])/2)
    centers = np.array(list(itertools.product(ra_centers, dec_centers)))

    # for each center, check if any fields are in its bounds
    keep = np.zeros(len(centers), dtype=bool)
    for ii,center in enumerate(centers):
        x, y = center
        nfields = sum(np.logical_and(np.abs(ra-x)<=0.5, np.abs(dec-y)<=0.5))
        if nfields > 0: keep[ii] = True

    return centers[keep], ra_grid, dec_grid


ra, dec = get_field_positions()

# define rough individual regions
reg1_dec = np.logical_and(dec > 34, dec < 40)
reg1 = np.logical_and(reg1_dec, ra < 62)
reg2_dec = np.logical_and(dec > 46, dec < 52)
reg2 = np.logical_and(reg2_dec, ra < 75)

# generate tiles
centers, ra_grid, dec_grid = get_tiles(ra[reg1], dec[reg1])
plt.scatter(ra[reg1], dec[reg1], c='k', s=2, alpha=0.5)
for ra_line in ra_grid:
    plt.axvline(ra_line, c='r', lw=2, linestyle='--', alpha=0.5)
for dec_line in dec_grid:
    plt.axhline(dec_line, c='r', lw=2, linestyle='--', alpha=0.5)
for center in centers:
    x, y = center
    plt.scatter(x, y, c='k', s=50, marker='x')
plt.xlabel("RA (deg)", fontsize=16)
plt.ylabel("Dec (deg)", fontsize=16)
plt.title("Tiling", fontsize=20)
plt.savefig("tiling_sample.png")
plt.show()

# centers, ra_grid, dec_grid = get_tiles(ra[reg2], dec[reg2])
# plt.scatter(ra, dec, c='k', s=2)
# for ra_line in ra_grid:
#     plt.axvline(ra_line, c='r', lw=2, linestyle='--', alpha=0.5)
# for dec_line in dec_grid:
#     plt.axhline(dec_line, c='r', lw=2, linestyle='--', alpha=0.5)
# for center in centers:
#     x, y = center
#     plt.scatter(x, y, c='k', s=50, marker='x')
# plt.show()
# 
# plt.scatter(ra[reg2], dec[reg2])
# plt.show()


