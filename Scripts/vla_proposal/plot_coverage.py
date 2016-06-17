""" Overplot the VLASS fields with the LIGO field """

import numpy as np
import matplotlib.pyplot as plt
from read_log import get_field_positions
from plotpilot import get_vlass_fields

ra_ligo, dec_ligo = get_field_positions()
fig = get_vlass_fields()

plt.scatter(ra_ligo, dec_ligo, s=2, c='k')
plt.title("VLASS Pilot Tiles and Search Field for GW151226")
plt.show()
