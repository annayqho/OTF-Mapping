""" Overplot the VLASS fields with the LIGO field """

import numpy as np
import matplotlib.pyplot as plt
from read_log import get_field_positions

ra_ligo, dec_ligo = get_field_positions()
