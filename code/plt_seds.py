

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import cmasher as cmr
import CFigTools.CustomFigure as CF
import gpscssmodels
import pandas as pd
import plotting_functions

# Paths and non-changing variables 
master_data_dir = "/data/raw_data/"
data_dir = "/data/LBA/catalogues/"
lba_targets = ["j0227-0621", "j0322-482", "j2239-451"]

mask_nms, mwa_13, mwa_14, xtra_nms = plotting_functions.make_nm_arrays()
lba_pop_pd = plotting_functions.read_fluxes(mask_nms)

# Created dictionary (src_dict) with array entries for each epoch/chunk of stuff 
src_dict = plotting_functions.make_src_dict(lba_pop_pd)
print(src_dict[lba_targets[0]].keys())

# Plotting of SEDs for each of the srcs 
for i in range(len(lba_targets)):
    plotting_functions.plt_sed("/data/LBA/plots/",src_dict[lba_targets[i]],papersize=True)