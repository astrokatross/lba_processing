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
# print(src_dict[lba_targets[0]].keys())

# Plotting of SEDs for each of the srcs
# TODO: Add the parameters to the src dictionary 
# TODO: Add function for plot (multiple function options for each region)
for i in range(len(lba_targets)):
    tar_dict = src_dict[lba_targets[i]]
    try:
        model_params = []
        for j in range(len(tar_dict["flx_lba_2ghz"])):
            model_params.append(plotting_functions.fit_lba(
                [tar_dict["flx_lba_2ghz"][j], tar_dict["flx_lba_8ghz"][j]]))            
            print(f"powlaw params for region {j}: {model_params[-1]}")
        
        model_params.append(plotting_functions.fit_lba(
            [np.sum(tar_dict["flx_lba_2ghz"]), np.sum(tar_dict["flx_lba_8ghz"])]
        ))
        print(f"integrated flux of LBA fitting parameters: {model_params[-1]}")
        model_params.append(plotting_functions.fit_lba(
            tar_dict["flx_atca_20"],
            fit_function=gpscssmodels.powlawbreak,
            frequency=[
                1.33,
                1.407,
                1.638,
                1.869,
                2.1,
                2.331,
                2.562,
                2.793,
                4.71,
                5.090,
                5.500,
                5.910,
                6.320,
                8.732,
                9.245,
                9.758,
                10.269,
            ],
        ))
        print(f"ATCA SED powlaw params: {model_params[-1]}")
    except TypeError:
        print("Only one region to fit")
        try: 
            model_params = []
            params = plotting_functions.fit_lba(
                [tar_dict["flx_lba_2ghz"], tar_dict["flx_lba_8ghz"]]
            )
            model_params.append(params)
            print(params)
        except ValueError:
            print("NaN or inf for LBA fluxes, check I have something to fit")
    plotting_functions.plt_sed("/data/LBA/plots/", tar_dict, model_params, papersize=False)

    
