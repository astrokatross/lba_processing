

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
frequencies = ["2GHz", "8GHz"]
lba_targets = ["j0227-0621", "j0322-482", "j2239-451"]
mwa_targets = ["GLEAM J022744-062106", "GLEAM J032237-482010", "GLEAM J223933-451414"]

mask_nms, mwa_13, mwa_14, xtra_nms = plotting_functions.make_nm_arrays()
lba_pop_pd = plotting_functions.read_fluxes(mask_nms)
src_dict = {}

# Created dictionary (src_dict) with array entries for each epoch/chunk of stuff 
# Array names: flx_mwa_13, flx_mwa_14, err_mwa_13, err_mwa_14, flx_supp, flx_lba_2ghz, flx_lba_8ghz, err_lba_2ghz, err_lba_8ghz, flx_atca_20, err_atca_20
for i in range(len(mwa_targets)):
    # Adding the mwa 13 and 14 fluxes 
    # TO DO: add the 2020 monitoring epochs too 
    tar = lba_targets[i]
    src = {}
    src["MWA Name"] = mwa_targets[i]
    src_13_fluxes = np.squeeze(lba_pop_pd.loc[lba_pop_pd["Name"]==mwa_targets[i]][mwa_13[0]].to_numpy())
    src_13_errs = np.squeeze(lba_pop_pd.loc[lba_pop_pd["Name"]==mwa_targets[i]][mwa_13[1]].to_numpy())
    src_14_fluxes = np.squeeze(lba_pop_pd.loc[lba_pop_pd["Name"]==mwa_targets[i]][mwa_14[0]].to_numpy())
    src_14_errs = np.squeeze(lba_pop_pd.loc[lba_pop_pd["Name"]==mwa_targets[i]][mwa_14[1]].to_numpy())
    src_14_errs = np.insert(src_14_errs, 0, [np.nan, np.nan, np.nan, np.nan])
    src_14_fluxes = np.insert(src_14_fluxes, 0, [np.nan, np.nan, np.nan, np.nan])
    src["flx_mwa_13"] = src_13_fluxes
    src["flx_mwa_14"] = src_14_fluxes
    src["err_mwa_13"] = src_13_errs
    src["err_mwa_14"] = src_14_errs

    # Adding the fluxes of supplementary surveys
    extra_fluxes = np.squeeze(lba_pop_pd.loc[lba_pop_pd["Name"]==mwa_targets[i]][xtra_nms].to_numpy())
    extra_fluxes[0] = extra_fluxes[0]*0.001
    extra_fluxes[2] = extra_fluxes[2]*0.001
    extra_fluxes[3] = extra_fluxes[3]*0.001
    extra_fluxes[5] = extra_fluxes[5]*0.001
    extra_fluxes[6] = extra_fluxes[6]*0.001
    extra_fluxes[7] = extra_fluxes[7]*0.001
    extra_fluxes[8] = extra_fluxes[8]*0.001
    src["flx_supp"] = extra_fluxes

    # Adding the LBA fluxes for each region 
    try:
        t = Table.read(f"{data_dir}{tar}_2GHz.fits")
        lba_2ghz_flux = np.squeeze(np.array(t["int_flux"]))
        lba_2ghz_rms = np.squeeze(np.array(t["local_rms"]))
    except FileNotFoundError:
        print(f"No catalogue for {tar} 2GHz")
        lba_2ghz_flux = np.nan
        lba_2ghz_rms = np.nan
    try:
        t = Table.read(f"{data_dir}{tar}_8GHz.fits")
        lba_8ghz_flux = np.squeeze(np.array(t["int_flux"]))
        lba_8ghz_rms = np.squeeze(np.array(t["local_rms"]))
    except FileNotFoundError:
        print(f"No catalogue for {tar} 8GHz")
        lba_8ghz_flux = np.nan
        lba_8ghz_rms = np.nan
    src["flx_lba_2ghz"] = lba_2ghz_flux
    src["flx_lba_8ghz"] = lba_8ghz_flux
    src["err_lba_2ghz"] = lba_2ghz_rms
    src["err_lba_8ghz"] = lba_8ghz_rms

    # Adding the overall 2020 atca fluxes 
    try:
        t = Table.read(f"{data_dir}{tar}_atca.fits")
        atca_flx = np.squeeze(np.array(t["ATCA"]))
        atca_err = np.array(np.sqrt((0.05*atca_flx)**2 + (0.0004**2)))
    except FileNotFoundError:
        print(f"No catalogue for {tar} atca")
        atca_flx = [np.nan * 17]
        atca_err = atca_flx
    atca_src = np.vstack((atca_flx,atca_err))
    src["flx_atca_20"] = atca_flx
    src["err_atca_20"] = atca_err

    src_dict[lba_targets[i]] = src

for i in range(len(lba_targets)):
    plotting_functions.plt_sed("/data/LBA/plots/",src_dict[lba_targets[i]],papersize=False)