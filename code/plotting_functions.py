from os import minor
from typing import Type
import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
import CFigTools.CustomFigure as CF
import pandas as pd 

def plt_sed(
    save_dir,
    src_dict,
    papersize=True,
    colors=cmr.take_cmap_colors(
        "cmr.gothic", 8, cmap_range=(0.15, 0.8), return_fmt="hex"
    ),
    epochnms=[
        "2013",
        "2014",
        "Jan20",
        "Mar20",
        "Apr20",
        "May20",
        "Jul20",
        "Sept20",
    ],
    frequency=np.array(
        [
            0.076,
            0.084,
            0.092,
            0.099,
            0.107,
            0.115,
            0.122,
            0.130,
            0.143,
            0.151,
            0.158,
            0.166,
            0.174,
            0.181,
            0.189,
            0.197,
            0.204,
            0.212,
            0.220,
            0.227,
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
            0.150,
            0.408,
            0.843,
            1.400,
            0.074,
            20,
            8.6,
            4.8,
            0.8875,
            2.4,
            8.3,
        ]
    ),
):
    # note of which frequencies are what: 
    # mwa = frequency[0:20]
    # atca = frequency[20:37]
    # xtra = frequency[37:46]
    # lba = [46:48]

    target = src_dict["MWA Name"].strip("GLEAM ")[0:7]
    lba_colors = cmr.take_cmap_colors("cmr.flamingo", 3, cmap_range=(0.15, 0.8), return_fmt="hex")
    if papersize is True:
        figsize=(3.5, 2.25)
        fontsize=8
        titlefontsize=12
        xlabelfontsize=8
        ylabelfontsize=8
        ticklabelfontsize=8
        majorticklength=2.5 
        minorticklength=2
        tickwidth=0.5
        ext = "pdf"
        print("Paper size")
        s=5
        elinewidth=0.5
        capsize=1
    else: 
        figsize=(15,10)
        fontsize=20
        titlefontsize=40
        xlabelfontsize=25
        ylabelfontsize=25
        ticklabelfontsize=20
        majorticklength=7.5 
        minorticklength=6
        tickwidth=1
        ext = "png"
        s=50
        elinewidth = 1
        capsize=1.5
    f = CF.sed_fig(figsize=figsize)
    # Plotting the extra surveys first 
    extra_markers = ["X", "s", "*", "p", "D", ">", "<", "^", "1"]
    extra_frequencies = frequency[37:46]
    for i in range(len(src_dict["flx_supp"])):
        xtra = src_dict["flx_supp"]
        if xtra[i] == 0.0:
            xtra[i] = np.nan
        f.plot_point(
            extra_frequencies[i],
            xtra[i],
            marker=extra_markers[i],
            marker_color="k",
            s=s,
            alpha=0.5,
        )
    
    # Plotting the mwa points 
    f.plot_spectrum(
        frequency[0:20],
        src_dict["flx_mwa_13"],
        src_dict["err_mwa_13"],
        marker="o",
        label="2013",
        marker_color=colors[0],
        s=s,
        elinewidth=elinewidth,
        capsize=capsize,
    )
    f.plot_spectrum(
        frequency[0:20],
        src_dict["flx_mwa_14"],
        src_dict["err_mwa_14"],
        marker="o",
        label="2014",
        marker_color=colors[1],
        s=s,
        elinewidth=elinewidth,
        capsize=capsize,
    )

    # Plotting ATCA spectrum 
    try:
        f.plot_spectrum(
            frequency[20:37],
            src_dict["flx_atca_20"],
            src_dict["err_atca_20"],
            marker="o",
            label="2020",
            marker_color="k",
            s=s,
            alpha=0.5,
            elinewidth=elinewidth,
            capsize=capsize,
        )
    except TypeError or ValueError:
        print("Not plotting ATCA fluxes, check that you've read in the correct values")
    # Plotting the LBA points 
    lba_2ghz = src_dict["flx_lba_2ghz"]
    err_lba_2ghz = src_dict["err_lba_2ghz"]
    try:
        for i in range(len(lba_2ghz)):
            f.plot_point(
                2.4,
                lba_2ghz[i],
                # err_lba_2ghz[i],
                marker="o",
                marker_color=lba_colors[i],
                s=s,
                # alpha=0.5,
                elinewidth=elinewidth,
                capsize=capsize,
            )
        lba_8ghz = src_dict["flx_lba_8ghz"]
        err_lba_8ghz = src_dict["err_lba_8ghz"]
        for i in range(len(lba_2ghz)):
            f.plot_point(
                8.3,
                lba_8ghz[i],
                # err_lba_8ghz[i],
                marker="o",
                marker_color=lba_colors[i],
                s=s,
                # alpha=0.5,
                elinewidth=elinewidth,
                capsize=capsize,
            )
    except TypeError or ValueError:
        print("Not plotting LBA, check there is a catalogue to plot")
        
    f.legend(loc="lower center", fontsize=fontsize)
    f.format(xunit="GHz",xlabelfontsize=xlabelfontsize, ylabelfontsize=ylabelfontsize, ticklabelfontsize=ticklabelfontsize,majorticklength=majorticklength, minorticklength=minorticklength, tickwidth=tickwidth)
    f.title(src_dict["MWA Name"], fontsize=titlefontsize)
    f.save(f"{save_dir}{target}_sed", ext=ext)
    plt.close()
    plt.clf()
    return


def make_nm_arrays():
    # Setting arrays for masks and extracting sources
    exts = [
        "107",
        "115",
        "123",
        "130",
        "143",
        "150",
        "158",
        "166",
        "174",
        "181",
        "189",
        "197",
        "204",
        "212",
        "220",
        "227",
    ]
    mwa_mask_names = ["Name", "S_076", "S_084", "S_092", "S_099","S_076_err", "S_084_err", "S_092_err", "S_099_err"]
    mwa_2013_fluxes = ["S_076", "S_084", "S_092", "S_099"]
    mwa_2014_fluxes = []
    mwa_2013_errors = ["S_076_err", "S_084_err", "S_092_err", "S_099_err"]
    mwa_2014_errors = []
    for ext in exts:
        mwa_mask_names.extend((f"S_{ext}_yr1", f"S_{ext}_yr2",f"local_rms_{ext}_yr1", f"local_rms_{ext}_yr2"))
        mwa_2013_fluxes.append(f"S_{ext}_yr1")
        mwa_2014_fluxes.append(f"S_{ext}_yr2")
        mwa_2013_errors.append(f"local_rms_{ext}_yr1")
        mwa_2014_errors.append(f"local_rms_{ext}_yr2")
    # Extras
    xtra_fluxes = [
        "S_tgss",
        "S_mrc",
        "S_sumss",
        "S_nvss",
        "S_vlssr",
        "S20",
        "S8",
        "S5",
        "total_flux_source",  # RACS
    ]
    mwa_mask_names = np.hstack((mwa_mask_names,xtra_fluxes))
    mwa_2013 = np.vstack((mwa_2013_fluxes,mwa_2013_errors))
    mwa_2014 = np.vstack((mwa_2014_fluxes,mwa_2014_errors))
    return mwa_mask_names, mwa_2013, mwa_2014, xtra_fluxes

def read_fluxes(mask_nms, master_data_dir="/data/raw_data/", mwa_targets = ["GLEAM J022744-062106", "GLEAM J032237-482010", "GLEAM J223933-451414"]):
    master_pop_pd = pd.read_csv(f"{master_data_dir}/master_pop_extended.csv")
    lba_mask = master_pop_pd["Name"].isin(mwa_targets)
    lba_master_pop = master_pop_pd[lba_mask]
    lba_pop = lba_master_pop[mask_nms]
    return lba_pop