from os import minor
from typing import Type
import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
import CFigTools.CustomFigure as CF
import pandas as pd 
from astropy.table import Table
import json
import gpscssmodels
from scipy.optimize import curve_fit

channel = ("69", "93", "121", "145", "169")
subchans_dict = {
    "69": ["072-080", "080-088", "088-095", "095-103"],
    "93": ["103-111", "111-118", "118-126", "126-134"],
    "121": ["139-147", "147-154", "154-162", "162-170"],
    "145": ["170-177", "177-185", "185-193", "193-200"],
    "169": ["200-208", "208-216", "216-223", "223-231"],
}
epochs = ["2020-04", "2020-05", "2020-07", "2020-10"]


def plt_sed(
    save_dir,
    src_dict,
    model_params,
    atca_model = gpscssmodels.powlawbreak,
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
            2.276875,
            8.416875,
        ]
    ),
):
    # note of which frequencies are what: 
    # mwa = frequency[0:20]
    # atca = frequency[20:37]
    # xtra = frequency[37:46]
    # lba = [46:48]

    target = src_dict["MWA Name"].strip("GLEAM ")[0:7]
    lba_colors = cmr.take_cmap_colors("cmr.freeze", 3, cmap_range=(0.15, 0.8), return_fmt="hex")
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

    # Plotting MWA 2020 Monitoring 
    for i in range(len(epochs)):
        f.plot_spectrum(
            frequency[0:20],
            src_dict[f"flx_mwa_{epochs[i]}"],
            src_dict[f"err_mwa_{epochs[i]}"],
            marker="o",
            label=epochs[i],
            marker_color=colors[i+4],
            s=s,
            alpha=1,
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
            label="ATCA 2020",
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
    lba_2ghz_sum = np.sum(lba_2ghz)
    err_lba_2ghz = src_dict["err_lba_2ghz"]
    lba_8ghz = src_dict["flx_lba_8ghz"]
    lba_8ghz_sum = np.sum(lba_8ghz)
    err_lba_8ghz = src_dict["err_lba_8ghz"]
    try:
        for i in range(len(lba_2ghz)):
            f.display_model(
                np.linspace(0.01, 25, num=10000),
                gpscssmodels.powlaw(np.linspace(0.01, 25, num=10000), *model_params[i]),
                color=lba_colors[i]
            )
    except TypeError: 
        f.display_model(
            np.linspace(0.01, 25, num=10000),
            gpscssmodels.powlaw(np.linspace(0.01, 25, num=10000), *model_params[0]),
            color=lba_colors[-1]
        )
    f.display_model(
        np.linspace(0.01, 25, num=10000),
        gpscssmodels.powlaw(np.linspace(0.01, 25, num=10000), *model_params[-2]),
        color='k',
        alpha=0.8,
    )
    f.display_model(
        np.linspace(0.01, 25, num=10000),
        atca_model(np.linspace(0.01, 25, num=10000), *model_params[-1]),
        color='k',
        alpha=0.5,
    )
    try:
        for i in range(len(lba_2ghz)):
            f.plot_point(
                frequency[-2],
                lba_2ghz[i],
                # err_lba_2ghz[i],
                marker="o",
                marker_color=lba_colors[i],
                s=s,
                # alpha=0.5,
                elinewidth=elinewidth,
                capsize=capsize,
            )
        for i in range(len(lba_2ghz)):
            f.plot_point(
                frequency[-1],
                lba_8ghz[i],
                # err_lba_8ghz[i],
                marker="o",
                marker_color=lba_colors[i],
                s=s,
                # alpha=0.5,
                elinewidth=elinewidth,
                capsize=capsize,
            )
    except ValueError:
        print("Not plotting LBA, check there is a catalogue to plot")
    except TypeError:
        f.plot_point(
            frequency[-1],
            lba_2ghz,
            # err_lba_2ghz[i],
            marker="o",
            marker_color=lba_colors[-1],
            s=s,
            label="LBA",
            # alpha=0.5,
            elinewidth=elinewidth,
            capsize=capsize,
        )
        f.plot_point(
            frequency[-2],
            lba_8ghz,
            # err_lba_8ghz[i],
            marker="o",
            marker_color=lba_colors[-1],
            s=s,
            # alpha=0.5,
            elinewidth=elinewidth,
            capsize=capsize,
        )
    f.plot_point(
        frequency[-2],
        lba_2ghz_sum,
        # err_lba_2ghz[i],
        marker="X",
        marker_color='k',
        s=s,
        label="integrated LBA",
        # alpha=0.5,
        elinewidth=elinewidth,
        capsize=capsize,
    )
    f.plot_point(
        frequency[-1],
        lba_8ghz_sum,
        # err_lba_8ghz[i],
        marker="X",
        marker_color='k',
        s=s,
        # alpha=0.5,
        elinewidth=elinewidth,
        capsize=capsize,
    )

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

def make_src_dict(lba_pop_pd, data_dir = "/data/LBA/catalogues/", lba_targets = ["j0227-0621", "j0322-482", "j2239-451"], mwa_targets = ["GLEAM J022744-062106", "GLEAM J032237-482010", "GLEAM J223933-451414"]):
    src_dict = {}
    for i in range(len(mwa_targets)):
        mask_nms, mwa_13, mwa_14, xtra_nms = make_nm_arrays()
        
        # MWA: Adding the 13 and 14 fluxes 
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

        # MWA: Adding 2020 monitoring fluxes
        mwa_tar = mwa_targets[i].strip("GLEAM ")[0:7]
        for j in range(len(epochs)):
            epoch = epochs[j]
            chan_flux = []
            err_chan_flux = []
            for c in range(len(channel)):
                subchans = subchans_dict[channel[c]]
                chan = channel[c]
                if epoch == "2020-04":
                    percentage = 0.05
                else:
                    percentage = 0.02
                for subchan in subchans:
                    try:
                        src_mwa_pd = pd.read_csv(
                            f"/data/MWA/{epoch}/{chan}/minimosaic/{mwa_tar}_{subchan}MHz_ddmod_scaled_comp_xmatch.csv"
                        )
                        mask = src_mwa_pd["Name"] == mwa_targets[i]
                        src_pd = src_mwa_pd[mask]
                        if src_pd.empty is False:
                            mwa_flux_chan = np.squeeze(src_pd["int_flux"].values)
                            mwa_errs_chan = np.squeeze(
                                np.sqrt(src_pd["local_rms"]) ** 2
                                + (percentage * mwa_flux_chan) ** 2
                            )
                            chan_flux.append(mwa_flux_chan)
                            err_chan_flux.append(mwa_errs_chan)
                        else:
                            chan_flux.append(np.nan)
                            err_chan_flux.append(np.nan)
                            pass
                    except (FileNotFoundError, KeyError):
                        chan_flux.append(np.nan)
                        err_chan_flux.append(np.nan)
                        pass

            src[f"flx_mwa_{epoch}"] = np.squeeze(chan_flux)
            src[f"err_mwa_{epoch}"] = np.squeeze(err_chan_flux)

        # SUPP: Adding the fluxes of supplementary surveys
        extra_fluxes = np.squeeze(lba_pop_pd.loc[lba_pop_pd["Name"]==mwa_targets[i]][xtra_nms].to_numpy())
        extra_fluxes[0] = extra_fluxes[0]*0.001
        extra_fluxes[2] = extra_fluxes[2]*0.001
        extra_fluxes[3] = extra_fluxes[3]*0.001
        extra_fluxes[5] = extra_fluxes[5]*0.001
        extra_fluxes[6] = extra_fluxes[6]*0.001
        extra_fluxes[7] = extra_fluxes[7]*0.001
        extra_fluxes[8] = extra_fluxes[8]*0.001
        src["flx_supp"] = extra_fluxes

        # LBA: Adding the LBA fluxes for each region 
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

        # ATCA: Adding the overall 2020 atca fluxes 
        try:
            t = Table.read(f"/data/LBA/catalogues/{tar}_atca.fits")
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
    return src_dict

def fit_lba(fluxes, fit_function= gpscssmodels.powlaw, frequency=[2.276875,8.416875]):
    if fit_function == gpscssmodels.powlawbreak:
        p0 = [1.6,0.03,7]
        params, _ = curve_fit(fit_function, frequency, fluxes, p0=p0)
    else: 
        params, _ = curve_fit(fit_function, frequency, fluxes)
    return params 