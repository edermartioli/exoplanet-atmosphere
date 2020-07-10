# -*- coding: iso-8859-1 -*-
"""
    Created on Jul 06 2020
    
    Description: python library for handling SPIRou data
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import numpy as np
import matplotlib.pyplot as plt

import astropy.io.fits as fits
from astropy.io import ascii
from scipy import constants

import os,sys

import time

#--- Load a spirou spectrum from e.fits or t.fits file (which are the default products at CADC)
# This function preserves the spectral order structure
def load_spirou_AB_efits_spectrum(input, nan_pos_filter=True, preprocess=False, apply_BERV_to_preprocess=True, source_rv=0., normalization_in_preprocess=True, normalize_blaze=True) :
    
    # open fits file
    hdu = fits.open(input)
    
    if input.endswith("e.fits") :
        WaveAB = hdu["WaveAB"].data
        FluxAB = hdu["FluxAB"].data
        #BlazeAB = hdu[9].data / np.median(hdu[9].data)
        if normalize_blaze :
            BlazeAB = hdu["BlazeAB"].data / np.nanmean(hdu["BlazeAB"].data)
        else :
            BlazeAB = hdu["BlazeAB"].data
        
        WaveC = hdu["WaveC"].data
        FluxC = hdu["FluxC"].data
        #BlazeC = hdu["BlazeC"].data / np.median(hdu["BlazeC"].data)
        BlazeC = hdu["BlazeC"].data / np.nanmean(hdu["BlazeC"].data)

    elif input.endswith("t.fits") :
        WaveAB = hdu["WaveAB"].data
        FluxAB = hdu["FluxAB"].data
        #BlazeAB = hdu[3].data / np.median(hdu[3].data)
        if normalize_blaze :
            BlazeAB = hdu["BlazeAB"].data / np.nanmean(hdu["BlazeAB"].data)
        else :
            BlazeAB = hdu["BlazeAB"].data
        Recon = hdu["Recon"].data
    else :
        print("ERROR: input file type not recognized")
        exit()

    WaveABout, FluxABout, BlazeABout = [], [], []
    WaveCout, FluxCout, BlazeCout = [], [], []
    Reconout = []
    for i in range(len(WaveAB)) :
        if nan_pos_filter :
            # mask NaN values
            nanmask = np.where(~np.isnan(FluxAB[i]))
            # mask negative and zero values
            negmask = np.where(FluxAB[i][nanmask] > 0)

            if len(WaveAB[i][nanmask][negmask]) :
                WaveABout.append(WaveAB[i][nanmask][negmask])
                FluxABout.append(FluxAB[i][nanmask][negmask])
                BlazeABout.append(BlazeAB[i][nanmask][negmask])
                if input.endswith("e.fits") :
                    WaveCout.append(WaveC[i][nanmask][negmask])
                    FluxCout.append(FluxC[i][nanmask][negmask])
                    BlazeCout.append(BlazeC[i][nanmask][negmask])
                elif input.endswith("t.fits") :
                    Reconout.append(Recon[i][nanmask][negmask])
            else :
                WaveABout.append(np.array([]))
                FluxABout.append(np.array([]))
                BlazeABout.append(np.array([]))
                if input.endswith("e.fits") :
                    WaveCout.append(np.array([]))
                    FluxCout.append(np.array([]))
                    BlazeCout.append(np.array([]))
                elif input.endswith("t.fits") :
                    Reconout.append(np.array([]))

        else :
            WaveABout.append(WaveAB[i])
            FluxABout.append(FluxAB[i])
            BlazeABout.append(BlazeAB[i])
            if input.endswith("e.fits") :
                WaveCout.append(WaveC[i])
                FluxCout.append(FluxC[i])
                BlazeCout.append(BlazeC[i])
            elif input.endswith("t.fits") :
                Reconout.append(Recon[i])

    loc = {}
    loc['filename'] = input
    loc['header0'] = hdu[0].header
    loc['header1'] = hdu[1].header

    loc['WaveAB'] = WaveABout
    loc['FluxAB'] = FluxABout
    loc['BlazeAB'] = BlazeABout
    
    if input.endswith("e.fits") :
        loc['WaveC'] = WaveCout
        loc['FluxC'] = FluxCout
        loc['BlazeC'] = BlazeCout
        loc['headerC'] = hdu['FluxC'].header

    elif input.endswith("t.fits") :
        loc['Recon'] = Reconout

    if preprocess :
        # Pre-process spectrum to normalize data, remove nans and zeros, apply BERV correction if requested, etc
        loc = pre_process(loc, apply_BERV=apply_BERV_to_preprocess, source_rv=source_rv, normalize=normalization_in_preprocess, nan_pos_filter=nan_pos_filter)

    return loc



def pre_process(spectrum, apply_BERV=True, source_rv=0., normalize=True, nan_pos_filter=True) :
    out_wl, out_flux, out_fluxerr = [], [], []
    
    if normalize :
        norm_chunk = get_chunk_data(spectrum, 1054., 1058., 6, rv_overscan = 0., source_rv=0.0, apply_BERV=False, cal_fiber=False, normalize=False, nan_pos_filter=True, plot=False)
        spectrum["normalization_factor"] = np.median(norm_chunk['flux'])
    
    out_continuum = []
    
    for order in range(len(spectrum['WaveAB'])) :
        if len(spectrum['WaveAB'][order]) :
            wl0, wlf = spectrum['WaveAB'][order][0], spectrum['WaveAB'][order][-1]
            loc = get_chunk_data(spectrum, wl0, wlf, order, rv_overscan = 0., source_rv=source_rv, cal_fiber=False, apply_BERV=apply_BERV, normalize=normalize, nan_pos_filter=nan_pos_filter, plot=False)
            out_wl.append(loc['wl'])
            out_flux.append(loc['flux'])
            out_fluxerr.append(loc['fluxerr'])
        else :
            out_wl.append(np.array([]))
            out_flux.append(np.array([]))
            out_fluxerr.append(np.array([]))

    spectrum['wl'] = out_wl
    spectrum['flux'] = out_flux
    spectrum['fluxerr'] = out_fluxerr

    return spectrum

#### Function to get chunk data #########
def get_chunk_data(spectrum, wl0, wlf, order, rv_overscan = 100.0, source_rv=0.0, apply_BERV=True, cal_fiber=False, normalize=True, nan_pos_filter=True, plot=False) :

    loc = {}
    loc['order'] = order

    loc['filename'] = spectrum['filename']
    loc['wl0'] = wl0
    loc['wlf'] = wlf

    wlc = (wl0 + wlf) / 2.
    
    # set overscan to avoid edge issues
     # in km/s
    loc['rv_overscan'] = rv_overscan

    #rv_overscan = 0.0 # in km/s
    dwl_1 = rv_overscan * wl0 / (constants.c / 1000.)
    dwl_2 = rv_overscan * wlf / (constants.c / 1000.)
    
    # get BERV from header
    if apply_BERV :
        BERV = spectrum['header1']['BERV']
    else :
        BERV = 0.

    # get DETECTOR GAIN and READ NOISE from header
    gain, rdnoise = spectrum['header0']['GAIN'], spectrum['header0']['RDNOISE']

    if cal_fiber :
        if nan_pos_filter :
            # mask NaN values
            nanmask = np.where(~np.isnan(spectrum['FluxAB'][order]))
            # mask negative and zero values
            negmask = np.where(spectrum['FluxAB'][order][nanmask] > 0)
            # set calibration fiber flux and wavelength vectors
            flux = spectrum['FluxC'][order][nanmask][negmask] / spectrum['BlazeC'][order][nanmask][negmask]
            wave = spectrum['WaveC'][order][nanmask][negmask]
            # calculate flux variance
            fluxerr = np.sqrt((spectrum['FluxAB'][order][nanmask][negmask] + (rdnoise * rdnoise / gain * gain) ) / spectrum['BlazeAB'][order][nanmask][negmask])
        else :
            # set calibration fiber flux and wavelength vectors
            flux = spectrum['FluxC'][order] / spectrum['BlazeC'][order]
            wave = spectrum['WaveC'][order]
            # calculate flux variance
            fluxerr = np.sqrt((spectrum['FluxAB'][order] + (rdnoise * rdnoise / gain * gain) ) / spectrum['BlazeAB'][order])

    else :
        if nan_pos_filter :
            # mask NaN values
            nanmask = np.where(~np.isnan(spectrum['FluxAB'][order]))
            # mask negative and zero values
            negmask = np.where(spectrum['FluxAB'][order][nanmask] > 0)
            # set science fiber flux and wavelength vectors
            flux = spectrum['FluxAB'][order][nanmask][negmask] / spectrum['BlazeAB'][order][nanmask][negmask]
            # apply BERV correction - Barycentric Earth Radial Velocity (BERV)
            wave = spectrum['WaveAB'][order][nanmask][negmask] * ( 1.0 + (BERV - source_rv) / (constants.c / 1000.) )
            # calculate flux variance
            fluxerr = np.sqrt(spectrum['FluxAB'][order][nanmask][negmask] + (rdnoise * rdnoise / gain * gain)) / spectrum['BlazeAB'][order][nanmask][negmask]
            if 'Recon' in spectrum.keys():
                recon = spectrum['Recon'][order][nanmask][negmask]
        else :
            # set science fiber flux and wavelength vectors
            flux = spectrum['FluxAB'][order] / spectrum['BlazeAB'][order]
            # apply BERV correction - Barycentric Earth Radial Velocity (BERV)
            wave = spectrum['WaveAB'][order] * ( 1.0 + (BERV - source_rv) / (constants.c / 1000.) )
            # calculate flux variance
            fluxerr = np.sqrt(spectrum['FluxAB'][order] + (rdnoise * rdnoise / gain * gain)) / spectrum['BlazeAB'][order]
            # get telluric absorption spectrum
            if 'Recon' in spectrum.keys():
                recon = spectrum['Recon'][order]
    
    # set wavelength masks
    wlmask = np.where(np.logical_and(wave > wl0 - dwl_1, wave < wlf + dwl_2))

    if len(flux[wlmask]) == 0 :
        loc['wl'] = np.array([])
        loc['flux'] = np.array([])
        loc['fluxerr'] = np.array([])
        if 'Recon' in spectrum.keys():
            loc['recon'] = np.array([])
        return loc

    # measure continuum and normalize flux if nomalize=True
    if normalize == True:
        # Calculate normalization factor
        normalization_factor = spectrum['normalization_factor']

        loc['normalization_factor'] = normalization_factor

        # normalize flux
        flux = flux / normalization_factor
        fluxerr = fluxerr / normalization_factor

    # mask data
    flux, fluxerr, wave = flux[wlmask], fluxerr[wlmask], wave[wlmask]
    if 'Recon' in spectrum.keys():
        recon = recon[wlmask]

    if plot :
        plt.plot(wave,flux)
        plt.errorbar(wave,flux,yerr=fluxerr,fmt='.')
        if 'Recon' in spectrum.keys():
            plt.plot(wave,flux*recon,'-',linewidth=0.3)
        plt.show()

    loc['order'] = order
    loc['wl'] = wave
    if cal_fiber :
        loc['flux'] = flux / np.max(flux)
        loc['fluxerr'] = fluxerr / np.max(flux)
    else :
        loc['flux'] = flux
        loc['fluxerr'] = fluxerr
        if 'Recon' in spectrum.keys():
            loc['recon'] = recon

    return loc
##-- end of function

def get_normalization_factor(flux_order, i_max, norm_window=50) :
    
    min_i = i_max - norm_window
    max_i = i_max + norm_window
    if min_i < 0 :
        min_i = 0
        max_i = min_i + 2 * norm_window
    if max_i >= len(flux_order) :
        max_i = len(flux_order) - 1

    # Calculate normalization factor as the median of flux within window around maximum signal
    normalization_factor = np.nanmedian(flux_order[min_i:max_i])

    return normalization_factor


# Functon below defines SPIRou spectral orders, useful wavelength ranges and the NIR bands
def spirou_order_mask():
    order_mask = [[0, 967, 980, 'Y'],
             [1, 977, 994, 'Y'],
             [2, 989, 1005.5,'Y'],
             [3, 1001., 1020,'Y'],
             [4, 1018, 1035,'Y'],
             [5, 1028, 1050,'Y'],
             [6, 1042.3, 1065,'Y'],
             [7, 1055, 1078,'Y'],
             [8, 1071.5, 1096,'Y'],
             [9, 1084.5, 1107,'Y'],
             [10, 1107, 1123,'J'],
             [11, 1123, 1140.2,'J'],
             [12, 1137.25, 1162,'J'],
             [13, 1150, 1178,'J'],
             [14, 1168, 1198,'J'],
             [15, 1186, 1216,'J'],
             [16, 1204, 1235,'J'],
             [17, 1223, 1255,'J'],
             [18, 1243, 1275,'J'],
             [19, 1263, 1297,'J'],
             [20, 1284, 1321,'J'],
             [21, 1306, 1344,'J'],
             [22, 1327.5, 1367.5,'J'],
             [23, 1350.1, 1392,'J'],
             [24, 1374.3, 1415,'J'],
             [25, 1399.7, 1443.6,'H'],
             [26, 1426, 1470.9,'H'],
             [27, 1453.5, 1499,'H'],
             [28, 1482, 1528.6,'H'],
             [29, 1512, 1557.5,'H'],
             [30, 1544, 1591.1,'H'],
             [31, 1576.6, 1623,'H'],
             [32, 1608.5, 1658.9,'H'],
             [33, 1643.5, 1695,'H'],
             [34, 1679.8, 1733,'H'],
             [35, 1718, 1772,'H'],
             [36, 1758, 1813.5,'H'],
             [37, 1800.7, 1856.5,'H'],
             [38, 1843.9, 1902, 'H'],
             [39, 1890, 1949.5,'H'],
             [40, 1938.4, 1999.5, 'H'],
             [41, 1989.5, 2052, 'K'],
             [42, 2043, 2108, 'K'],
             [43, 2100, 2166, 'K'],
             [44, 2160, 2228,'K'],
             [45, 2223.5, 2293.6,'K'],
             [46, 2291, 2363,'K'],
             [47, 2362, 2437.2,'K'],
             [48, 2439, 2510,'K']]
             
    outorders, wl0, wlf, colors = [], [], [], []
    for order in order_mask:
        outorders.append(order[0])
        wl0.append(order[1])
        wlf.append(order[2])
        colors.append(order[3])
    
    loc = {}
    loc['orders'] = outorders
    loc['wl0'] = wl0
    loc['wlf'] = wlf
    loc['colors'] = colors
    return loc

def read_rvdata(rvfile) :
    
    rvdata = ascii.read(rvfile,data_start=2)
    rvbjds = np.array(rvdata['rjd']) + 2400000.
    rvs, rverrs = np.array(rvdata["vrad"]), np.array(rvdata["svrad"])
    
    return rvbjds, rvs, rverrs


def load_array_of_spirou_spectra(inputedata, rvfile="", verbose=False, plot=False) :

    loc = {}
    loc["input"] = inputedata

    # Read input RVs from file

    if rvfile == "":
        if verbose:
            print("Reading RV file...")
        loc["rvfile"] = rvfile
        rvbjds, rvs, rverrs = np.zeros(len(inputedata)), np.zeros(len(inputedata)), np.zeros(len(inputedata))
    else :
        rvbjds, rvs, rverrs = read_rvdata(rvfile)

        if len(rvs) != len(inputedata):
            print("WARNING: size of RVs is different than number of input *e.fits files")
            print("*** Ignoring input RVs ***")
            rvbjds, rvs, rverrs = np.zeros(len(inputedata)), np.zeros(len(inputedata)), np.zeros(len(inputedata))
    #---

    orders = spirou_order_mask()
    loc["orders"] = orders

    output = []
    spectra = []

    for i in range(len(inputedata)) :
        # create output CCF file name based on input spectrum file name
        outfilename = inputedata[i].split(".")[-2] + "_s.fits"
        output.append(outfilename)

        if verbose :
            print("Loadinng spectrum:",inputedata[i],"{0}/{1}".format(i,len(inputedata)-1))

        # Load SPIRou *e.fits spectrum file
        spectrum = load_spirou_AB_efits_spectrum(inputedata[i], nan_pos_filter=True, preprocess=True, apply_BERV_to_preprocess=False, normalization_in_preprocess=True)

        # set source RVs
        spectrum["source_rv"] = rvs[i]
        spectrum["rvfile"] = rvfile
        spectrum['RV'] = rvs[i]
        spectrum['RVERR'] = rverrs[i]

        spectrum['BJD_mid'] = spectrum['header1']["BJD"] + (spectrum['header0']['MJDEND'] - spectrum['header0']['MJD-OBS']) / 2.
        spectrum['BERV'] = float(spectrum['header1']['BERV'])
        spectrum['snr32'] = spectrum['header1']['SNR32']
        spectrum['airmass'] = spectrum['header0']['AIRMASS']

        wl0 = 0
        out_wl, out_flux, out_fluxerr, out_order = np.array([]), np.array([]), np.array([]), np.array([])
        NORDERS = 0
        for order in range(len(orders['orders'])) :
            mask = ~np.isnan(spectrum['flux'][order])
            mask &= spectrum['wl'][order] > orders['wl0'][order]
            mask &= spectrum['wl'][order] < orders['wlf'][order]

            if len(spectrum['wl'][order][mask]) :
                wl, flux = spectrum['wl'][order][mask], spectrum['flux'][order][mask]
                fluxerr = spectrum['fluxerr'][order][mask]
    
                wl0 = wl[-1]

                if plot :
                    p = plt.plot(wl, flux)
                    color = p[0].get_color()
                    plt.plot(spectrum['wl'][order], spectrum['flux'][order], color=color, lw=0.3, alpha=0.6)

                out_wl = np.append(out_wl,wl)
                out_flux = np.append(out_flux,flux)
                out_fluxerr = np.append(out_fluxerr,fluxerr)
                order_vec = np.full_like(wl,float(NORDERS))
                out_order = np.append(out_order,order_vec)
                NORDERS += 1

        if plot :
            plt.xlabel(r"wavelength [nm]")
            plt.xlabel(r"flux")
            plt.show()
            exit()

        spectrum['NORDERS'] = NORDERS

        spectrum['out_wl'] = out_wl
        spectrum['out_flux'] = out_flux
        spectrum['out_fluxerr'] = out_fluxerr
        spectrum['out_order'] = out_order

        spectra.append(spectrum)

    loc["output"] = output
    loc["spectra"] = spectra

    return loc


def save_spectrum_to_fits(spectrum, output, wl0=0., wlf=1e50) :

    header = fits.Header()
    
    header.set('ORIGIN', "spiroulib.save_spectrum_to_fits()")
    header.set('FILENAME', output)
    header.set('UTCSAVED', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
    
    header.set('BJDMID', spectrum['BJD_mid'], "Barycentric Julian Date at middle of exposure")
    header.set('BERV', spectrum['BERV'], "Barycentric Earth Radial Velocity correction (km/s)")
    header.set('SNR32', spectrum['snr32'], "Signal-to-Noise ratio for order 32")
    header.set('AIRMASS', spectrum['airmass'], "Airmass")
    header.set('RV', spectrum['RV'], "Radial velocity of source (km/s)")
    header.set('RVERR', spectrum['RVERR'], "Radial velocity of source (km/s)")
    header.set('RVFILE', spectrum['rvfile'], "Radial velocity data file")
    
    header.set('NORDERS', spectrum['NORDERS'], "Number of spectral orders")
    
    if wl0 > 0 :
        header.set('WL0', wl0, "Initial wavelength [nm]")
    if wlf < 1e50 :
        header.set('WLF', wlf, "Final wavelength [nm]")
    
    header.set('TTYPE1', "WAVE")
    header.set('TUNIT1', "NM")
    
    header.set('TTYPE2', "FLUX")
    header.set('TUNIT2', "COUNTS")
    
    header.set('TTYPE3', "FLUXERR")
    header.set('TUNIT3', "COUNTS")
    
    header.set('TTYPE4', "ORDER")
    header.set('TUNIT4', "NUMBER")

    wlmask = ((spectrum['out_wl'] > wl0) & (spectrum['out_wl'] < wlf))

    #minorder, maxorder = np.min(spectrum['out_order'][wlmask]), np.max(spectrum['out_order'][wlmask])
    #wlmask &= spectrum['out_order'] > minorder
    #wlmask &= spectrum['out_order'] < maxorder

    outhdulist = []
    
    primary_hdu = fits.PrimaryHDU(header=header)
    outhdulist.append(primary_hdu)
    
    hdu_wl = fits.ImageHDU(data=spectrum['out_wl'][wlmask], name="WAVE")
    outhdulist.append(hdu_wl)
    
    hdu_flux = fits.ImageHDU(data=spectrum['out_flux'][wlmask], name="FLUX")
    outhdulist.append(hdu_flux)
    
    hdu_fluxerr = fits.ImageHDU(data=spectrum['out_fluxerr'][wlmask], name="FLUXERR")
    outhdulist.append(hdu_fluxerr)
    
    hdu_order = fits.ImageHDU(data=spectrum['out_order'][wlmask], name="ORDER")
    outhdulist.append(hdu_order)
    
    mef_hdu = fits.HDUList(outhdulist)
    mef_hdu.writeto(output, overwrite=True)


def save_spirou_spectra_to_fits(dataset, wl0=0., wlf=1e50, overwrite=False, verbose=False) :

    for i in range(len(dataset["input"])) :
        
        output = dataset["output"][i]
        
        if os.path.exists(output) and not overwrite :
            print("File",output," already exists, skipping ...")
            continue
        
        if verbose :
            print("Saving spectrum:",dataset["input"][i]," to output=", output, "{0}/{1}".format(i,len(dataset["input"])-1))

        spectrum = dataset["spectra"][i]

        save_spectrum_to_fits(spectrum, output, wl0=wl0, wlf=wlf)
