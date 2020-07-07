# -*- coding: iso-8859-1 -*-
"""
    Created on Jul 05 2020
    
    Description: python library for handling the library of exoplanet atmosphere models
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import os, sys

import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

import json

import scipy.signal as signal
import scipy.interpolate as sint
from scipy.interpolate import interp1d

from copy import copy, deepcopy

def get_spectrum_info_from_fits(filename) :
    
    """
        Get information in model FITS header
        
        :param filename: string, input exoplanet atmosphere model FITS file name

        :return loc: dict, dictionary containing information on input model
        """

    loc = {}
    hdr = fits.getheader(filename,0)
    # store filename and filepath
    loc['filename'] = os.path.basename(filename)
    loc['filepath'] = os.path.abspath(filename)
    
    for key in hdr.keys() :
        loc[key] = hdr[key]
    
    return loc


def get_best_model_file(modeldb, T_EQU=1100, R_PL=1.25, M_PL=1.81, species='all') :
    """
        Get file path of best model in input database
        
        :param modeldb: string, json database file name
        :param T_EQU: float, input equilibrium temperature [K]
        :param R_PL: float, input planet radius [RJup]
        :param M_PL: float, input planet mass [MJup]
        :param species: string, input molecular species

        :return filepath: string, file path for best model matching input params
        """

    # load json database file containig all models in the library
    try :
        with open(modeldb, 'r') as f:
            datastore = json.load(f)
    except :
        print("ERROR: could not open models database file ",modeldb)
        exit()

    mind_T_EQU = 1.e20
    mind_R_PL, mind_M_PL = 1.e20, 1.e20

    minkey = None

    # loop over all entried in the database to get best matched model file
    for key in datastore.keys() :
        d_T_EQU = np.abs(datastore[key]['T_EQU'] - T_EQU)
        d_R_PL = np.abs(datastore[key]['R_PL'] - R_PL)
        d_M_PL = np.abs(datastore[key]['M_PL'] - M_PL)
    
        if d_T_EQU <= mind_T_EQU and d_R_PL <= mind_R_PL and d_M_PL <= mind_M_PL and (species in datastore[key]['filename']):
            mind_T_EQU = d_T_EQU
            mind_R_PL, mind_M_PL = d_R_PL, d_M_PL
            minkey = key

    # return best model file path
    return datastore[minkey]['filepath']


def load_spectrum_from_fits(filename, wlmin=900., wlmax=2500.) :

    # initialize model spectrum with header info
    loc = get_spectrum_info_from_fits(filename)
    
    # load fits data
    hdu = fits.open(filename)
    hdr = hdu[0].header
    
    wl_0 = hdr["CRVAL1"]
    wl_f = hdr["CRVALEND"]
    wl_step = hdr["CDELT1"]
    
    wl_num = len(hdu["TRANSMISSION"].data)
    
    # create wavelength array
    wl = np.geomspace(wl_0*1000., wl_f*1000., wl_num)
    
    wlmask = np.where(np.logical_and(wl > wlmin, wl < wlmax))
    loc['wl'] = wl[wlmask]

    if "TRANSMISSION" in hdu :
        transmission = hdu["TRANSMISSION"].data
        #max_transm = np.max(transmission[wlmask])
        #loc['transmission'] = (max_transm - transmission[wlmask]) / np.max((max_transm - transmission[wlmask]))
        loc['transmission'] = transmission[wlmask]
        
    if "EMISSION" in hdu :
        emission = hdu["EMISSION"].data
        loc['emission'] = emission[wlmask]

    return loc


def species_colors() :
    colors = {}
    colors["h2o"]='#1f77b4'
    colors["co2"]='#ff7f0e'
    colors["ch4"]='#2ca02c'
    colors["o2"]='#d62728'
    colors["co"]='#9467bd'
    colors["k"]='#8c564b'
    colors["n2o"]='#e377c2'
    colors["no"]='#7f7f7f'
    colors["so2"]='#bcbd22'
    colors["no2"]='#17becf'
    colors["nh3"]='blue'
    colors["hno3"]='green'
    colors["o3"]='grey'

    return colors


"""
def calculate_model(modeldb, observed_wl, resolution=0., teff=5777., logg=4.4374, feh=0.0, normalize=False) :

    loc = {}
    loc["modeldb"] = modeldb
    loc["resolution"] = resolution
    loc["teff"] = teff
    loc["logg"] = logg
    loc["feh"] = feh

    loc["model_file"] = get_BTSettl_best_model_file(modeldb, teff=teff)
    
    model_spectrum = get_BTSettl_spectrum_from_fits(loc["model_file"])
    
    wl0 = observed_wl[0] * (1.0 - 2.0 / (constants.c / 1000.))
    wlf = observed_wl[-1] * (1.0 + 2.0 / (constants.c / 1000.))

    if wlf > model_spectrum['wl'][-1] :
        print("WARNING: Requested wavelengths limits are outside permitted boundaries for stellar model")
        print(wl0, ">", model_spectrum['wl'][0], "and", wlf, "<", model_spectrum['wl'][-1])
        print("Setting wlf=2480 nm")
        wlf = 2499.
        mask = observed_wl < wlf * (1.0 - 2.0 / (constants.c / 1000.))
        observed_wl = observed_wl[mask]
    
    wlmask = model_spectrum['wl'] > wl0
    wlmask &= model_spectrum['wl'] < wlf
    
    model_chunk = {}
    model_chunk['wl'] = model_spectrum['wl'][wlmask]
    model_chunk['flux'] = model_spectrum['flux'][wlmask]
    model_chunk['fluxerr'] = model_spectrum['fluxerr'][wlmask]

    resampled_model_chunk = {}
    resampled_model_chunk['wl'] = deepcopy(observed_wl)
    resampled_model_chunk['fluxerr'] = np.zeros_like(observed_wl)

    resampled_model_chunk['flux'] = espectrolib.interp_spectrum(observed_wl, model_chunk['wl'], model_chunk['flux'], kind='linear')

    if normalize :
        cont_model = fit_continuum(resampled_model_chunk['wl'][mask], resampled_model_chunk['flux'][mask], function='polynomial', order=3, nit=5, rej_low=2.0, rej_high=2.5, grow=1, med_filt=0, percentile_low=0., percentile_high=100.,min_points=10, xlabel="", ylabel="", plot_fit=True, verbose=False)
    
        resampled_model_chunk['flux'][mask] /= cont_model

    if resolution != 0. :
        resampled_model_chunk = convolve_spectrum(resampled_model_chunk, to_resolution=resolution)
    
    loc["wl"] = resampled_model_chunk["wl"]
    loc["flux"] = resampled_model_chunk["flux"]
    loc["fluxerr"] = resampled_model_chunk["fluxerr"]

    return loc
"""

def fit_continuum(wav, spec, function='polynomial', order=3, nit=5, rej_low=2.0,
    rej_high=2.5, grow=1, med_filt=0, percentile_low=0., percentile_high=100.,
                  min_points=10, xlabel="", ylabel="", plot_fit=True, verbose=False):
    """
    Continuum fitting re-implemented from IRAF's 'continuum' function
    in non-interactive mode only but with additional options.

    :Parameters:
    
    wav: array(float)
        abscissa values (wavelengths, velocities, ...)

    spec: array(float)
        spectrum values

    function: str
        function to fit to the continuum among 'polynomial', 'spline3'

    order: int
        fit function order:
        'polynomial': degree (not number of parameters as in IRAF)
        'spline3': number of knots

    nit: int
        number of iteractions of non-continuum points
        see also 'min_points' parameter

    rej_low: float
        rejection threshold in unit of residul standard deviation for point
        below the continuum

    rej_high: float
        same as rej_low for point above the continuum

    grow: int
        number of neighboring points to reject

    med_filt: int
        median filter the spectrum on 'med_filt' pixels prior to fit
        improvement over IRAF function
        'med_filt' must be an odd integer

    percentile_low: float
        reject point below below 'percentile_low' percentile prior to fit
        improvement over IRAF function
        "percentile_low' must be a float between 0. and 100.

    percentile_high: float
        same as percentile_low but reject points in percentile above
        'percentile_high'
        
    min_points: int
        stop rejection iterations when the number of points to fit is less than
        'min_points'

    plot_fit: bool
        if true display two plots:
            1. spectrum, fit function, rejected points
            2. residual, rejected points

    verbose: bool
        if true fit information is printed on STDOUT:
            * number of fit points
            * RMS residual
    """
    mspec = np.ma.masked_array(spec, mask=np.zeros_like(spec))
    # mask 1st and last point: avoid error when no point is masked
    # [not in IRAF]
    mspec.mask[0] = True
    mspec.mask[-1] = True
    
    mspec = np.ma.masked_where(np.isnan(spec), mspec)
    
    # apply median filtering prior to fit
    # [opt] [not in IRAF]
    if int(med_filt):
        fspec = signal.medfilt(spec, kernel_size=med_filt)
    else:
        fspec = spec
    # consider only a fraction of the points within percentile range
    # [opt] [not in IRAF]
    mspec = np.ma.masked_where(fspec < np.percentile(fspec, percentile_low),
        mspec)
    mspec = np.ma.masked_where(fspec > np.percentile(fspec, percentile_high),
        mspec)
    # perform 1st fit
    if function == 'polynomial':
        coeff = np.polyfit(wav[~mspec.mask], spec[~mspec.mask], order)
        cont = np.poly1d(coeff)(wav)
    elif function == 'spline3':
        knots = wav[0] + np.arange(order+1)[1:]*((wav[-1]-wav[0])/(order+1))
        spl = sint.splrep(wav[~mspec.mask], spec[~mspec.mask], k=3, t=knots)
        cont = sint.splev(wav, spl)
    else:
        raise(AttributeError)
    # iteration loop: reject outliers and fit again
    if nit > 0:
        for it in range(nit):
            res = fspec-cont
            sigm = np.std(res[~mspec.mask])
            # mask outliers
            mspec1 = np.ma.masked_where(res < -rej_low*sigm, mspec)
            mspec1 = np.ma.masked_where(res > rej_high*sigm, mspec1)
            # exlude neighbors cf IRAF's continuum parameter 'grow'
            if grow > 0:
                for sl in np.ma.clump_masked(mspec1):
                    for ii in range(sl.start-grow, sl.start):
                        if ii >= 0:
                            mspec1.mask[ii] = True
                    for ii in range(sl.stop+1, sl.stop+grow+1):
                        if ii < len(mspec1):
                            mspec1.mask[ii] = True
            # stop rejection process when min_points is reached
            # [opt] [not in IRAF]
            if np.ma.count(mspec1) < min_points:
                if verbose:
                    print("  min_points %d reached" % min_points)
                break
            mspec = mspec1
            if function == 'polynomial':
                coeff = np.polyfit(wav[~mspec.mask], spec[~mspec.mask], order)
                cont = np.poly1d(coeff)(wav)
            elif function == 'spline3':
                knots = wav[0] + np.arange(order+1)[1:]*((wav[-1]-wav[0])/(order+1))
                spl = sint.splrep(wav[~mspec.mask], spec[~mspec.mask], k=3, t=knots)
                cont = sint.splev(wav, spl)
            else:
                raise(AttributeError)
    # compute residual and rms
    res = fspec-cont
    sigm = np.std(res[~mspec.mask])
    if verbose:
        print("  nfit=%d/%d" %  (np.ma.count(mspec), len(mspec)))
        print("  fit rms=%.3e" %  sigm)
    # compute residual and rms between original spectrum and model
    # different from above when median filtering is applied
    ores = spec-cont
    osigm = np.std(ores[~mspec.mask])
    if int(med_filt) and verbose:
        print("  unfiltered rms=%.3e" %  osigm)
    # plot fit results
    if plot_fit:
        # overplot spectrum and model + mark rejected points
        fig1 = pl.figure(1)
        ax1 = fig1.add_subplot(111)
        ax1.plot(wav[~mspec.mask], spec[~mspec.mask],
            c='tab:blue', lw=1.0)
        # overplot median filtered spectrum
        if int(med_filt):
            ax1.plot(wav[~mspec.mask], fspec[~mspec.mask],
                c='tab:cyan', lw=1.0)
        ax1.scatter(wav[mspec.mask], spec[mspec.mask], s=20., marker='d',
        edgecolors='tab:gray', facecolors='none', lw=0.5)
        ax1.plot(wav, cont, ls='--', c='tab:orange')
        if nit > 0:
            # plot residuals and rejection thresholds
            fig2 = pl.figure(2)
            ax2 = fig2.add_subplot(111)
            ax2.axhline(0., ls='--', c='tab:orange', lw=1.)
            ax2.axhline(-rej_low*sigm, ls=':')
            ax2.axhline(rej_high*sigm, ls=':')
            ax2.scatter(wav[mspec.mask], res[mspec.mask],
                s=20., marker='d', edgecolors='tab:gray', facecolors='none',
                lw=0.5)
            ax2.scatter(wav[~mspec.mask], ores[~mspec.mask],
                marker='o', s=10., edgecolors='tab:blue', facecolors='none',
                lw=.5)
            # overplot median filtered spectrum
            if int(med_filt):
                ax2.scatter(wav[~mspec.mask], res[~mspec.mask],
                    marker='s', s=5., edgecolors='tab:cyan', facecolors='none',
                    lw=.2)
        if xlabel != "" :
            pl.xlabel(xlabel)
        if ylabel != "" :
            pl.ylabel(ylabel)
        pl.show()
    return cont


def __get_fwhm(lambda_peak, from_resolution, to_resolution):
    """
    Calculate the FWHM of the gaussian needed to convert
    a spectrum from one resolution to another at a given wavelength point.
    """
    if from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are equal or bigger than original")
    from_delta_lambda = (1.0*lambda_peak) / from_resolution
    to_delta_lambda = (1.0*lambda_peak) / to_resolution
    fwhm = np.sqrt(to_delta_lambda**2 - from_delta_lambda**2)
    return fwhm

def __fwhm_to_sigma(fwhm):
    """
    Calculate the sigma value from the FWHM.
    """
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    return sigma

def __convolve_spectrum(waveobs, flux, err, to_resolution, from_resolution=None):
    """
        Spectra resolution smoothness/degradation. Procedure:

        1) Define a bin per measure which marks the wavelength range that it covers.
        2) For each point, identify the window segment to convolve by using the bin widths and the FWHM.
        3) Build a gaussian using the sigma value and the wavelength values of the spectrum window.
        4) Convolve the spectrum window with the gaussian and save the convolved value.

        If "from_resolution" is not specified or its equal to "to_resolution", then the spectrum
        is convolved with the instrumental gaussian defined by "to_resolution".

        If "to_resolution" is specified, the convolution is made with the difference of
        both resolutions in order to degrade the spectrum.
    """
    if from_resolution is not None and from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are bigger than original")

    total_points = len(waveobs)
    convolved_flux = np.zeros(total_points)
    convolved_err = np.zeros(total_points)

    last_reported_progress = -1

    # Consider the wavelength of the measurements as the center of the bins
    # Calculate the wavelength distance between the center of each bin
    wave_distance = waveobs[1:] - waveobs[:-1]
    # Define the edge of each bin as half the wavelength distance to the bin next to it
    edges_tmp = waveobs[:-1] + 0.5 * (wave_distance)
    # Define the edges for the first and last measure which where out of the previous calculations
    first_edge = waveobs[0] - 0.5*wave_distance[0]
    last_edge = waveobs[-1] + 0.5*wave_distance[-1]
    # Build the final edges array
    edges = np.array([first_edge] + edges_tmp.tolist() + [last_edge])

    # Bin width
    bin_width = edges[1:] - edges[:-1]          # width per pixel

    # FWHM of the gaussian for the given resolution
    if from_resolution is None:
        # Convolve using instrumental resolution (smooth but not degrade)
        fwhm = waveobs / to_resolution
    else:
        # Degrade resolution
        fwhm = __get_fwhm(waveobs, from_resolution, to_resolution)
    sigma = __fwhm_to_sigma(fwhm)
    # Convert from wavelength units to bins
    fwhm_bin = fwhm / bin_width

    # Round number of bins per FWHM
    nbins = np.ceil(fwhm_bin) #npixels

    # Number of measures
    nwaveobs = len(waveobs)

    # In theory, len(nbins) == len(waveobs)
    for i in np.arange(len(nbins)):
        current_nbins = 2 * nbins[i] # Each side
        current_center = waveobs[i] # Center
        current_sigma = sigma[i]

        # Find lower and uper index for the gaussian, taking care of the current spectrum limits
        lower_pos = int(max(0, i - current_nbins))
        upper_pos = int(min(nwaveobs, i + current_nbins + 1))

        # Select only the flux values for the segment that we are going to convolve
        flux_segment = flux[lower_pos:upper_pos+1]
        err_segment = err[lower_pos:upper_pos+1]
        waveobs_segment = waveobs[lower_pos:upper_pos+1]

        # Build the gaussian corresponding to the instrumental spread function
        gaussian = np.exp(- ((waveobs_segment - current_center)**2) / (2*current_sigma**2)) / np.sqrt(2*np.pi*current_sigma**2)
        gaussian = gaussian / np.sum(gaussian)

        # Convolve the current position by using the segment and the gaussian
        if flux[i] > 0:
            # Zero or negative values are considered as gaps in the spectrum
            only_positive_fluxes = flux_segment > 0
            weighted_flux = flux_segment[only_positive_fluxes] * gaussian[only_positive_fluxes]
            current_convolved_flux = weighted_flux.sum()
            convolved_flux[i] = current_convolved_flux
        else:
            convolved_err[i] = 0.0

        if err[i] > 0:
            # * Propagate error Only if the current value has a valid error value assigned
            #
            # Error propagation considering that measures are dependent (more conservative approach)
            # because it is common to find spectra with errors calculated from a SNR which
            # at the same time has been estimated from all the measurements in the same spectra
            #
            weighted_err = err_segment * gaussian
            current_convolved_err = weighted_err.sum()
            #current_convolved_err = np.sqrt(np.power(weighted_err, 2).sum()) # Case for independent errors
            convolved_err[i] = current_convolved_err
        else:
            convolved_err[i] = 0.0

    return waveobs, convolved_flux, convolved_err


def convolve_spectrum(spectrum, to_resolution, from_resolution=None):
    """
    Spectra resolution smoothness/degradation.

    If "from_resolution" is not specified or its equal to "to_resolution", then the spectrum
    is convolved with the instrumental gaussian defined by "to_resolution".

    If "from_resolution" is specified, the convolution is made with the difference of
    both resolutions in order to degrade the spectrum.
    """
    if from_resolution is not None and from_resolution <= to_resolution:
        raise Exception("This method cannot deal with final resolutions that are bigger than original")

    waveobs, flux, err = __convolve_spectrum(spectrum['wl'], spectrum['flux'], spectrum['fluxerr'], to_resolution, from_resolution=from_resolution)
    convolved_spectrum = {}
    convolved_spectrum['wl'] = waveobs
    convolved_spectrum['flux'] = flux
    convolved_spectrum['fluxerr'] = err
    return convolved_spectrum


# function to interpolate spectrum
def interp_spectrum(wl_out, wl_in, flux_in, kind='cubic') :
    wl_in_copy = deepcopy(wl_in)
    
    # create interpolation function for input data
    f = interp1d(wl_in_copy, flux_in, kind=kind)
    
    # create mask for valid range of output vector
    mask = wl_out > wl_in[0]
    mask &= wl_out < wl_in[-1]
    
    flux_out = np.full_like(wl_out, np.nan)
    
    # interpolate data
    flux_out[mask] = f(wl_out[mask])
    return flux_out
