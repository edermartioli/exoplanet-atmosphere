"""
    Created on July 6 2020
    
    Description: Example to interpolate transmission and emission spectra of HD 189733b 
    with 1 spectra containing all main molecules in the atmosphere, and 4 other
    spectra containing the 4 most abundant ones individually (CO, H2O, CH4, K). 
    Multi-variable interpolation with variables: T_equ, i, R_pl, M_pl
    
    @author: Emilie Bessette <bessette@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import numpy as np
from scipy.interpolate import griddata 

import astropy.io.fits as fits
import numpy as np
import json
from scipy.interpolate import griddata 

import os, sys

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


def interpolate_model_linear(y1, y2, T1, T2, T_EQU):
    """
        Get interpolated spectrum for unique T_EQU. Can input transmission
        or emission data into parameters y1 and y2, to get corresponding transmission
        or emission interpolated spectrum.
        
        :param y1: numpy array, y-axis values of lower bound
        :param y2: numpy array, y-axis values of upper bound
        :param T1: float, lower bound equilibrium temperature [K]
        :param T2: float, upper bound equilibrium temperature [K]
        :param T_EQU: float, input equilibrium temperature [K]
        :return model_interp: numpy array, y-axis values interpolated
        """
    
    points = np.array([T1, T2])
    bounds = np.array([y1, y2])
    model_interp = griddata(points, bounds, (T_EQU), method='linear')
    return model_interp

def interpolate_model(modeldb, T_EQU, AB, R_PL, M_PL, species, wlmin=900., wlmax=2500.) :
    """
        Get file path of best model in input database
        
        :param modeldb: string, json database file name
        :param T_EQU: float, input equilibrium temperature [K]
        :param AB: float, input abundance for selected species log([X]/[H2])
        :param R_PL: float, input planet radius [RJup]
        :param M_PL: float, input planet mass [MJup]
        :param species: string, input molecular species
        :return loc: dictionary with interpolated results
        """

    # load json database file containig all models in the library
    try :
        with open(modeldb, 'r') as f:
            datastore = json.load(f)
    except :
        print("ERROR: could not open models database file ",modeldb)
        exit()

    mind_ab = 1.e20
    mind_R_PL, mind_M_PL = 1.e20, 1.e20

    minkey = None

    # loop over all entried in the database to get best matched model file
    for key in datastore.keys() :
        
        if species == 'all' :
            ab_key = 'AB_{0}'.format(datastore[key]['SELECSPC'])
        else :
            ab_key = 'AB_{0}'.format(species)
    
        d_ab = np.abs(datastore[key][ab_key] - AB)
        d_R_PL = np.abs(datastore[key]['RPJUP'] - R_PL)
        d_M_PL = np.abs(datastore[key]['MPJUP'] - M_PL)
    
        if d_ab <= mind_ab and d_R_PL <= mind_R_PL and d_M_PL <= mind_M_PL and (species in datastore[key]['filename']):
            mind_ab = d_ab
            mind_R_PL, mind_M_PL = d_R_PL, d_M_PL
            minkey = key

    nearest_ab = datastore[minkey]['AB_'+species]
    nearest_R_PL = np.around(datastore[minkey]['RPJUP'],3)
    nearest_M_PL = np.around(datastore[minkey]['MPJUP'],3)
    
    minkey_new=[]
    for key in datastore.keys():
        if datastore[key]['AB_'+species] == nearest_ab and np.around(datastore[key]['RPJUP'],3) == nearest_R_PL \
            and np.around(datastore[key]['MPJUP'],3) == nearest_M_PL:
            minkey_new.append(key)
            
    for i, key in enumerate(minkey_new):
        if datastore[key]['TEQ']>T_EQU:
            T_EQU_upp_path = key
            T_EQU_upp = datastore[T_EQU_upp_path]['TEQ']
            T_EQU_low_path = minkey_new[i-1]
            T_EQU_low = datastore[T_EQU_low_path]['TEQ']
            break
    
    hdu_low=fits.open(T_EQU_low_path)  
    hdu_upp=fits.open(T_EQU_upp_path)  
    model_interp_transm = interpolate_model_linear(hdu_low[1].data, hdu_upp[1].data, T_EQU_low, T_EQU_upp, T_EQU)
    model_interp_emiss = interpolate_model_linear(hdu_low[2].data, hdu_upp[2].data, T_EQU_low, T_EQU_upp, T_EQU)
    
    hdr = hdu_low[0].header
    wl_0 = hdr["CRVAL1"]
    wl_f = hdr["CRVALEND"]
    wl_step = hdr["CDELT1"]
    wl_num = len(hdu_low["TRANSMISSION"].data)
    
    # create wavelength array
    wl = np.geomspace(wl_0*1000., wl_f*1000., wl_num)
    
    wlmask = np.where(np.logical_and(wl > wlmin, wl < wlmax))
    loc = get_spectrum_info_from_fits(T_EQU_low_path)  
    loc['wl'] = wl[wlmask]
    loc['TEQ'] = T_EQU

    if "TRANSMISSION" in hdu_low :
        transmission = model_interp_transm
        loc['transmission'] = transmission[wlmask]
    
    if "EMISSION" in hdu_low :
        emission = model_interp_emiss
        loc['emission'] = emission[wlmask]

    return loc
