# -*- coding: iso-8859-1 -*-
"""
    Created on Jul 07 2020
    
    Description: python library for model exoplanet atmospheres using petitRADTRANS
    
    @authors:  Emilie Bessette <bessette@iap.fr>, Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import matplotlib.pyplot as plt
import numpy as np
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
import astropy.io.fits as fits
from copy import deepcopy

def save_to_fits(output, loc) :
    # create primary fits extension
    primary_hdu = fits.PrimaryHDU()
    
    # create pointer to primary header
    header = primary_hdu.header

    # Define x-axis coordinates (wavelength)
    header.set('CUNIT1', "MICRONS", 'x-axis wavelength units')
    header.set('CTYPE1', "WAVE")
    header.set('CDELT1', 'Evenly spaced log_10', 'wavelength step, can use numpy.geomspace()')
    header.set('CRVAL1', loc["WLEN_MICRONS"][0], 'initial wavelength')
    header.set('CRVALEND', loc["WLEN_MICRONS"][1], 'final wavelength')
    header.set('CNUM', len(loc['transmission']), 'number of samples to generate')
    header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    
    # set up primary header keywords
    header.set('ORIGIN', "petitRADTRANS v. 67b06e0d")
    header.set('RESOMODE', loc['MODE'], 'Resolution mode')
    header.set('P0', loc['P0'], 'Reference pressure [bar]')
    header.set('GAMMA', loc['GAMMA'], 'Ratio between optical and IR opacity')
    header.set('KAPPA', loc['KAPPA'], 'Atmospheric opacity coefficient in the IR')
    header.set('TINT', loc['TINT'], 'Planetary internal temperature [K]')

    for i in range(len(loc["CONTINUUM_OPACITIES"])) :
        keyword = 'CONT{0:03d}'.format(i)
        header.set(keyword, loc["CONTINUUM_OPACITIES"][i], 'Continuum opacities')

    for i in range(len(loc["RAYLEIGH_SPECIES"])) :
        keyword = 'RS_{0:03d}'.format(i)
        header.set(keyword, loc["RAYLEIGH_SPECIES"][i], 'Rayleigh scattering species')

    selecspc_index = loc['VARIABLE_SPECIES_INDEX']
    header.set('SELECSPC', loc['SPECIES'][selecspc_index], 'Selected species')

    for i in range(len(loc["SPECIES"])) :
        header.set(loc["SPECIES"][i], True, 'Presence of molecule')
        ab_keyword = 'AB_{0}'.format(loc["SPECIES"][i])
        header.set(ab_keyword, loc["ABUNDANCES"][i], 'Abundance of molecule w.r.t. H2, log10([X]/[H2])')
        op_keyword = 'OP_{0}'.format(loc["SPECIES"][i])
        header.set(op_keyword, loc["OPACITY_FILES"][i], 'Opacity file for molecule')

    if 'AB_He' not in header.keys() :
        header.set('AB_He', loc["AB_He"], 'Abundance of molecule [molar fraction]')

    if 'AB_H2' not in header.keys() :
        header.set('AB_H2', loc["AB_H2"], 'Abundance of molecule [molar fraction]')

    header.set('PSTART', loc["P_START"], 'First value of pressure range (logscale)')
    header.set('PSTOP', loc["P_STOP"], 'Last value of pressure range (logscale)')
    header.set('PSCALE', loc["P_SCALE"], 'pressure step, can use numpy.linspace()')
    header.set('PNUM', loc["P_NSAMPLES"], 'Number of samples to generate')
    header.set('TMOD', loc["P-T_Model"], 'Guillot P-T model')
    #varying parameters
    header.set('RPJUP', loc['RP'], 'Planetary radius [R_jupiter]')
    header.set('MPJUP', loc['MP'] , 'Planetary mass [M_jupiter]')
    header.set('RPCGS', loc['RP_CGS'], 'Planetary radius [cm]')
    header.set('MPCGS', loc['MP_CGS'] , 'Planetary mass [g]')
    header.set('GRAVITY', loc['GRAVITY_CGS'], 'Planetary gravity at R_pl [g.cm-2]')

    header.set('TEQ', loc['TEQ'], 'Planetary equilibrium temperature [K]')

    # set up transmission extension
    transm_header = fits.Header()
    transm_header.set('CUNIT1', "MICRONS", 'wavelength units')
    transm_header.set('CTYPE1', "WAVE")
    transm_header.set('CDELT1', 'Evenly spaced log_10', 'wavelength step')
    transm_header.set('CRVAL1', loc["WLEN_MICRONS"][0], 'initial wavelength')
    transm_header.set('CRVALEND', loc["WLEN_MICRONS"][1], 'final wavelength')
    transm_header.set('CNUM', len(loc['transmission']), 'number of samples to generate')
    transm_header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    transm_header.set('YAXISU', 'Transit radius [R_Jup]', 'unit for y-axis')
    hdu_transmission = fits.ImageHDU(data=loc['transmission'], name="TRANSMISSION", header=transm_header)
    
    # set up emission extension
    emiss_header = fits.Header()
    emiss_header.set('CUNIT1', "MICRONS", 'wavelength units')
    emiss_header.set('CTYPE1', "WAVE")
    emiss_header.set('CDELT1', 'Evenly spaced log_10', 'wavelength step')
    emiss_header.set('CRVAL1', loc["WLEN_MICRONS"][0], 'initial wavelength')
    emiss_header.set('CRVALEND', loc["WLEN_MICRONS"][1], 'final wavelength')
    emiss_header.set('CNUM', len(loc['emission']), 'number of samples to generate')
    emiss_header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    emiss_header.set('YAXISU', 'Planet flux [10e-6 erg cm^-2 s^-1 Hz^-1]', 'unit for y-axis')
    hdu_emission = fits.ImageHDU(data=loc['emission'], name="EMISSION", header=emiss_header)

    # bundle up all extensions in an HDU list
    mef_hdu = fits.HDUList([primary_hdu, hdu_transmission, hdu_emission])
    
    # write output
    mef_hdu.writeto(output, overwrite=True) #change to correct directory + filename


def opacity_files(species, mode='lbl') :
    
    loc = {}
    
    if mode == 'lbl':
        loc['H2O'] = 'H2O_main_iso'
        loc['CH4'] = 'CH4_main_iso'
        loc['CO'] = 'CO_main_iso'
        #loc['CO'] = 'CO_all_iso'
        loc['CO2'] = 'CO2_main_iso'
        loc['NH3'] = 'H2O_main_iso'
        loc['H2'] = 'H2_main_iso'
        loc['H2S'] = 'H2S_main_iso'
        loc['C2H2'] = 'C2H2_main_iso'
        loc['HCN'] = 'HCN_main_iso'
        loc['O3'] = 'O3_main_iso'
        loc['PH3'] = 'PH3_main_iso'
        loc['TiO'] = 'TiO_all_iso'
        loc['FeH'] = 'FeH_main_iso'

    elif mode == 'c-k':
        loc['C2H2'] = 'C2H2'
        loc['CH4'] = 'CH4'
        loc['CO'] = 'CO'
        loc['CO2'] = 'CO2'
        loc['FeH'] = 'FeH'
        loc['H2'] = 'H2'
        loc['H2O'] = 'H2O'
        loc['H2S'] = 'H2S'
        loc['HCN'] = 'HCN'
        loc['HDO'] = 'HDO'
        loc['NH3'] = 'NH3'
        loc['O3'] = 'O3'
        loc['OH'] = 'OH'
        loc['PH3'] = 'PH3'
        loc['TiO'] = 'TiO'
            
    loc['SiO'] = 'SiO_main_iso'
    loc['VO'] = 'VO'
    loc['K'] = 'K'
    loc['Na'] = 'Na'

    out_opacity_files = []
    
    for item in species :
        if item in loc:
            out_opacity_files.append(loc[item])
    
    return out_opacity_files


def calculate_model(loc) :
    
    '''
        Description: main wrapper function to run petitRADTRANS to calculate atmosphere model
        '''
    
    WLEN = [loc['WL0']/1000., loc['WLF']/1000.]
    loc["WLEN_MICRONS"] = WLEN
  
    loc["RP_CGS"] = loc['RP'] * nc.r_jup_mean #nominal hot Jupiter radius should be 1.25*nc.r_jup_mean
    loc["MP_CGS"] = loc['MP'] * nc.m_jup #nominal hot Jupiter mass should be 1.81*nc.r_jup_mean
    loc["GRAVITY_CGS"] = (nc.G * loc["MP_CGS"]) / (loc["RP_CGS"]**2)
    
    loc["OPACITY_FILES"] = opacity_files(loc['SPECIES'], mode=loc['MODE'])

    # Create the Ratrans object here
    atmosphere = Radtrans(line_species=loc["OPACITY_FILES"], \
                      rayleigh_species=loc['RAYLEIGH_SPECIES'], \
                      continuum_opacities=loc['CONTINUUM_OPACITIES'], \
                      mode=loc['MODE'], \
                      wlen_bords_micron = loc["WLEN_MICRONS"])

    # Create pressure and temperature distribution: Guillot
    loc["P_SCALE"] = 'log_10'
    pressures = np.logspace(loc["P_START"], loc["P_STOP"], loc["P_NSAMPLES"])

    loc["P-T_Model"] = 'nc.guillot_global()'
    temperature = nc.guillot_global(pressures, loc['KAPPA'], loc['GAMMA'], loc['GRAVITY_CGS'] , loc['TINT'], loc['TEQ'])

    atmosphere.setup_opa_structure(pressures)
    
    loc_abundances = {}
    # Define abundances of line species and mean molecular weight of atmosphere
    for i in range(len(loc["OPACITY_FILES"])) :
        ab_percent = 10**(loc['ABUNDANCES'][i] + np.log10(loc['AB_H2']))
        loc_abundances[loc["OPACITY_FILES"][i]] = ab_percent * np.ones_like(temperature)
    # Define abundances of H2 and He and mean molecular weight of atmosphere
    loc_abundances["He"] = loc['AB_He'] * np.ones_like(temperature)
    loc_abundances["H2"] = loc['AB_H2'] * np.ones_like(temperature)

    loc["loc_abundances"] = loc_abundances

    MMW = 2.3 * np.ones_like(temperature)
    
    # calculation of transmission spectrum
    atmosphere.calc_transm(temperature, loc_abundances, loc["GRAVITY_CGS"], MMW, R_pl=loc["RP_CGS"], P0_bar=loc['P0'])

    # calculation of emission spectrum
    atmosphere.calc_flux(temperature, loc_abundances, loc["GRAVITY_CGS"], MMW)

    # create transmission data vector
    loc['transmission'] = atmosphere.transm_rad / nc.r_jup_mean
    
    # create emission data vector
    loc['emission'] = atmosphere.flux / 1e-6
    
    # create wavelength array
    loc['wl'] = nc.c/atmosphere.freq/1e-4
    
    return loc


def get_exoatmos_variables(p) :

    variables = {}
    
    for var in p['VARIABLES'] :
        if var == 'ABUNDANCES' :
            if len(p['SPECIES']) != len(p['ABUNDANCES']) :
                print("ERROR: number of species must be equal to the number of input abundances")
                exit()
            
            for i in range(len(p['ABUNDANCES'])) :
                if type(p['ABUNDANCES'][i]) is list :
                    keyword = "AB_{}".format(p['SPECIES'][i])
                    variables[keyword] = p['ABUNDANCES'][i]
        else :
            if type(p[var]) is list :
                variables[var] = p[var]

    return variables


def get_parameters_array(p, variables, abundance_index=0) :
    # in this version it only considers the following variables:
    #  Teq, Abundance, Rp, Mp

    # Equilibrium temperature
    if 'TEQ' in variables :
        Teq_array = np.linspace(p['TEQ'][0], p['TEQ'][1], p['TEQ'][2])
    else :
        Teq_array = [p['TEQ']]

    # Radius
    if 'RP' in variables :
        rp_array = np.linspace(p['RP'][0], p['RP'][1], p['RP'][2])
    else :
        rp_array = [p['RP']]

    # Radius
    if 'MP' in variables :
        mp_array = np.linspace(p['MP'][0], p['MP'][1], p['MP'][2])
    else :
        mp_array = [p['MP']]

    # consider only abundance of first element
    keyword = "AB_{}".format(p['SPECIES'][abundance_index])
    if keyword in variables :
        ab_array = np.linspace(p['ABUNDANCES'][abundance_index][0], p['ABUNDANCES'][abundance_index][1], p['ABUNDANCES'][abundance_index][2])
    else :
        ab_array = [p['ABUNDANCES'][abundance_index]]


    total_nfiles = len(rp_array)*len(mp_array)*len(Teq_array)*len(ab_array)
    print("The library will contain a total of {} files.".format(total_nfiles))

    p_array = []

    for T_equ in Teq_array :
        for ab in ab_array :
            for rp in rp_array :
                for mp in mp_array :
                    p_loc = deepcopy(p)
                    output = "{0}/{1}_T{2:.0f}_ab{3:.2f}_Rp{4:.3f}_Mp{5:.3f}.fits".format(p['LIBDIR'],p['SPECIES'][abundance_index],T_equ,ab,rp,mp)
                    p_loc['FILENAME'] = output
                    
                    p_loc['TEQ'] = T_equ
                    p_loc['RP'] = rp
                    p_loc['MP'] = mp
                    p_loc['ABUNDANCES'][abundance_index] = ab

                    p_array.append(p_loc)
                    
                    print("Output model: ",output)

    return p_array

def plot_model(model, type="transmission") :
    
    if type == "transmission" :
        plt.plot(model['wl']*1000., model['transmission'], label=model["SPECIES"])
        plt.ylabel(r"Planet radius [R$_{\rm Jup}$]")
    elif type == "emission" :
        plt.plot(model['wl']*1000., model['transmission'], label=model["SPECIES"])
        plt.ylabel(r"Planet flux [10e-6 erg cm^-2 s^-1 Hz^-1]")
    else :
        print("ERROR: type not supported. Choose transmission or emission")
        exit()

    plt.legend()
    plt.xlabel(r"$\lambda$ [nm]")
    plt.show()
