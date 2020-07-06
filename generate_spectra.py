"""
    Created on July 3 2020
    
    Description: Example to generate transmission and emission spectra of HD 189733b 
    with 1 spectra containing all main molecules in the atmosphere, and 4 other
    spectra containing the 4 most abundant ones individually (CO, H2O, CH4, K)
    
    @author: Emilie Bessette <bessette@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import numpy as np
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
import astropy.io.fits as fits


WLEN = [0.95, 2.5]
R_pl = 1.138*nc.r_jup_mean #nominal hot Jupiter radius should be 1.25*nc.r_jup_mean
M_pl = 1.142*nc.m_jup #nominal hot Jupiter mass should be 1.81*nc.r_jup_mean
gravity = (nc.G*M_pl)/(R_pl**2)
P0 = 0.01
i = 1.0

# Create the Ratrans object here
atmosphere = Radtrans(line_species=['CH4_main_iso','CO_main_iso','H2O_main_iso','K','NH3_main_iso','Na','CO2_main_iso'], \
                    rayleigh_species=['H2','He'], \
                    continuum_opacities = ['H2-H2', 'H2-He'], \
                    mode='lbl', \
                    wlen_bords_micron = WLEN)

#%%
# Create pressure and temperature distribution: Guillot
pressures = np.logspace(-6, 2, 100)

kappa_IR = 0.01
gamma = 0.4
T_int = 500.
T_equ = 1100.

temperature = nc.guillot_global(pressures, kappa_IR, gamma, gravity, T_int, T_equ)

#%%
atmosphere.setup_opa_structure(pressures)

# Define abundances of line species and mean molecular weight of atmosphere
abundances = {}
abundances['CH4_main_iso'] = 2.6e-6 * np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['CO_main_iso'] = 4.5e-4 * np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['H2O_main_iso'] = 3.5e-4 * np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['K'] = 2.6e-6 * np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['NH3_main_iso'] = 6.0e-7 * np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['Na'] = 1.5e-7 * np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['CO2_main_iso'] = 1.5e-7 * np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['He'] = 0.19* np.ones_like(temperature)*(1-(i-1)*0.01556)
abundances['H2'] = 0.81* np.ones_like(temperature)*(1-(i-1)*0.01556)

MMW = 2.3 * np.ones_like(temperature) 

atmosphere.calc_transm(temperature, abundances, gravity, MMW, R_pl=R_pl, P0_bar=P0)

#%%

def save_to_fits(loc) :
    # create primary fits extension
    primary_hdu = fits.PrimaryHDU()
    
    # create pointer to primary header
    header = primary_hdu.header

    # Define x-axis coordinates (wavelength)
    header.set('CUNIT1', "MICRONS", 'x-axis wavelength units')
    header.set('CTYPE1', "WAVE")
    header.set('CDELT1', 'Evenly spaced log_10', 'wavelength step, can use numpy.geomspace()')
    header.set('CRVAL1', loc['wl0'], 'initial wavelength')
    header.set('CRVALEND', loc['wlf'], 'final wavelength')
    header.set('CNUM', len(atmosphere.transm_rad/nc.r_jup_mean), 'number of samples to generate')
    header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    
    # set up primary header keywords
    header.set('ORIGIN', "petitRADTRANS v. 67b06e0d")
    header.set('RESO', "lbl", 'Resolution mode')
    header.set('P0', P0, 'Reference pressure [bar]')
    header.set('gamma', gamma, 'Ratio between optical and IR opacity')
    header.set('kappa_IR', kappa_IR, 'Atmospheric opacity in the IR wavelengths')
    header.set('T_int', T_int, 'Planetary internal temperature [K]')
    
    if 'H2' in atmosphere.rayleigh_species: 
        header.set('H2', True, 'Presence of molecule')
        header.set('AB_H2', 0.81, 'Abundance of molecule times i [molar fraction]')
    if 'He' in atmosphere.rayleigh_species:
        header.set('He', True, 'Presence of molecule')
        header.set('AB_He', 0.19, 'Abundance of molecule times i [molar fraction]')
    if 'CO_main_iso' in atmosphere.line_species:
        header.set('CO', True, 'Presence of molecule')
        header.set('AB_CO', 4.5e-4, 'Abundance of molecule times i [molar fraction]')
    if 'H2O_main_iso' in atmosphere.line_species:
        header.set('H2O', True, 'Presence of molecule')
        header.set('AB_H2O', 3.5e-4, 'Abundance of molecule times i [molar fraction]')
    if 'CH4_main_iso' in atmosphere.line_species:
        header.set('CH4', True, 'Presence of molecule')
        header.set('AB_CH4', 2.6e-6, 'Abundance of molecule times i [molar fraction]')
    if 'NH3_main_iso' in atmosphere.line_species:
        header.set('NH3', True, 'Presence of molecule')
        header.set('AB_NH3', 6.0e-7, 'Abundance of molecule times i [molar fraction]')
    if 'CO2_main_iso' in atmosphere.line_species:
        header.set('CO2', True, 'Presence of molecule')
        header.set('AB_CO2', 1.5e-7, 'Abundance of molecule times i [molar fraction]')
    if 'Na' in atmosphere.line_species:
        header.set('Na', True, 'Presence of molecule')
        header.set('AB_Na', 1.5e-7, 'Abundance of molecule times i [molar fraction]')
    if 'K' in atmosphere.line_species:
        header.set('K', True, 'Presence of molecule')
        header.set('AB_K', 2.6e-6, 'Abundance of molecule times i [molar fraction]')
    
    header.set('CONT', "H2-H2, H2-He", 'Continuum opacities')
    header.set('PSTART', 1e-6, 'First value of pressure range')
    header.set('PSTOP', 1e2, 'Last value of pressure range')
    header.set('PSCALE', 'log_10', 'pressure step, can use numpy.linspace()')
    header.set('PNUM', 100, 'Number of samples to generate')
    header.set('TMOD', 'nc.guillot_global()', 'Guillot P-T model')
    
    #varying parameters
    header.set('R_pl', R_pl/nc.r_jup_mean, 'Planetary radius [R_jupiter]')
    header.set('M_pl', M_pl/nc.m_jup , 'Planetary mass [M_jupiter]')
    header.set('T_equ', T_equ, 'Planetary equilibrium temperature [K]')
    header.set('i', i, 'Intensity of abundance of molecule')

    # set up transmission extension
    transm_header = fits.Header()
    transm_header.set('CUNIT1', "MICRONS", 'wavelength units')
    transm_header.set('CTYPE1', "WAVE")
    transm_header.set('CDELT1', 'Evenly spaced log_10', 'wavelength step')
    transm_header.set('CRVAL1', loc['wl0'], 'initial wavelength')
    transm_header.set('CRVALEND', loc['wlf'], 'final wavelength')
    transm_header.set('CNUM', len(atmosphere.transm_rad/nc.r_jup_mean), 'number of samples to generate')
    transm_header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    transm_header.set('YAXISU', 'Transit radius [R_Jup]', 'unit for y-axis')
    hdu_transmission = fits.ImageHDU(data=loc['transmission'], name="TRANSMISSION", header=transm_header)
    
    # set up emission extension
    emiss_header = fits.Header()
    emiss_header.set('CUNIT1', "MICRONS", 'wavelength units')
    emiss_header.set('CTYPE1', "WAVE")
    emiss_header.set('CDELT1', 'Evenly spaced log_10', 'wavelength step')
    emiss_header.set('CRVAL1', loc['wl0'], 'initial wavelength')
    emiss_header.set('CRVALEND', loc['wlf'], 'final wavelength')
    emiss_header.set('CNUM', len(atmosphere.transm_rad/nc.r_jup_mean), 'number of samples to generate')
    emiss_header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    emiss_header.set('YAXISU', 'Planet flux [10e-6 erg cm^-2 s^-1 Hz^-1]', 'unit for y-axis')
    hdu_emission = fits.ImageHDU(data=loc['emission'], name="EMISSION", header=emiss_header)

    # bundle up all extensions in an HDU list
    mef_hdu = fits.HDUList([primary_hdu, hdu_transmission, hdu_emission])
    
    # write output
    mef_hdu.writeto('/filepath/filename.fits', overwrite=True) #change to correct directory + filename


# create dictionary container to save data
loc = {}

# define basic wavelength grid
loc['wl0'] = WLEN[0]
loc['wlf'] = WLEN[1]

# create transmission data vector
loc['transmission'] = atmosphere.transm_rad/nc.r_jup_mean
# create emission data vector
loc['emission'] = atmosphere.flux/1e-6

# save data to output fits
save_to_fits(loc)
