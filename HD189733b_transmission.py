Python 2.7.12 (default, Aug 22 2019, 16:36:40) 
[GCC 5.4.0 20160609] on linux2
Type "copyright", "credits" or "license()" for more information.
>>> """
    Created on July 3 2020
    
    Description: Example to generate transmission spectra of HD 189733b 
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
R_pl = 1.25*nc.r_jup_mean #radius of HD189733b is 1.138*R_J
M_pl = 1.81*nc.m_jup #mass of HD189733b is 1.142*M_J
gravity = (nc.G*M_pl)/(R_pl**2)
P0 = 0.01
i = 1.0

# Create the Ratrans object here
atmosphere = Radtrans(line_species=['CH4_main_iso','CO_main_iso','H2O_main_iso','K','NH3_main_iso','Na','CO2_main_iso'], \
                    rayleigh_species=['H2','He'], \
                    continuum_opacities = ['H2-H2', 'H2-He'], \
                    mode='lbl', \
                    wlen_bords_micron = WLEN)

# atmosphere = Radtrans(line_species=['CO_main_iso'], \
#                     rayleigh_species=['H2','He'], \
#                     continuum_opacities = ['H2-H2', 'H2-He'], \
#                     mode='lbl', \
#                     wlen_bords_micron = WLEN)

# atmosphere = Radtrans(line_species=['H2O_main_iso'], \
#                     rayleigh_species=['H2','He'], \
#                     continuum_opacities = ['H2-H2', 'H2-He'], \
#                     mode='lbl', \
#                     wlen_bords_micron = WLEN)
    
# atmosphere = Radtrans(line_species=['CH4_main_iso'], \
#                     rayleigh_species=['H2','He'], \
#                     continuum_opacities = ['H2-H2', 'H2-He'], \
#                     mode='lbl', \
#                     wlen_bords_micron = WLEN)
    
# atmosphere = Radtrans(line_species=['K'], \
#                     rayleigh_species=['H2','He'], \
#                     continuum_opacities = ['H2-H2', 'H2-He'], \
#                     mode='lbl', \
#                     wlen_bords_micron = WLEN)

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
abundances['CH4_main_iso'] = 2.6e-6 * np.ones_like(temperature)*i
abundances['CO_main_iso'] = 4.5e-4 * np.ones_like(temperature)*i
abundances['H2O_main_iso'] = 3.5e-4 * np.ones_like(temperature)*i
abundances['K'] = 2.6e-6 * np.ones_like(temperature)*i
abundances['NH3_main_iso'] = 6.0e-7 * np.ones_like(temperature)*i
abundances['Na'] = 1.5e-7 * np.ones_like(temperature)*i
abundances['CO2_main_iso'] = 1.5e-7 * np.ones_like(temperature)*i
abundances['He'] = 0.19* np.ones_like(temperature)*i
abundances['H2'] = 0.85* np.ones_like(temperature)*i

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
    header.set('CNUM', 968, 'number of samples to generate')
    header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    
    # set up primary header keywords
    header.set('ORIGIN', "petitRADTRANS v. 67b06e0d")
    header.set('RESO', "lbl", 'Resolution mode')
    header.set('P0', 0.01, 'Reference pressure [bar]')
    header.set('gamma', 0.4, 'Ratio between optical and IR opacity')
    header.set('kappa_IR', 0.01, 'Atmospheric opacity in the IR wavelengths')
    header.set('T_int', 500, 'Planetary internal temperature [K]')
    
    if 'H2' in atmosphere.rayleigh_species: 
        header.set('H2', True, 'Presence of molecule')
        header.set('AB_H2', 0.85, 'Abundance of molecule times k [molar fraction]')
    if 'He' in atmosphere.rayleigh_species:
        header.set('He', True, 'Presence of molecule')
        header.set('AB_He', 0.19, 'Abundance of molecule times k [molar fraction]')
    if 'CO_main_iso' in atmosphere.line_species:
        header.set('CO', True, 'Presence of molecule')
        header.set('AB_CO', 4.5e-4, 'Abundance of molecule times k [molar fraction]')
    if 'H2O_main_iso' in atmosphere.line_species:
        header.set('H2O', True, 'Presence of molecule')
        header.set('AB_H2O', 3.5e-4, 'Abundance of molecule times k [molar fraction]')
    if 'CH4_main_iso' in atmosphere.line_species:
        header.set('CH4', True, 'Presence of molecule')
        header.set('AB_CH4', 2.6e-6, 'Abundance of molecule times k [molar fraction]')
    if 'NH3_main_iso' in atmosphere.line_species:
        header.set('NH3', True, 'Presence of molecule')
        header.set('AB_NH3', 6.0e-7, 'Abundance of molecule times k [molar fraction]')
    if 'CO2_main_iso' in atmosphere.line_species:
        header.set('CO2', True, 'Presence of molecule')
        header.set('AB_CO2', 1.5e-7, 'Abundance of molecule times k [molar fraction]')
    if 'Na' in atmosphere.line_species:
        header.set('Na', True, 'Presence of molecule')
        header.set('AB_Na', 1.5e-7, 'Abundance of molecule times k [molar fraction]')
    if 'K' in atmosphere.line_species:
        header.set('K', True, 'Presence of molecule')
        header.set('AB_K', 2.6e-6, 'Abundance of molecule times k [molar fraction]')
    
    header.set('CONT', "H2-H2, H2-He", 'Continuum opacities')
    header.set('PRANGE', 'np.logspace(-6,2,100)', 'Pressure range for Guillot')
    header.set('TMOD', 'nc.guillot_global()', 'Guillot P-T model')
    
    #varying parameters
    header.set('R_pl', 1.25, 'Planetary radius [R_jupiter]')
    header.set('M_pl', 1.81, 'Planetary mass [M_jupiter]')
    header.set('T_equ', 1100, 'Planetary equilibrium temperature [K]')
    header.set('i', 1, 'Intensity of abundance of molecule')

    # set up transmission extension
    transm_header = fits.Header()
    transm_header.set('CUNIT1', "MICRONS", 'wavelength units')
    transm_header.set('CTYPE1', "WAVE")
    transm_header.set('CDELT1', 'Evenly spaced log_10', 'wavelength step')
    transm_header.set('CRVAL1', loc['wl0'], 'initial wavelength')
    transm_header.set('CRVALEND', loc['wlf'], 'final wavelength')
    transm_header.set('CNUM', 968, 'number of samples to generate')
    transm_header.set('CAXIS', 'nc.c/atmosphere.freq/1e-4', 'expression for x-axis')
    hdu_transmission = fits.ImageHDU(data=loc['transmission'], name="TRANSMISSION", header=transm_header)

    # bundle up all extensions in an HDU list
    mef_hdu = fits.HDUList([primary_hdu, hdu_transmission])
    
    # write output
    mef_hdu.writeto('/filepath/filename.fits', overwrite=True) #change to correct directory + filename


# create dictionary container to save data
loc = {}

# define basic wavelength grid
loc['wl0'] = WLEN[0]
loc['wlf'] = WLEN[1]

# create transmission data vector
loc['transmission'] = atmosphere.transm_rad/nc.r_jup_mean

# save data to output fits
save_to_fits(loc)
