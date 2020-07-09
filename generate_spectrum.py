# -*- coding: iso-8859-1 -*-

"""
    Created on July 6 2020
    
    Description: This routine generates an exoplanet atmosphere model using petitRADTRANS
                 for a given set of input parameters
    
    @author: Emilie Bessette <bessette@iap.fr>, Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python generate_spectrum.py --output='H2O_T1100_ab-3.5_Rp1.14_Mp1.14.fits' --teq=1100 --rp=1.14 --mp=1.14 --species='H2O CO2 CH4 CO' --abundances="-3. -6 -7.5 -3" -p
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

import exoatmoslib, exoatmos_params
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-o", "--output", dest="output", help="Output model FITS filename",type='string',default='')
parser.add_option("-T", "--teq", dest="teq", help="Equilibrium temperature [K]",type='string',default='')
parser.add_option("-R", "--rp", dest="rp", help="Planet radius [RJup]",type='string',default='')
parser.add_option("-m", "--mp", dest="mp", help="Planet mass [MJup]",type='string',default='')
parser.add_option("-s", "--species", dest="species", help="Select molecule species",type='string',default='')
parser.add_option("-a", "--abundances", dest="abundances", help="Mass fraction of species",type='string',default='')
parser.add_option("-0", "--wl0", dest="wl0", help="Start wavelength",type='string',default='950.')
parser.add_option("-f", "--wlf", dest="wlf", help="Final wavelength",type='string',default='2500.')
parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with generate_spectrum.py -h ")
    sys.exit(1)

if options.verbose:
    print('Output model FITS filename: ', options.output)
    print('Equilibrium temperature [K]: ', options.teq)
    print('Planet radius [RJup]: ', options.rp)
    print('Planet mass [MJup]: ', options.mp)
    print('Select molecule species: ', options.species)
    print('Mass fraction of species: ', options.abundances)
    print('Start wavelength: ', options.wl0)
    print('Final wavelength: ', options.wlf)

# load default parameters
p = exoatmos_params.load_exoatmos_lib_parameters()

# One can either use the default parameters above
# or change some of the parameters as follows:
# Equilibrium temperature [K]
p['TEQ'] = float(options.teq)

# Species for the atmospheric composition of models
species = (options.species).split(" ")
p['SPECIES'] = species

# Log of relative molar fraction abundances: log([X]/[H2])
ab_strings = (options.abundances).split(" ")
abundances = []
for i in range(len(ab_strings)) :
    abundances.append(float(ab_strings[i]))
p['ABUNDANCES'] = abundances

# planet mass in [MJup]
p['MP'] = float(options.rp)
# planet radius in [RJup]
p['RP'] = float(options.mp)

# generate model
model = exoatmoslib.calculate_model(p)

if options.plot :
    # plot model
    exoatmoslib.plot_model(model, type="transmission")

if options.output != "" :
    # save model to fits
    exoatmoslib.save_to_fits(options.output, model)
