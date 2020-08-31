# -*- coding: iso-8859-1 -*-

"""
    Created on July 6 2020
    
    Description: This routine shows an example on how to obtain an exoplanet atmosphere model
                 from our exoatmos library given the input parameters: Teq, Rp, Mp, abundance, species.
                 One can use either the nearest approach, which select an existing model with
                 parameters closest to the input parameters, or one can use interpolated models,
                 which gets the closest Rp, Mp, and abundances, and performs a linear interpolation
                 between the two closest Teq values.
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python select_model.py --db=mini-lib/db.json --teq=1132 --ab=-3.4 --rp=1.25 --mp=1.81 --species='H2O' --wl0=1200 --wlf=1600 -vn
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

import models_lib
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-d", "--db", dest="db", help="Input database json file",type='string',default='')
parser.add_option("-T", "--teq", dest="teq", help="Equilibrium temperature [K]",type='string',default='')
parser.add_option("-a", "--ab", dest="ab", help="Abundance of species X, log([X]/[H2])",type='string',default='')
parser.add_option("-R", "--rp", dest="rp", help="Planet radius [RJup]",type='string',default='')
parser.add_option("-m", "--mp", dest="mp", help="Planet mass [MJup]",type='string',default='')
parser.add_option("-s", "--species", dest="species", help="Select molecule species",type='string',default='all')
parser.add_option("-0", "--wl0", dest="wl0", help="Start wavelength",type='string',default='900.')
parser.add_option("-f", "--wlf", dest="wlf", help="Final wavelength",type='string',default='2500.')
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-n", action="store_true", dest="nearest", help="Get nearest model instead of interpolation in Teq", default=False)
parser.add_option("-p", action="store_true", dest="plot", help="Plot selected spectrum", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with select_model.py -h ")
    sys.exit(1)

if options.verbose:
    print('Input database json file: ', options.db)
    print('Equilibrium temperature [K]: ', options.teq)
    print('Abundance of species X, log([X]/[H2]): ', options.ab)
    print('Planet radius [RJup]: ', options.rp)
    print('Planet mass [MJup]: ', options.mp)
    print('Select molecule species: ', options.species)
    print('Start wavelength: ', options.wl0)
    print('Final wavelength: ', options.wlf)

if options.nearest :
    best_model_filepath = models_lib.get_best_model_file(options.db, T_EQU=float(options.teq), AB=float(options.ab), R_PL=float(options.rp), M_PL=float(options.mp), species=options.species)
    print("Selected model: ", best_model_filepath)
    model = models_lib.load_spectrum_from_fits(best_model_filepath, wlmin=float(options.wl0), wlmax=float(options.wlf))
else :
    model = models_lib.get_interpolated_model(options.db, T_EQU=float(options.teq), AB=float(options.ab), R_PL=float(options.rp), M_PL=float(options.mp), species=options.species, wlmin=float(options.wl0), wlmax=float(options.wlf), return_wl=True, return_emission=False)
    print("Model calculated by linear interpolation of the following selected models: ", model['T_EQU_low_path'],"and", model['T_EQU_upp_path'])

if options.plot :
    plt.plot(model['wl'], model['transmission'], label=options.species)
    plt.legend()
    plt.xlabel(r"wavelength (nm)")
    plt.ylabel(r"transit radius R$_{\rm Jup}$")
    #plt.ylabel(r"relative transmission")
    plt.show()
