# -*- coding: iso-8859-1 -*-

"""
    Created on July 6 2020
    
    Description: This routine selects the exoplanet atmosphere model closest to the input parameters: Teq, Rp, Mp, Abundances..
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python select_model.py --db=Model-library/db.json --teq=1100 --rp=1.25 --mp=1.81 --species='h2o' --wl0=1200 --wlf=1600
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

best_model_filepath = models_lib.get_best_model_file(options.db, T_EQU=float(options.teq), AB=float(options.ab), R_PL=float(options.rp), M_PL=float(options.mp), species=options.species)

print("Selected model: ", best_model_filepath)

if options.species == 'all' :
    model_all = models_lib.load_spectrum_from_fits(best_model_filepath, wlmin=float(options.wl0), wlmax=float(options.wlf))
    plt.plot(model_all['wl'], model_all['transmission'], label="All")

    colors = models_lib.species_colors()
    species_list = ['H2O', 'CO', 'CH4', 'CO2']
    
    for species in species_list :
        model_filepath = models_lib.get_best_model_file(options.db, T_EQU=float(options.teq), R_PL=float(options.rp), M_PL=float(options.mp), species=species)
        model_species = models_lib.load_spectrum_from_fits(model_filepath, wlmin=float(options.wl0), wlmax=float(options.wlf))
        plt.plot(model_species['wl'], model_species['transmission'], color=colors[species], label=species, alpha=0.5)
else :
    model = models_lib.load_spectrum_from_fits(best_model_filepath, wlmin=float(options.wl0), wlmax=float(options.wlf))
    plt.plot(model['wl'], model['transmission'], label=options.species)

plt.legend()
plt.xlabel(r"wavelength (nm)")
plt.ylabel(r"transit radius R$_{\rm Jup}$")
#plt.ylabel(r"relative transmission")
plt.show()
