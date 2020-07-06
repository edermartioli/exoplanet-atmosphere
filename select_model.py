# -*- coding: iso-8859-1 -*-

"""
    Created on July 6 2020
    
    Description: This routine selects the exoplanet atmosphere model closest to the input parameters: Teq, Rp, Mp, Abundances..
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python select_model.py --db=Model-library/db.json --teq=1100 --rp=1.25 --mp=1.81
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

import models_lib

parser = OptionParser()
parser.add_option("-d", "--db", dest="db", help="Input database json file",type='string',default='')
parser.add_option("-T", "--teq", dest="teq", help="Equilibrium temperature [K]",type='string',default='')
parser.add_option("-R", "--rp", dest="rp", help="Planet radius [RJup]",type='string',default='')
parser.add_option("-m", "--mp", dest="mp", help="Planet mass [MJup]",type='string',default='')
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with create_stellar_model_db.py -h ")
    sys.exit(1)

if options.verbose:
    print('Input database json file: ', options.db)
    print('Equilibrium temperature [K]: ', options.teq)
    print('Planet radius [RJup]: ', options.rp)
    print('Planet mass [MJup]: ', options.mp)

best_model_filepath = models_lib.get_best_model_file(options.db, T_EQU=float(options.teq), R_PL=float(options.rp), M_PL=float(options.mp))

print("Selected model: ", best_model_filepath)
