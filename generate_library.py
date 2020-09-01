# -*- coding: iso-8859-1 -*-

"""
    Created on July 7 2020
    
    Description: This routine generates a library of exoplanet atmosphere models using petitRADTRANS
    
    @author: Emilie Bessette <bessette@iap.fr>, Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python generate_library.py --libdir=./mini-lib/
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

import numpy as np

import models_lib
import exoatmoslib
import exoatmos_params
#import exoatmos_params_full_lib as exoatmos_params
import json

parser = OptionParser()
parser.add_option("-d", "--libdir", dest="libdir", help="Directory to save all files in the library",type='string',default='')
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with generate_library.py -h ")
    sys.exit(1)

if options.verbose:
    print('Directory to save model files in the library: ', options.libdir)

# Check if output directory exists, if not then create a new one
if not os.path.exists(options.libdir):
    print("Inexistent libdir path. Creating new directory:",options.libdir)
    os.makedirs(options.libdir)

# load user defined parameters
p = exoatmos_params.load_exoatmos_lib_parameters(options.libdir, variable_parameters=True)

# get variable parameters
variables = exoatmoslib.get_exoatmos_variables(p)

# Uncomment below to check variable parameters and ranges
#print(variables)

# set array of parameter dicts being one for each model in the library
p_array = exoatmoslib.get_parameters_array(p, variables, abundance_index=0)

"""
    # uncomment this part to check full set of parameters to input in each
    # model calculation for all library spectra
for par in p_array :
    print(p_array)
    print("\n")
"""

models = {}

for p_loc in p_array :
    if options.verbose :
        print("Generating model file:", p_loc['FILENAME'])
    
    # Calculate model using petitRADTRANS
    model = exoatmoslib.calculate_model(p_loc)

    # save model to FITS file
    exoatmoslib.save_to_fits(p_loc['FILENAME'], model)

    # append model file to database
    models[p_loc['FILENAME']] = models_lib.get_spectrum_info_from_fits(p_loc['FILENAME'])

# define database filename
dbfilename = "{}/db.json".format(options.libdir)
if options.verbose :
    print("Saving database into file:",dbfilename)

# flush all entries to the database
with open(dbfilename, 'w') as json_file:
    json.dump(models, json_file)
