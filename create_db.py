# -*- coding: iso-8859-1 -*-

"""
    Created on Jul 5 2020
    
    Description: This routine creates a data base for the library of exoplanet atmosphere models
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python create_db.py --pattern=Model-library/*.fits --output=Model-library/db.json -v
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

import glob
import json

import models_lib

parser = OptionParser()
parser.add_option("-p", "--pattern", dest="pattern", help="Input data pattern",type='string',default='')
parser.add_option("-o", "--output", dest="output", help="Output json database file",type='string',default='')
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with create_db.py -h ")
    sys.exit(1)

if options.verbose:
    print('Input data pattern: ', options.pattern)
    print('Output json database file: ', options.output)

# make list of data files
if options.verbose:
    print("Loading list of spectrum models ...")
inputdata = sorted(glob.glob(options.pattern))

models = {}
for model_file in inputdata :
    filename = os.path.basename(model_file)
    if options.verbose :
        print("Loding model file:", filename)
    models[filename] = models_lib.get_spectrum_info_from_fits(model_file)

if options.verbose :
    print("Saving database into file:", options.output)

with open(options.output, 'w') as json_file:
    json.dump(models, json_file)
