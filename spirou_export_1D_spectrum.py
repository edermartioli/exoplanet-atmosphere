# -*- coding: iso-8859-1 -*-
"""
    Created on July 3 2020
    
    Description: This routine reads SPIRou reduced spectra and exports the data to 1D spectrum in *.fits format, suitable for the trasmission spectroscopy analysis
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    # example using pattern to select files with relative path:
    python spirou_export_1D_spectrum.py --epattern=Data/HD189733/*e.fits --rvfile=Data/HD189733/HD189733.rdb wl0=1600 wlf=1700 -p

    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys

import glob

import matplotlib.pyplot as plt
import numpy as np

import astropy.io.fits as fits
from astropy.io import ascii
from scipy import constants

import time

import spiroulib

parser = OptionParser()
parser.add_option("-e", "--epattern", dest="epattern", help="Spectral efits data pattern",type='string',default="*e.fits")
parser.add_option("-r", "--rvfile", dest="rvfile", help="Input file with RVs",type='string',default="")
parser.add_option("-1", "--wl0", dest="wl0", help="Initial wavelength [nm]",type='string',default="0.")
parser.add_option("-2", "--wlf", dest="wlf", help="Final wavelength [nm]",type='string',default="1e50")
parser.add_option("-o", action="store_true", dest="overwrite", help="overwrite output fits", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)
parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with spirou_export_1D_spectrum.py -h ")
    sys.exit(1)

if options.verbose:
    print('Spectral e.fits data pattern: ', options.epattern)
    print('Input file with RVs: ', options.rvfile)
    print('Initial wavelength [nm]: ', options.wl0)
    print('Final wavelength [nm]: ', options.wlf)

# make list of efits data files
if options.verbose:
    print("Creating list of e.fits spectrum files...")
inputedata = sorted(glob.glob(options.epattern))

dataset = spiroulib.load_array_of_spirou_spectra(inputedata, rvfile=options.rvfile, verbose=options.verbose, plot=options.plot)

spiroulib.save_spirou_spectra_to_fits(dataset, overwrite=options.overwrite, wl0=float(options.wl0), wlf=float(options.wlf), verbose=options.verbose)

