# -*- coding: iso-8859-1 -*-
"""
    Created on Jul 06 2020
    
    Description: python library for exoplanetary quantities
    
    @author: Eder Martioli <martioli@iap.fr>
    
    Institut d'Astrophysique de Paris, France.
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

import numpy as np
import json


def get_exoplanet(dbfile, exoplanet="HD 189733 b") :
    """
        Function to return a given exoplanet entry in the
        json database obtained from
        https://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html
        """
    try :
        with open(dbfile, 'r') as f:
            datastore = json.load(f)
        for i in range(len(datastore)) :
            if datastore[i]['pl_name'] == exoplanet :
                return datastore[i]
        print("WARNING: could not find exoplanet ",exoplanet, " in database file",dbfile)
    except :
        print("ERROR: could not open exoplanet database file ",options.input)
        exit()
