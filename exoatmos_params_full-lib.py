#---------------------------------------------------------------------
#   define Library parameters
# -----------------------------------------------------------------------------
def load_exoatmos_lib_parameters(libdir="", variable_parameters=False) :
    
    #initialize parameters dictionary
    p = {}
    
    ####### SET LIBRARY DIRECTORY ############
    p['LIBDIR'] = libdir
    #-----------------------------------------------------

    ##################### VARIABLES ######################
    # Define list of parameters considered as POSSIBLE variables
    p['VARIABLES'] = ['TEQ', 'AB_He', 'AB_H2', 'ABUNDANCES', 'MP', 'RP', 'TINT', 'P0', 'KAPPA', 'GAMMA']
    # Note: in fact only TEQ, ABUNDANCE (for 1 species), RP, and MP are variables
    #       but the list above may be used in future implementations to allow
    #       other variables in the library grid.
    
    # Set index below for the selected species which the abundance is variable
    p['VARIABLE_SPECIES_INDEX'] = 0
    # Note: the index above is because we consider only one species to vary abundance at a time
    #       One can include more species in the model, but the abundances must be constant
    #-----------------------------------------------------
    
    ############### RESOLUTION MODE ######################
    # Define resolution mode:
    # high-resolution: p['mode'] = 'lbl'
    p['MODE'] = 'lbl'
    # low-resolution: p['mode'] = 'c-k'
    #p['MODE'] = 'c-k'
    #-----------------------------------------------------

    ########### SPECTRAL RANGE (in nm) ####################
    p['WL0'] = 950.
    p['WLF'] = 2500.
    #-----------------------------------------------------
    
    ########### EQUILIBRIUM TEMPERATURE [K] ###############
    # Equilibrium temperature [K]
    # For ranges insert 3 values: [ini, end, nsamples]
    if variable_parameters :
        p['TEQ'] = [500, 2500, 21]
    else :
        # For constant temperature insert a single value: Teq
        p['TEQ'] = 1100
    ### Follow the same format for all variables below
    #-----------------------------------------------------
    
    ######### CHEMICAL COMPOSITION #############
    # Define abundance of Helium
    p['AB_He'] = 0.24
    # Define abundance of H2
    p['AB_H2'] = 0.74
    
    # Whether or not to generate separate models for each species
    p['SEPARATE_SPECIES'] = True
    # Species for the atmospheric composition of models
    #p['SPECIES'] = ['H2O', 'CO2', 'CO', 'CH4'] # to include other species in the model
    #p['SPECIES'] = ['CO2']
    #p['SPECIES'] = ['CO']
    #p['SPECIES'] = ['CH4']
    p['SPECIES'] = ['H2O']
    
    # Log of molar fraction abundances: log([X]/[H2])
    if variable_parameters :
        #p['ABUNDANCES'] = [[-6,-2, 3], -6, -3, -8] # to include other species in the model
        #p['ABUNDANCES'] = [[-10,-2,10]]
        #p['ABUNDANCES'] = [[-10,-2,10]]
        #p['ABUNDANCES'] = [[-10,-2,10]]
        p['ABUNDANCES'] = [[-8,-2,13]]
    else :
        p['ABUNDANCES'] = [-3.5]
    #############################################
    # Define species included in the Rayleigh scattering
    p['RAYLEIGH_SPECIES'] = ['H2','He']
    # Define pairs of molecules contributing to the continuum opacities
    p['CONTINUUM_OPACITIES'] = ['H2-H2', 'H2-He']
    #-------------------------------------------------------
    
    ######### PLANET PHYSICAL PARAMETERS ##########
    # planet mass in [MJup]
    p['MP'] = 1.814
    # planet radius in [RJup]
    p['RP'] = 1.260
    # planet internal temperature in [K]
    p['TINT'] = 500.
    # atmospheric pressure [bar] at "surface" defined by planet radius
    p['P0'] = 0.01
    # cross section.  in the near-infrared kappa_IR=0.01
    p['KAPPA'] = 0.01
    # Define gamma, which is a free parameter for the T-P profile adopted in the Guillot model
    p['GAMMA'] = 0.4
    # Define pressure range (log scale) to calculate T-P profile
    p['P_START'] = -6
    p['P_STOP'] = 2
    p['P_NSAMPLES'] = 100
    #-------------------------------------------------------
    
    return p
