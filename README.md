# exoplanet-atmosphere
Library of exoplanet atmosphere models

The following depencies should be installed:

`numpy`, `scipy`, `astropy`, `matplotlib`, `optparse`, `json`, `copy`, ..

To generate a model spectrum using petitRADTRANS:

```
import exoatmoslib, exoatmos_params
import matplotlib.pyplot as plt 

p = exoatmos_params.load_exoatmos_lib_parameters()

# One can either use the default parameters above
# or change some of the parameters as follows:
# Equilibrium temperature [K]
p['TEQ'] = 1100
# Species for the atmospheric composition of models
p['SPECIES'] = ['H2O', 'CO2', 'CH4', 'CO']
# Log of relative molar fraction abundances: log([X]/[H2])
p['ABUNDANCES'] = [-3., -6, -7.5, -3]
# planet mass in [MJup]
p['MP'] = 1.142
# planet radius in [RJup]
p['RP'] = 1.138

# generate model
model = exoatmoslib.calculate_model(p)

# save model to FITS file
exoatmoslib.save_to_fits('model_example.fits', model)
```

The library of models is saved in the sub-directory `Model-library`. Here we make available only a mini-version of this library to allow testing of routines to manage the library. The full library has about XX Gb and it is stored at the IAP exo-atmos server.

To start using the library, one needs to generate a local database. Every time the model files are updated, it is required to generate again the database. The database can be created by running the following command:

```
python create_db.py --pattern=Model-library/*.fits --output=Model-library/db.json -v
```

Then one can test access to the library by the following command:

```
python select_model.py --db=Model-library/db.json --teq=1100 --rp=1.25 --mp=1.81 --wl0=950 --wlf=2500
```

![Alt text](Figures/select_model_all.png?raw=true "Title")
