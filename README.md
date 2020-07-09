# exoplanet-atmosphere
Library of exoplanet atmosphere models

The following depencies should be installed:

`numpy`, `scipy`, `astropy`, `matplotlib`, `optparse`, `json`, `copy`, ..

Below is an example to generate a model spectrum using petitRADTRANS:

```
python generate_spectrum.py --output='H2O_T1100_ab-3.5_Rp1.14_Mp1.14.fits' --teq=1100 --rp=1.14 --mp=1.14 --species='H2O CO2 CH4 CO' --abundances="-3. -6 -7.5 -3" -p
```

The following routine generates a library of models:

```
python generate_library.py --libdir=./mini-lib/
```

Here we make available a mini-version of a library to allow testing of routines to manage the library. The full library is stored at the IAP exo-atmos server.

One can test access to the library by the following command:

```
python select_model.py --db=mini-lib/db.json --teq=1100 --rp=1.14 --mp=1.14 --species=H2O --ab=-3
```

![Alt text](Figures/select_model_all.png?raw=true "Title")
