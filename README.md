# exoplanet-atmosphere
Library of exoplanet atmosphere models

The following depencies should be installed:

`numpy`, `scipy`, `astropy`, `matplotlib`, `optparse`, `json`, `copy`, `petitRADTRANS`, ...

Below is an example to generate a model spectrum using petitRADTRANS:

```
python generate_spectrum.py --output='H2O_T1100_ab-3.5_Rp1.14_Mp1.14.fits' --teq=1100 --rp=1.14 --mp=1.14 --species='H2O CO2 CH4 CO' --abundances="-3. -6 -7.5 -3" -p
```
![Alt text](Figures/generate_spectrum_example.png?raw=true "Title")

The following routine generates a library of models using the parameters in `exoatmos_params.py`:

```
python generate_library.py --libdir=./mini-lib/
```

Here we make available a mini-version of the library to allow testing of routines to manage the library. The full library is stored at the IAP exo-atmos server.

One can test access to the library by the following command:

```
python select_model.py --db=mini-lib/db.json --teq=1132 --ab=-3.4 --rp=1.25 --mp=1.81 --species='H2O' -vn
```

which returns the following:

```
Selected model:  mini-lib/H2O_T1100_ab-3.00_Rp1.138_Mp1.142.fits
```

![Alt text](Figures/select_model_H2O.png?raw=true "Title")
