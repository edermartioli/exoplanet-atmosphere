# exoplanet-atmosphere

Welcome to the exoatmos_lib GitHub page. 
This README file will serve as a small documentation, and the Python files needed to generate and use the library  of exoplanet atmosphere models are present above. 

The exoatmos_lib library consists of modelled transmission and emission spectra of exoplanets, for clear atmospheres, for modulable ranges of equilibrium temperature, presence and abundance of molecular species, planetary radius and mass. The spectra in the library are of high (λ/Δλ=10^6) resolution and are modelled using the `petitRADTRANS` package.

This GitHub page contains the following: 
- mini-lib: folder containing 9 FITS files (with both transmission and emission spectra data inside) modelled in low resolution mode (or 'c-k' mode of (λ/Δλ=1000)).
- mini-lib_lbl: folder containing 9 FITS files (with both transmission and emission spectra data inside) modelled in high resolution mode (or 'lbl' mode of (λ/Δλ=10^6)).



GUIDE:
The following depencies should be installed:

`numpy`, `scipy`, `astropy`, `matplotlib`, `optparse`, `json`, `copy`, `petitRADTRANS`, `glob`, `os, sys`.

Below is an example to generate and plot a model spectrum using petitRADTRANS:

```
python generate_spectrum.py --output='H2O_T1100_ab-3.5_Rp1.14_Mp1.14.fits' --teq=1100 --rp=1.14 --mp=1.14 --species='H2O CO2 CH4 CO' --abundances="-3. -6 -7.5 -3" -p
```
![Alt text](Figures/generate_spectrum_example.png?raw=true "Title")

The following routine generates a library called "mini-lib" of models using the parameters in `exoatmos_params.py`:

```
python generate_library.py --libdir=./mini-lib/
```

Here we make a mini-version of the library available, to allow for testing of routines to manage the library. The full library is stored on the IAP exoatmos server at the filepath '/mnt/globalNS/tmp/bessette/exoatmos_lib/'.

One can test access to the mini-library by the following command:

```
python select_model.py --db=mini-lib/db.json --teq=1132 --ab=-3.4 --rp=1.25 --mp=1.81 --species='H2O' -vn
```

which returns the following:

```
Selected model:  mini-lib/H2O_T1100_ab-3.00_Rp1.138_Mp1.142.fits
```

![Alt text](Figures/select_model_H2O.png?raw=true "Title")
