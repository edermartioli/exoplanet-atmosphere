# exoplanet-atmosphere
Library of exoplanet atmosphere models

The following depencies should be installed:

`numpy`, `scipy`, `astropy`, `matplotlib`, `optparse`, `json`, `copy`, ..

To generate a model spectrum using petitRADTRANS:

```
python HD189733b_transmission.py ...
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
