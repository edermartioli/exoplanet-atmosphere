# exoplanet-atmosphere
Library of exoplanet atmosphere models

The following depencies should be installed:

`numpy`, `scipy`, `astropy`, `matplotlib`, `optparse`, `copy`, ..

To generate a model spectrum:

```
python ../../spirou_pol.py --exp1=2329699e.fits --exp2=2329700e.fits --exp3=2329701e.fits 
--exp4=2329702e.fits --lsdmask=../../lsd_masks/marcs_t5000g50_all --output=2329699_pol.fits 
--output_lsd=2329699_lsd.fits -p -s -L
```
