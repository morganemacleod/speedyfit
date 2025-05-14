"""
MM: 5/14/25
Created to group fits files from
https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas
into a single fits that can be used by speedyfit
"""

from glob import glob
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits



### NOTE: CHANGE THESE TO USE DIFFERENT FOLDER
file_list = glob("/Users/morganmacleod/CODE/speedyfit/sed_models/ckm05/*.fits")
filename_out = "sed_models/castellikurucz_z-0.5_seds.fits"


def logg_from_str(gg):
    "eg g05 = 0.5"
    return np.round(float(gg[1:])*0.1,1)


hdu_list = []

primary_hdu = fits.PrimaryHDU()
primary_hdu.header['ORIGIN'] = 'MM, 5/14/25, aggregating STIS files from Castelli and Kurucz'
primary_hdu.header['HISTORY'] = "File created by F.R.Boffi ATLAS9 model atmospheres by Castelli and Kurucz (2004). Wavelength is in Angstrom. Fluxes tabulated in units of erg/s/cm^2/A (after converting original units into flam, as described in README file and the SYNPHOT manual) and are surface fluxes. To transform to observed fluxes multiply by (R/D)^2 where R is the radius of the star and D the distance. Each column in the table represents the spectrum of a star for the same metallicity and effective temperature but different gravity. Gravities range from log_g = +0.0 (g00 in the column header) to log_g = +5.0 (g50 in the column header)."

for i,ff in enumerate(file_list):
    hdu = fits.open(ff)
    # fill in the primary header
    if i==0:
        primary_hdu.header['LOG_Z'] = hdu[0].header['LOG_Z']
        hdu_list.append(primary_hdu)
    
    # header
    teff = hdu[0].header['TEFF']
    
    # extract the data columns
    wave = hdu[1].data['WAVELENGTH']
    
    
    for j,gg in enumerate(['g00','g05','g10','g15','g20','g25','g30','g35','g40','g45','g50']):
        logg = logg_from_str(gg)
        flux = hdu[1].data[gg]

        # only include ones with data
        if np.max(flux)>0:
            # make an astropy table
            table = Table([wave,flux],names=['wavelength','flux'])
            # Create a new BinTableHDU
            table_hdu = fits.BinTableHDU(data=table, name=f'T{teff:05}_logg{logg:03}')
            print(f'T{teff:05}_logg{logg:03}')
            # Add metadata to this HDU's header
            table_hdu.header['TEFF'] = teff
            table_hdu.header['LOGG'] = logg
            table_hdu.header['SRCFILE'] = ff.split('/')[-1]

            hdu_list.append(table_hdu)

        
        plt.loglog(wave,flux)
    hdu.close()
    
# Create the full HDUList and write to a new file
hdulist = fits.HDUList(hdu_list)
# Assume hdulist is your HDUList with a PrimaryHDU + multiple BinTableHDUs
primary_hdu = hdulist[0]  # Keep the primary HDU
# Sort the rest (BinTableHDUs) by HDU name
sorted_bintable_hdus = sorted(hdulist[1:], key=lambda hdu: hdu.name)
# Rebuild HDUList with sorted HDUs
sorted_hdulist = fits.HDUList([primary_hdu] + sorted_bintable_hdus)
print('saving to:',filename_out)
sorted_hdulist.writeto(filename_out, overwrite=True)
