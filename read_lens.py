
import numpy as np
from astropy.io import fits

name = "/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/fits11.28.fits"

def get_simulated_lens(n):
    with fits.open(name) as hdul:
        # image = FitsRawFormatterBase.readImage(hdul)
        m = n
        while True:
            try:
                image = hdul[2*m-1].data
                #image = np.log10(hdul[3*m-1].data)
                print('lens min max: ', np.min(hdul[2*m-1].data), np.max(hdul[2*m-1].data))
                masked_image = hdul[2*m].data
            except:
                m = m + 1
                print('m :', m)
                continue
            break

        
        return image
    hdul.close()