
import numpy as np
from astropy.io import fits

name = "/data/a.saricaoglu/lenstronomy/lenstronomy-data/mock_lenses/12.21/1623_simulation.fits"

def get_simulated_lens(n):
    with fits.open(name) as hdul:
        # image = FitsRawFormatterBase.readImage(hdul)
        m = n
        while True:
            try:
                image = hdul[m].data
                #image = np.log10(hdul[3*m-1].data)
                print('lens min max: ', np.min(hdul[m].data), np.max(hdul[m].data))
                noisy_image = hdul[m+1].data
            except:
                m = m + 1
                print('m :', m)
                continue
            break

        
        return image, noisy_image
    hdul.close()