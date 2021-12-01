#!/usr/bin/env python3

import miso_spice_kernels


def compute_spice_pointing(spice_fits):
    ''' Compute SPICE pointing

    Parameters
    ==========
    spice_fits : astropy.io.fits.HDUList
        SPICE L2 FITS for which to compute the pointing

    Returns
    =======
    pointing : TODO: specify
    '''

    import eipy; eipy.sh()

    t_beg = spice_fits[0].header['DATE-BEG']
    t_end = spice_fits[0].header['DATE-END']
    spicek = miso_spice_kernels.SpiceKernels()
    orbit = spicek.getOrbit('SOLO', t_beg, t_end)
