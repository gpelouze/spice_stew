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

    timestamps = spice_fits[9].data['TIMAQUTC'][0, 0, 0, 0]

    spicek = miso_spice_kernels.SpiceKernels()
    orbit = spicek.getOrbit('SOLO', timestamps)
    pos, rot, ecl = spicek.getObjectPositions('SOLO', timestamps)
