#!/usr/bin/env python3

import glob
import os

import astropy.units as u
import numpy as np
import spiceypy


def compute_spice_pointing(timestamps):
    ''' Compute SPICE pointing

    Parameters
    ==========
    timestamps : list of str, size (nt,)
        List of timestamps at which to get the pointing. These timestamps are
        strings describing an epoch, parsable by spiceypy.str2et (e.g.
        YYYY-MM-DDTHH:MM:SS.SSS)..

    Returns
    =======
    Tx : array of shape (nt,)
        Helioprojective westward angle.
    Ty : array of shape (nt,)
        Helioprojective northward angle.
    roll : array of shape (nt,)
        Spacecraft roll angle.
    '''

    # Load SPICE kernels
    spiceypy.kclear()
    kernels_folder = os.environ['SPICE_KERNELS_SOLO']
    mk_file = sorted(glob.glob(os.path.join(kernels_folder, 'mk', 'solo_ANC_soc-flown-mk_*.tm')))[-1]
    spiceypy.furnsh(mk_file)

    # Compute pointing
    Tx = []
    Ty = []
    roll = []
    for et in spiceypy.str2et(timestamps):
        rot_mat = spiceypy.pxform('SOLO_SUN_RTN', 'SOLO_SRF', et)
        r1, r2, r3 = spiceypy.m2eul(rot_mat, 1, 2, 3)
        Tx.append(r3 + np.pi)
        Ty.append(- r2)
        roll.append(- r1)
    Tx = u.Quantity(Tx, 'rad').to('arcsec')
    Ty = u.Quantity(Ty, 'rad').to('arcsec')
    roll = u.Quantity(roll, 'rad').to('deg')

    return Tx, Ty, roll
