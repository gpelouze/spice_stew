#!/usr/bin/env python3

import glob
import os

import astropy.units as u
import numpy as np
import spiceypy


class SpiceSpicePointing():
    def __init__(self, kernels_folder=None, clean_start=True):
        ''' Compute the pointing of the SPICE spectrometer using SPICE kernels

        Parameters
        ==========
        kernels_folder : str or None (default: None)
            Folder containing SOLO SPICE kernels. If None, use the environment
            variable $SPICE_KERNELS_SOLO.
        clean_start : bool (default: True)
            If True, clear existing SPICE kernels at startup.
        '''
        if clean_start:
            self.clear_kernels()
        if kernels_folder is None:
            kernels_folder = os.environ['SPICE_KERNELS_SOLO']
        self.load_kernels(kernels_folder)

    def load_kernels(self, kernels_folder):
        ''' Load SPICE kernels for SOLO

        Parameters
        ==========
        kernels_folder : str
            Folder containing SOLO SPICE kernels.
        '''
        mk_file = sorted(glob.glob(os.path.join(kernels_folder, 'mk', 'solo_ANC_soc-flown-mk_*.tm')))[-1]
        spiceypy.furnsh(mk_file)

    def clear_kernels(self):
        ''' Clear SPICE kernels '''
        spiceypy.kclear()

    def _parse_timestamps(self, timestamps):
        if np.issubdtype(timestamps.dtype, np.datetime64):
            timestamps = np.datetime_as_string(timestamps)
        return spiceypy.str2et(timestamps)

    def compute_pointing(self, timestamps):
        ''' Compute SPICE pointing

        Parameters
        ==========
        timestamps : list of str or np.datetime64, size (nt,)
            List of timestamps at which to compute the pointing. These
            timestamps are strings describing an epoch, parsable by
            spiceypy.str2et (e.g. YYYY-MM-DDTHH:MM:SS.SSS)..

        Returns
        =======
        Tx : array of shape (nt,)
            Helioprojective westward angle.
        Ty : array of shape (nt,)
            Helioprojective northward angle.
        roll : array of shape (nt,)
            Spacecraft roll angle.
        '''
        Tx = []
        Ty = []
        roll = []
        for et in self._parse_timestamps(timestamps):
            rot_mat = spiceypy.pxform('SOLO_SUN_RTN', 'SOLO_SRF', et)
            r2, r3, r1 = spiceypy.m2eul(rot_mat, 2, 3, 1)
            Ty.append(- r2 + np.pi)
            Tx.append(r3)
            roll.append(r1 + np.pi)
        Tx = u.Quantity(Tx, 'rad').to('arcsec')
        Ty = u.Quantity(Ty, 'rad').to('arcsec')
        roll = u.Quantity(roll, 'rad').to('deg')

        return Tx, Ty, roll
