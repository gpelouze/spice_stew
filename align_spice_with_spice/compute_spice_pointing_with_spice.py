#!/usr/bin/env python3

import glob
import os
import re

import spiceypy


class SolarOrbiterSpiceKernels():
    ''' Compute the position of Solo and other bodies using SPICE kernels

    Adapted from MISO/mission_planning/models/spice_kernels.py
    '''

    def __init__(self):
        self.kernels_folder = os.environ['SPICE_KERNELS_SOLO']
        # The mk file loads all the kernel files we need
        self.mk_filename = sorted(glob.glob(os.path.join(self.kernels_folder, 'mk', 'solo_ANC_soc-flown-mk_*.tm')))[-1]
        spiceypy.kclear()
        spiceypy.furnsh(self.mk_filename)

    def get_version(self):
        ''' Get SPICE kernel version '''
        mk_name = os.path.basename(self.mk_filename)
        m = re.match('^solo_ANC_soc-(flown)-mk_(.*).tm$', mk_name)
        if m:
            return ' '.join(m.groups())
        else:
            return mk_name

    def get_rotation_matrix(self, object_id, frame, et):
        ''' Get the rotation matrix between an object and the HCI frame

        Parameters
        ==========
        object_id : str
            Target body name. Can be 'EARTH', 'SOLO', ...
        frame : str
            Name of the frame to transform to.
        et : float
            Epoch at which to get the rotation matrix.

        Returns
        =======
        pform : ndarray of shape (3, 3)
            The rotation matrix.
        '''
        if object_id == 'SOLO':
            try:
                return spiceypy.pxform('SOLO_SRF', frame, et)
            except spiceypy.exceptions.SpiceNOFRAMECONNECT:
                return spiceypy.pxform('SOLO_SUN_RTN', frame, et)
        else:
            return spiceypy.pxform('IAU_' + object_id, frame, et)

    def get_positions(self, object_id, timestamps, obs_id='SUN'):
        ''' Get the position of an object in the HCI frame

        Parameters
        ==========
        object_id : str
            Target body name. Can be 'EARTH', 'SOLO', ...
        timestamps : list of str
            List of timestamps at which to get the orbit. These timestamps are
            strings describing an epoch, parsable by spiceypy.str2et (e.g.
            YYYY-MM-DDTHH:MM:SS.SSS)..
        obs_id : str (default: 'SUN')
            The reference observing body name.

        Returns
        =======
        positions : list of dict
            List of positions for the object in the HCI frame, in km.
            Each item of the list is a dict with keys 'x', 'y', and 'z'.
        '''
        positions = []
        for et in spiceypy.str2et(timestamps):
            pos, lightTimes = spiceypy.spkpos(object_id, et, 'HCI', 'NONE', obs_id)
            positions.append({'x': pos[0], 'y': pos[1], 'z': pos[2]})
        return positions

    def get_rotations(self, object_id, timestamps):
        ''' Get the rotation of an object in the HCI frame

        Parameters
        ==========
        object_id : str
            Target body name. Can be 'EARTH', 'SOLO', ...
        timestamps : list of str
            List of timestamps at which to get the orbit. These timestamps are
            strings describing an epoch, parsable by spiceypy.str2et (e.g.
            YYYY-MM-DDTHH:MM:SS.SSS)..

        Returns
        =======
        rotations : list of dict
            List of rotation vectors for the object in the HCI frame.
            Each item of the list is a dict with keys 'x', 'y', and 'z'.
        '''
        rotations = []
        axis = [1., 0., 0.]
        for et in spiceypy.str2et(timestamps):
            pform = self.get_rotation_matrix(object_id, 'HCI', et)
            rot = spiceypy.mxv(pform, axis)
            rotations.append({'x': rot[0], 'y': rot[1], 'z': rot[2]})
        return rotations


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

    spicek = SolarOrbiterSpiceKernels()
    pos = spicek.get_positions('SOLO', timestamps)
    rot = spicek.get_rotations('SOLO', timestamps)
