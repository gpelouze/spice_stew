#!/usr/bin/env python3

import argparse
import copy
import glob
import itertools
import os
import re
import warnings

from astropy import wcs
from astropy.io import fits
from tqdm import tqdm
import astropy.units as u
import numpy as np
import scipy.interpolate as si
import spiceypy

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    warnings.warn('Could not import matplotlib, visualisation will not work')


__all__ = ['SpiceSpicePointing', 'correct_spice_pointing']


class SpiceFilename(dict):
    re_spice_L123_filename = re.compile(
        r'''
        solo
        _(?P<level>L[123])
        _spice
            (?P<concat>-concat)?
            -(?P<slit>[wn])
            -(?P<type>(ras|sit|exp))
            (?P<db>-db)?
            (?P<int>-int)?
        _(?P<time>\d{8}T\d{6})
        _(?P<version>V\d{2})
        _(?P<SPIOBSID>\d+)-(?P<RASTERNO>\d+)
        (?P<ext>\..*)
        ''',
        re.VERBOSE
        )

    def __init__(self, filename):
        """ Represent SPICE filename

        Parameters
        ----------
        filename : str
            eg: /a/b/solo_L2_spice-n-ras_20210914T025031_V03_67109159-000.fits
        """
        super().__init__(self._parse_filename(filename))
        self._original_filename = filename

        self._check_filename_reconstruction()

    def to_str(self, **kw):
        """ Generate filename string

        Parameters
        ----------
        **kw :
            If specified, these keywords are replaced in the generated filename
            If not specified, the default values are used. Supported
            keywords are:
            - level ('L1', 'L2', ...)
            - concat ('-concat' or '')
            - slit ('n' or 'w')
            - type ('ras', 'sit', 'exp')
            - db ('-db' or '')
            - int ('-int' or '')
            - time (YYYYMMDDTHHMMSS)
            - version ('V03', ...)
            - SPIOBSID ('67109159', ...)
            - RASTERNO ('000', ...)
            - ext ('.fits', '.fits.gz', ...)

        Returns
        -------
        filename : str
            New filename with interpolated keywords.
        """
        if kw:
            props = self.copy()
            props.update(**kw)
        else:
            props = self
        basename = ('solo_{level}_spice{concat}-{slit}-{type}{db}{int}'
                    '_{time}_{version}_{SPIOBSID}-{RASTERNO}{ext}')
        filename = os.path.join(
            props['path'],
            basename.format(**props),
            )
        return filename

    def _parse_filename(self, filename):
        path, basename = os.path.split(filename)
        m = self.re_spice_L123_filename.match(basename)
        if m is None:
            raise ValueError(f'could not parse SPICE filename: {basename}')
        d = m.groupdict()
        for k, v in d.items():
            if v is None:
                d[k] = ''
        d['path'] = path
        return d

    def _check_filename_reconstruction(self):
        filename = self.to_str()
        msg = (f"Filename reconstruction error: "
               f"{self._original_filename} -> {filename}")
        assert filename == self._original_filename, msg

    def __repr__(self):
        return f'<{type(self).__name__} {self}>'

    def __str__(self):
        return self.to_str()


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
        mk_files_pattern = os.path.join(
            kernels_folder,
            'mk',
            'solo_ANC_soc-flown-mk_*.tm')
        mk_files = sorted(glob.glob(mk_files_pattern))
        if len(mk_files) == 0:
            raise ValueError(f'No SPICE kernels found in {kernels_folder}')
        mk_file = mk_files[-1]
        spiceypy.furnsh(mk_file)

    def clear_kernels(self):
        ''' Clear SPICE kernels '''
        spiceypy.kclear()

    def _parse_timestamps(self, timestamps):
        if np.issubdtype(timestamps.dtype, np.datetime64):
            timestamps = np.datetime_as_string(timestamps)
        return spiceypy.str2et(timestamps)

    def _ang_to_pipi(self, a):
        ''' Convert any angle to the range ]-pi, pi]

        Parameters
        ==========
        a : float or astropy quantity
            Angle to convert. If it has no quantity, assumes it is in radians.

        Returns
        =======
        a : float or astropy quantity
            Converted angle.
        '''
        if hasattr(a, 'unit'):
            pi = u.Quantity(np.pi, 'rad').to(a.unit)
        else:
            pi = np.pi
        return - ((- a + pi) % (2*pi) - pi)

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
            rot_mat = spiceypy.pxform('SOLO_SOLAR_MHP', 'SOLO_SPICE_LW_OPT', et)
            (a, b, c), (d, e, f), (g, h, i) = rot_mat
            rot_mat = np.array([
                [-c, -a, -b],
                [-f, -d, -e],
                [i, g, h],
                ])
            r, b, a = spiceypy.m2eul(rot_mat, 1, 2, 3)
            Tx.append(-a)
            Ty.append(-b)
            roll.append(-r)
        Tx = u.Quantity(Tx, 'rad').to('arcsec')
        Ty = u.Quantity(Ty, 'rad').to('arcsec')
        roll = u.Quantity(roll, 'rad').to('deg')

        return Tx, Ty, roll


def remap_spice_hdu(hdu, solo_Tx, solo_Ty, solo_roll, sum_wvl=False):
    ''' Remap a SPICE spectral cube to corrected coordinates

    Parameters
    ==========
    hdu : astropy.io.fits HDU
        SPICE L2 FITS HDU to remap. (If the HDU is not of 'image' type, return
        it without modification.)
    solo_Tx, solo_Ty, solo_roll : 1D arrays
        Array containing the helioprojective longitude (Tx) and latitude (Ty)
        pointed by SOLO, as well as the spacecraft roll, computed from SPICE
        kernels. All arrays must be of shape (nt,), where nt is the number of
        slit positions (for rasters) or exposures (for sit and stares) in the
        FITS.
    sum_wvl : bool (default: False)
        If True, sum along wavelength axis to generate a quicklook image.

    Returns
    =======
    hdu : astropy.io.fits.ImageHDU
        Aligned SPICE 'L2a' HDU
    '''
    if (not hdu.is_image) or (hdu.name == 'WCSDVARR'):
        return hdu
    # Get wcs coordinates from fits
    w = wcs.WCS(hdu.header)
    iy, ix = np.indices(hdu.data.shape[2:])  # x and y indices
    iD = np.zeros_like(ix)
    it = np.zeros_like(ix)
    Tx, Ty, _, _ = w.pixel_to_world(ix, iy, iD, it)
    pi = u.Quantity(np.pi, 'rad')
    Tx = (Tx + pi) % (2*pi) - pi
    Ty = (Ty + pi) % (2*pi) - pi

    # Correct wcs coordinates using the spice data
    delta_Tx = solo_Tx - np.mean(solo_Tx)
    delta_Ty = solo_Ty - np.mean(solo_Ty)
    new_Tx = Tx - delta_Tx
    new_Ty = Ty - delta_Ty

    # Remap to new coordinates within time and/or wvl slices
    points = (
        Tx.to('arcsec').value.flatten(),
        Ty.to('arcsec').value.flatten(),
        )
    xi_interp = np.moveaxis(np.stack((
        new_Tx.to('arcsec').value,
        new_Ty.to('arcsec').value,
        )), 0, -1)
    nt, nD, _, _ = hdu.data.shape
    itD = itertools.product(range(nt), range(nD))

    new_hdu = hdu.copy()
    if sum_wvl:
        # Integrated intensity
        img = np.nansum(hdu.data, axis=1)  # Sum over wavelengths
        img = np.squeeze(img)  # Collapse 1-depth axis (t or X)
        interp = si.LinearNDInterpolator(points, img.flatten())
        new_img = interp(xi_interp)
        new_hdu.data = new_img.reshape(1, 1, *new_img.shape)
    else:
        # Full slices
        for it, iD in tqdm(itD, desc=f'Remapping {hdu.name}', total=nt*nD):
            img = hdu.data[it, iD]
            interp = si.LinearNDInterpolator(points, img.flatten())
            new_img = interp(xi_interp)
            new_hdu.data[it, iD] = new_img
    new_hdu.update_header()
    new_hdu.header.add_history('spice_stew.py')
    new_hdu.add_datasum()
    new_hdu.add_checksum()
    return new_hdu


def get_spice_timestamps(hdulist):
    ''' Get the acquisition timestamps for each slit position or exposure

    Parameters
    ==========
    hdulist : astropy.fits.io.HDUList
        SPICE L2 FITS HDU list.

    Returns
    =======
    timestamps : array of size (nt,)
        Timestamps for each slit position (rasters) or exposure (sit and
        stares)
    '''
    # extract timestamps from binary table HDU
    timestamps = hdulist['VARIABLE_KEYWORDS'].data['TIMAQUTC'][0, 0, 0, 0]
    timestamps = np.array([np.datetime64(t) for t in timestamps])
    # extract exposure time from primary HDU header
    t_exp = np.timedelta64(int(1e3*hdulist[0].header['XPOSURE']), 'ms')
    # add half-exposure-time to get the center of the exposure
    return timestamps + t_exp / 2


class PlotResults():
    def plot_pointing(self, timestamps, Tx, Ty, roll, spice_name, filename):
        t = [np.datetime64(t) for t in timestamps]
        t_ref = np.min(t).astype('datetime64[s]')
        t = (t - t_ref).astype('timedelta64[s]').astype(float) / 3600  # hour
        plt.clf()
        Tx_ref = np.mean(Tx).to('arcsec')
        Ty_ref = np.mean(Ty).to('arcsec')
        Tx = Tx - Tx_ref
        Ty = Ty - Ty_ref
        plt.plot(t, Tx, label=f'$X {Tx_ref.value:+.1f}$ arcsec')
        plt.plot(t, Ty, label=f'$Y {Ty_ref.value:+.1f}$ arcsec')
        plt.legend()
        plt.title(spice_name, fontsize=12)
        plt.xlabel(f'Time since {t_ref} [h]')
        plt.ylabel('Relative pointing angle [arcsec]')
        plt.savefig(filename)
        return t_ref

    def plot_hdu(self, hdu, ax, spice_name):
        img = hdu.data
        if img.ndim > 2:
            img = np.nansum(img, axis=1)  # Sum over wavelengths
            img = np.squeeze(img)  # Collapse 1-depth axis (t or X)
        def nanptp(a):
            return np.nanmax(a) - np.nanmin(a)
        img = img - np.nanmin(img) + 0.01*nanptp(img)
        norm = mpl.colors.LogNorm(
            vmin=np.nanpercentile(img, 1),
            vmax=np.nanpercentile(img, 99),
            )
        im = ax.imshow(img, origin='lower', norm=norm, aspect=1/4)
        plt.title(f'{spice_name}\n{hdu.name}', fontsize=12)
        plt.xlabel('X [px]')
        plt.ylabel('Y [px]')
        plt.colorbar(im)

    def plot_hdulist(self, hdulist, spice_name, filename):
        with PdfPages(filename) as pdf:
            for hdu in hdulist:
                if hdu.is_image and (hdu.name != 'WCSDVARR'):
                    plt.clf()
                    self.plot_hdu(hdu, plt.gca(), spice_name)
                    pdf.savefig()


def correct_spice_pointing(spice_spice_pointing, filename, output_dir,
                           overwrite=False, plot_results=False,
                           sum_wvl=False, windows=None):
    ''' Correct the pointing in a SPICE level 2 FITS

    Parameters
    ==========
    spice_spice_pointing : SpiceSpicePointing
        SPICE kernels interface for the SPICE spectrometer
    filename : str
        SPICE L2 FITS
    output_dir : str
        directory where the corrected FITS and plots are saved
    overwrite : bool (default: False)
        overwrite output FITS if it already exists
    plot_results : bool (default: False)
        generate plots to visualize the results
    sum_wvl : bool (default: False)
        If True, sum along wavelength axis to generate a quicklook image.
    windows : list of str or None (default: None)
        Windows to keep. If None, keep all windows.

    The aligned fits are saved in <output_dir> under the name
    <solo_L2_spice_..._remapped.fits> when sum_wvl is False, or
    <solo_L2_spice_..._remapped_img.fits>. when sum_wvl is True.

    Returns
    =======
    output_fits : str
        Path to the saved FITS
    '''
    os.makedirs(output_dir, exist_ok=True)
    # filename operations
    filename = SpiceFilename(filename)
    filename_output = copy.copy(filename)
    filename_output['path'] = output_dir
    if sum_wvl:
        filename_output['level'] = 'L2r_quicklook'
    else:
        filename_output['level'] = 'L2r'
    output_fits = filename_output.to_str()

    if os.path.isfile(output_fits) and not overwrite:
        print(f'Aligned file exists: {output_fits}, exiting')
        return output_fits

    # open FITS
    hdulist = fits.open(filename.to_str())
    timestamps = get_spice_timestamps(hdulist)

    # determine pointing
    Tx, Ty, roll = spice_spice_pointing.compute_pointing(timestamps)

    # interpolate data
    new_hdulist = fits.HDUList(hdus=[])
    if windows is None:
        windows = [hdu.name for hdu in hdulist
                   if hdu.is_image and (hdu.name != 'WCSDVARR')]
    windows = windows + [hdu.name for hdu in hdulist
                         if not hdu.is_image or (hdu.name == 'WCSDVARR')]
    for win in windows:
        hdu = hdulist[win]
        new_hdu = remap_spice_hdu(hdu, Tx, Ty, roll, sum_wvl=sum_wvl)
        new_hdulist.append(new_hdu)

    # save data
    new_hdulist.writeto(output_fits, overwrite=True)

    # generate plots
    if plot_results:
        p = PlotResults()
        spice_name = os.path.basename(hdulist.filename())
        p.plot_pointing(
            timestamps, Tx, Ty, roll, spice_name,
            filename_output.to_str(ext='_plot_TxTy.pdf'),
            )
        p.plot_hdulist(
            hdulist, spice_name,
            filename_output.to_str(ext='_plot_not_aligned.pdf'),
            )
        p.plot_hdulist(
            new_hdulist, spice_name,
            filename_output.to_str(ext='_plot_aligned.pdf'),
            )

    return output_fits


def main():
    p = argparse.ArgumentParser(
        description=('Correct the pointing of the SPICE spectrometer '
                     'using SPICE kernels'),
        )
    p.add_argument('files', nargs='+',
                   help='SPICE L2 FITS to align')
    p.add_argument('-O', '--output-dir', required=True,
                   help='output directory (required)')
    p.add_argument('-p', '--plot-results', action='store_true',
                   help='generate plots to visualize the results')
    p.add_argument('--sum-wvl', action='store_true',
                   help='save wavelength-integrated images')
    args = p.parse_args()

    spice_spice_pointing = SpiceSpicePointing()

    for filename in args.files:
        correct_spice_pointing(
            spice_spice_pointing,
            filename,
            args.output_dir,
            plot_results=args.plot_results,
            sum_wvl=args.sum_wvl,
            )


if __name__ == '__main__':
    main()
