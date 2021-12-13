#!/usr/bin/env python3

import argparse
import itertools
import os

from astropy import wcs
from astropy.io import fits
from tqdm import tqdm
import astropy.units as u
import numpy as np
import scipy.interpolate as si

from compute_spice_pointing_with_spice import SpiceSpicePointing


def remap_spice_hdu(hdu, solo_Tx, solo_Ty, solo_roll):
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
        slit positions in the FITS.

    Returns
    =======
    hdu : astropy.io.fits.ImageHDU
        Aligned SPICE 'L2a' HDU
    '''
    if not hdu.is_image:
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
    for it, iD in tqdm(itD, desc=f'Remapping {hdu.name}', total=nt*nD):
        img = hdu.data[it, iD]
        interp = si.LinearNDInterpolator(points, img.flatten())
        new_img = interp(xi_interp)
        new_hdu.data[it, iD] = new_img
    return new_hdu


class PlotResults():
    def plot_pointing(self, timestamps, Tx, Ty, roll, filename):
        t = [np.datetime64(t) for t in timestamps]
        plt.clf()
        plt.plot(t, Tx - Tx[0], label='$\\theta_x$')
        plt.plot(t, Ty - Ty[0], label='$\\theta_y$')
        plt.legend()
        plt.gcf().autofmt_xdate()
        plt.xlabel('Time')
        plt.ylabel('$\\Delta\\theta$ [arcsec]')
        plt.savefig(filename)

    def plot_hdu(self, hdu, ax):
        img = hdu.data
        if img.ndim > 2:
            img = np.nansum(img, axis=1)  # Sum over wavelengths
            img = np.squeeze(img)  # Collapse 1-depth axis (t or X)
        img = img - np.nanmin(img)
        norm = mpl.colors.LogNorm(
            vmin=np.max([1, np.nanpercentile(img, 1)]),
            vmax=np.nanpercentile(img, 99),
            )
        im = ax.imshow(img, origin='lower', norm=norm, aspect=1/4)
        plt.title(hdu.name)
        plt.colorbar(im)

    def plot_hdulist(self, hdulist, filename):
        with PdfPages(filename) as pdf:
            for hdu in hdulist:
                if hdu.is_image:
                    plt.clf()
                    self.plot_hdu(hdu, plt.gca())
                    pdf.savefig()


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('file', nargs='+',
                   help=('SPICE L2 FITS to align.'))
    p.add_argument('-O', '--output-dir',
                   help=('output directory, if different '
                         'from that of the input files'))
    p.add_argument('--plot-results', action='store_true',
                   help='plot results')
    args = p.parse_args()
    if args.output_dir is not None:
        os.makedirs(args.output_dir, exist_ok=True)
    else:
        raise NotImplementedError  # TODO: handle this case

    spice_spice_pointing = SpiceSpicePointing()

    for filename in args.file:
        basename = os.path.splitext(os.path.basename(filename))[0]
        hdulist = fits.open(filename)
        timestamps = hdulist[-1].data['TIMAQUTC'][0, 0, 0, 0]
        Tx, Ty, roll = spice_spice_pointing.compute_pointing(timestamps)
        new_hdulist = fits.HDUList(hdus=[])
        for hdu in hdulist:
            new_hdulist.append(remap_spice_hdu(hdu, Tx, Ty, roll))
        new_hdulist.writeto(
            f'{args.output_dir}/{basename}_remapped.fits',
            overwrite=True)

        if args.plot_results:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
            p = PlotResults()
            p.plot_pointing(timestamps, Tx, Ty, roll, f'{args.output_dir}/{basename}_plot_TxTy.pdf')
            p.plot_hdulist(hdulist, f'{args.output_dir}/{basename}_original.pdf')
            p.plot_hdulist(new_hdulist, f'{args.output_dir}/{basename}_remapped.pdf')
