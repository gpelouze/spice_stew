#!/usr/bin/env python3

import argparse
import os

from astropy.io import fits
from tqdm import tqdm

from compute_spice_pointing_with_spice import SpiceSpicePointing


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('file', nargs='+',
                   help=('SPICE L2 FITS to align.'))
    p.add_argument('-O', '--output-dir',
                   help=('output directory, if different '
                         'from that of the input files'))
    args = p.parse_args()
    if args.output_dir is not None:
        os.makedirs(args.output_dir, exist_ok=True)
    else:
        raise NotImplementedError  # TODO: handle this case

    spice_spice_pointing = SpiceSpicePointing()

    for filename in args.file:
        fits_file = fits.open(filename)
        timestamps = fits_file[-1].data['TIMAQUTC'][0, 0, 0, 0]
        Tx, Ty, roll = spice_spice_pointing.compute_pointing(timestamps)
