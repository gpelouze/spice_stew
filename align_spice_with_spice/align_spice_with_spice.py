#!/usr/bin/env python3

import argparse
import os

from astropy.io import fits
from tqdm import tqdm

from compute_spice_pointing_with_spice import compute_spice_pointing

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

    fits_files = [fits.open(f) for f in tqdm(args.file, desc='Opening files')]
    timestamps = [f[-1].data['TIMAQUTC'][0, 0, 0, 0] for f in fits_files]
    pointing = [compute_spice_pointing(ts) for ts in timestamps]
