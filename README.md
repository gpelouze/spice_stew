# SPICE stew: Spice Pointing Inference and Correction Exploiting the Spice Toolkit to Eliminate Wobbles

Tools to correct the pointing of SolarOrbiter/SPICE using data from the SPICE
toolkit.

## Installation

### Prerequisite: Solar Orbiter SPICE Kernels

The [SPICE kernels for Solar Orbiter][solo-spice-kernels] are required in order
to compute the pointing of SPICE.

The easiest way to install them and keep them up to date is to clone [ESA's git
repository][solo-spice-kernels-git] (~1 GiB) to a location of your choice:

```
git clone https://repos.cosmos.esa.int/socci/scm/spice_kernels/solar-orbiter.git /path/to/spice_kernels/SOLAR-ORBITER/
```

and replace `PATH_VALUES = = ( '..' )` in all `mk/*.tm` files with
`PATH_VALUES = ( '/path/to/spice_kernels/SOLAR-ORBITER/kernels' )` (see full
instruction in the repository's [README][solo-spice-kernels-git]).

Then, set the `$SPICE_KERNELS_SOLO` environment variable to point to this path.
This can be done by adding the following line to your `~/.bashrc`:

```
export SPICE_KERNELS_SOLO=/path/to/spice_kernels/SOLAR-ORBITER/kernels
```

[solo-spice-kernels]: https://www.cosmos.esa.int/web/spice/solar-orbiter
[solo-spice-kernels-git]: https://repos.cosmos.esa.int/socci/projects/SPICE_KERNELS/repos/solar-orbiter


### SPICE stew

To install the SPICE stew package, run:

```
pip install spice_stew
```


## Usage

### From the command line

Use the `spice_stew` command in the terminal, eg:

```
spice_stew /archive/SOLAR-ORBITER/SPICE/fits/level2/2021/09/14/solo_L2_spice-n-ras_20210914T025031_V06_67109159-000.fits -O ~/spice_stew_test/ --plot-results
```

Run `spice_stew -h` for more options.


### In Python

Use the `spice_stew.correct_spice_pointing` function, eg:

```
import spice_stew
ssp = spice_stew.SpiceSpicePointing()
spice_stew.correct_spice_pointing(
    ssp,
    '/archive/SOLAR-ORBITER/SPICE/fits/level2/2021/09/14/solo_L2_spice-n-ras_20210914T025031_V03_67109159-000.fits',
    '~/spice_stew_test/',
    plot_results=True,
    )
```


## Troubleshooting

- `SPICE(NOSUCHFILE)The first file '../ck/[...].bc' specified by
  KERNELS_TO_LOAD in the file [...]/mk/[...]..tm could not be located.` is
  raised when the `PATH_VALUES` have not been updated in the meta-kernels
  during the [installation of the kernels][prerequisite-solar-orbiter-spice-kernels].

- `SpiceNOFRAMECONNECT: At epoch ..., there is insufficient information
  available to transform from reference frame -144991 (SOLO_SUN_RTN) to
  reference frame -144000 (SOLO_SRF).` is raised when the requested date is not
  in the as-flown SPICE kernels. This happens either because the kernels have
  not been updated, or because a date in the future is being requested. In
  the 1st case, you can update the kernels with `git pull`.


## Reference / License

This package is released under a MIT open source licence. See `LICENSE.txt`.
