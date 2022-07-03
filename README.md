# KMAG - Khurana Magnetic Field Model (Saturn, 2012)

## Introduction

KMAG is a global magnetic field of the Saturn’s field (Khurana et al., 2006) that includes modules that specify:
    (1) the internal spherical harmonic field,
    (2) the ring current and the magnetotail current system,
    (3) the field from the radial current system which reinforces corotation on the outflowing plasma,
    (4) the shielding fields from an axially symmetric but hinged magnetopause and (5) the interconnection magnetic field from the solar wind IMF.

The current sheet field is based on models of disk-shaped current sheets by Tsyganenko and Peredo (1994). The tilt, and hinging of the current sheet is based on the general deformation technique of Tsyganenko (1998, 2002a, 2002b). The shielding fields are derived using a model of the magnetopause by Arridge et al. (2006) constructed from magnetopause crossings observed by both Cassini and Voyager. Saturn’s internal field uses magnetic moments derived from the Cassini Saturn Orbital Insertion (SOI) measurements (Dougherty, 2005). KMAG is based on a previous Euler-potential model of Jupter's magnetospheric field (Khurana, 1997).

This python package is a wrapper for the Fortran KMAG code included in the `src` directory, compiled using the `f2py` package.

## Installation

### Requirements
 Make sure your system has Fortran installed to successfully compile the source code:

 * Linux: `apt install gfortran`
 * Mac: `brew install gcc`

To install the package, simply run the following command in the main directory of the package:

```bash
pip install .
```

## Usage
The main function that computes the magnetic field strength is
called kmag(). It requires the following parameters to be passed through:

```python
rlt,brm,btm,bpm = kmag(time,epoch,r,theta,phi,by_imf,bz_imf,dp)
```
where
 * time (real) - a double precision variable denoting number of seconds from an epoch
 * epoch (string) - a five letter character which can be either 'ctime',
   'etime' or 'J2000'/upper or lower case
 * r, theta, phi (real) -  position of the field in a right-handed System III coordinate system (THETA and PHI are in radians)
 * by_imf, bz_imf (real) - Y and Z component of the interplanetary magnetic field in right-handed KSM coordinate system
 * rlt (real) - local time
 * brm, btm, bpm (real) - model magnetic field in a spherical SYSTEM III 
   coordinate system (BRM, BTM and BPM are in nT)

The supported coordinate systems are:
* S3C System III Cartesian (right-handed)
* KSO Kronian-Sun-Orbital
* KSM Kronian-Sun-Magnetic (So far same as KSO)
* DIP Dipole (cartesian)

### Example

See the example code given in `test/test.py` for an example of how to use the model. Run it with a shell command:

```bash
python test/test.py
```

The code will check if the model is working correctly and print the magnetic field for a selected reference position.

## References:
* Arridge, C. S., Achilleos, N., Dougherty, M. K., Khurana, K. K., & Russell, C. T. (2006). Modeling the size and shape of Saturn’s magnetopause with variable dynamic pressure. Journal of Geophysical Research, 111(A11), A11227. https://doi.org/10.1029/2005JA011574

* Dougherty, M. K. (2005). Cassini magnetometer observations during Saturn orbit insertion. Science, 307(5713), 1266–1270. https://doi.org/10.1126/science.1106098

* Khurana, K. K., Arridge, C. S., Schwarzl, H., & Dougherty, M. K. (2006). A model of Saturn’s magnetospheric field based on latest Cassini observations. In AGU Spring Meeting Abstracts (Vol. 2007, pp. P44A-01).

* Khurana, Krishan K. (1997). Euler potential models of Jupiter’s magnetospheric field. Journal of Geophysical Research: Space Physics, 102(A6), 11295–11306. https://doi.org/10.1029/97JA00563

* Tsyganenko, N. A. (1998). Modeling of twisted/warped magnetospheric configurations using the general deformation method. Journal of Geophysical Research: Space Physics, 103(A10), 23551–23563. https://doi.org/10.1029/98JA02292

* Tsyganenko, N. A. (2002a). A model of the near magnetosphere with a dawn-dusk asymmetry 1. Mathematical structure. Journal of Geophysical Research: Space Physics, 107(A8), SMP 12-1-SMP 12-15. https://doi.org/10.1029/2001JA000219

* Tsyganenko, N. A. (2002b). A model of the near magnetosphere with a dawn-dusk asymmetry 2. Parameterization and fitting to observations. Journal of Geophysical Research: Space Physics, 107(A8), SMP 10-1-SMP 10-17. https://doi.org/10.1029/2001JA000220

* Tsyganenko, N. A., & Peredo, M. (1994). Analytical models of the magnetic field of disk-shaped current sheets. Journal of Geophysical Research, 99(A1), 199. https://doi.org/10.1029/93JA02768

## License
[MIT](https://choosealicense.com/licenses/mit/)
