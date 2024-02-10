Radbelt: Work in progress to refactor the AE-8/AP-8 Van Allen belt model.

### Status

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![GitHub release](https://img.shields.io/github/release/jacobwilliams/radbelt.svg)](https://github.com/jacobwilliams/radbelt/releases/latest)
[![CI Status](https://github.com/jacobwilliams/radbelt/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/radbelt/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/radbelt/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/radbelt)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/radbelt)](https://github.com/jacobwilliams/radbelt/commits/master)

### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library and test cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `radbelt` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
radbelt = { git="https://github.com/jacobwilliams/radbelt.git" }
```

### Documentation

The latest API documentation can be found [here](https://jacobwilliams.github.io/radbelt/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### Original source

* The original sourcecode was hosted at GSFC "Modelweb", which no longer exists, but an archive can be found at the [Internet Archive](https://web.archive.org/web/20210318113325/https://ccmc.gsfc.nasa.gov/pub/modelweb/). It is presumed to be in the public domain. Reference: National Space Science Data Center, Data set PT-11B, Mar 1996. Dieter Bilitza, GSFC/NSSDC code 633, Greenbelt, MD 20771, tel. (301) 286-0190, dbilitza@pop600.gsfc.nasa.gov

### See also

* [NASA ModelWebArchive](https://git.smce.nasa.gov/ccmc-share/modelwebarchive)
* [An Astropy-friendly wrapper for the AE-8/AP-8 Van Allen belt model](https://github.com/nasa/radbelt)
* [pyIGRF](https://github.com/rilma/pyIGRF)
* https://github.com/lanl/RAM-SCB/blob/master/srcExternal/igrf.f
* https://github.com/space-physics/igrf/blob/main/src/igrf/fortran/igrf13.f
* https://web.archive.org/web/20210318113325/https://ccmc.gsfc.nasa.gov/pub/modelweb/
* Model parameters can be computed and plotted online at http://nssdc.gsfc.nasa.gov/space/model/ [broken link]

### Test case

See the `radbelt_test.f90` and `test.py` files:

Code | Runtime (sec) | Cases per second
--- | --- | ---
[Python](https://github.com/nasa/radbelt) version                  | 3.514   |    409
Fortran Function version (`get_flux()`) | 1.622   |   1198
Fortran Class version (`radbelt_type%get_flux()`)   | 0.017   | 112259

The main difference in speed from using the class method is that the data files are only read once, rather than each time the function is called (which is done in the other two versions).

### Brief description

These empirical models describe the differential or
integral, omnidirectional fluxes of electrons (AE-8) and protons
(AP-8) in the inner and outer radiation belts (electrons: L=1.1
to 11, protons: L=1.1 to 7) for two epochs representing solar
maximum (1970) and minimum (1964) conditions. The energy spectrum
ranges from 0.1 to 400 MeV for the protons and from 0.04 to 7 MeV
for the electrons. AE-8 and AP-8 are the most recent ones in a
series of models established by J. Vette and his colleges at NSSDC
starting in the early sixties. The models are based on almost all
available satellite data. It is IMPORTANT that the models maps for
solar maximum are used with a magnetic field model for epoch=1970
and for solar minimum for epoch=1964.

For each epoch and particle the model consists of a three-
dimensional table of (logarithm of) particle fluxes in energy, L-value,
and B/B0 (magnetic field strength normalized to the equator). The program
MODEL finds the particle fluxes for given energy, L-value and B/B0 by
interpolating in energy (subroutine TRARA1) and in L * B/B0 space (TRARA2).
The program RADBELT produces tables of integral or differential fluxes
for different energies varying with L or B/B0.

The coefficient files are provided in ASCII (*.asc) format:

Description | Filename | Size (KB)
--- | --- | ---
AE-8, epoch 1970, solar maximum       |  	ae8max.asc | 84
AE-8, epoch 1964, solar minimum	 	   |    ae8min.asc | 81
AP-8, epoch 1970, solar maximum		   |    ap8max.asc | 101
AP-8, epoch 1964, solar minimum  	 	 |    ap8min.asc | 102

In March 1995 the earlier used compressed model maps AP8MIC and AP8MAC
were replaced with the full maps AP8MIN/MAX with the help of D. Heynderickx
(BIRA, Brussel, Belgium) and A. Beliaev (INP/MSU, Moscow, Russia). Heynderickx
and Beliaev (1995) had found and corrected a small error in the AP8MIN map;
two lines had been exchanged.


### References

* G.W. Singley, and J.I. Vette, The AE-4 Model of the Outer Radiation
  Zone Electron Environment, NSSDC/WDC-A-R&S 72-06, 1972.
* M.J. Teague, and J.I. Vette, A Model of the Trapped Electron
  Population for Solar Minimum (AE-5), NSSDC/WDC-A-R&S 74-03, 1974.
* M.J. Teague, K.W. Chan, and J.I. Vette, AE-6: A Model Environment
  of Trapped Electrons for Solar Maximum, NSSDC/WDC-A-R&S 76-04, 1976
* D.W. Sawyer, and J.I. Vette, AP-8 Trapped Proton Environment for
  Solar Maximum and Minimum, NSSDC/WDC-A-R&S 76-06, 1976.
* J.I. Vette, K.W. Chan, and M.J. Teague, Problems in Modeling the
  Earth's Trapped Radiation Environment, AFGL-TR-78-0130, 1978.
* K.W. Chan, M.J. Teague, N.J. Schofield, and J.I. Vette, Modeling of
  Electron Time Variation in the Radiation Belts, p. 121-149, in:
  Quantitative Modeling of Magnetospheric Processes, W.P. Olson
  (ed.), geophysical monograph 21, American Geophysical Union, 1979.
* M.T. Teague, N.J. Schofield, K.W. Chan, and J.I. Vette, A Study of
  Inner Zone Electron Data and their Comparison with Trapped
  Radiation Models, NSSDC/WDC-A-R&S 79-06, 1979.
* J.I. Vette, The AE-8 Trapped Electron Model Environment,
  NSSDC/WDC-A-R&S 91-24, 1991.
* J.I. Vette, The NASA/National Space Science Data Center Trapped
  Radiation Environment Model Program (1964-1991), NSSDC/WDC-A-R&S
  91-29, 1991.
* D. Heynderickx and A. Beliaev, J. Spacecraft and Rockets 32, 190-192, 1995.
