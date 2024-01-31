National Space Science Data Center      Data set  PT-11B         Mar 1996 
=========================================================================

```
NSSDC-Id     PT-11B (combines PT-11A, MT-18A, MT-18B, MT-2AA, MT-2AB)

NAME: 	     Models of trapped particle fluxes AE-8 (electrons) and
	     AP-8 (protons) in inner and outer radiation belts

SOURCE:      Dieter Bilitza, GSFC/NSSDC code 633, Greenbelt,
             MD 20771, tel. (301) 286-0190
             dbilitza@pop600.gsfc.nasa.gov

CONTENT:     12 files					      *.*   blocks
	     	FORTRAN source code:  
	     driver program with interface          	RADBELT.FOR   48
             subroutines, functions 		 	 TRMFUN.FOR   30

		model tables binary and ASCII:
             AE-8, epoch 1964, solar minimum	 	 AE8MIN.BIN  104
	     AE-8, epoch 1970, solar maximum         	 AE8MAX.BIN  107
	     AP-8, epoch 1964, solar minimum  	 	 AP8MIN.BIN   53
	     AP-8, epoch 1970, solar maximum		 AP8MAX.BIN   52
             ASCII version of AE8MIN.BIN             	 AE8MIN.ASC  163
             ASCII version of AE8MAX.BIN             	 AE8MAX.ASC  168
             ASCII version of AP8MIN.BIN            	 AP8MIN.BIN   83
             ASCII version of AP8MAX.BIN            	 AP8MAX.BIN   81

	     user manual with examples  		RADBELT.LOG   48
             This file                                AAAREADME.DOC    9
```

### BRIEF DESCRIPTION:

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

The coefficient files are provided in VAX/VMS binary (*.bin) and
ASCII (*.asc) format. For all systems other than VMS the use of the ASCII
files is recommended. If using the ASCII coefficient one needs to slightly 
modify the RADBELT.FOR program as described in a comment statement in
RADBELT (this comments are found after the OPEN statement for the coefficient 
file).

In March 1995 the earlier used compressed model maps AP8MIC and AP8MAC
were replaced with the full maps AP8MIN/MAX with the help of D. Heynderickx 
(BIRA, Brussel, Belgium) and A. Beliaev (INP/MSU, Moscow, Russia). Heynderickx 
and Beliaev (1995) had found and corrected a small error in the AP8MIN map; 
two lines had been exchanged.

### AVAILABILITY:

(1) FORTRAN source code from this directory.
(2) Model parameters can be computed and plotted online at http://nssdc.gsfc.nasa.gov/space/model/ .


### See also

* [NASA ModelWebArchive](https://git.smce.nasa.gov/ccmc-share/modelwebarchive)

### REFERENCES:

G.W. Singley, and J.I. Vette, The AE-4 Model of the Outer Radiation
  Zone Electron Environment, NSSDC/WDC-A-R&S 72-06, 1972.

M.J. Teague, and J.I. Vette, A Model of the Trapped Electron
  Population for Solar Minimum (AE-5), NSSDC/WDC-A-R&S 74-03, 1974.

M.J. Teague, K.W. Chan, and J.I. Vette, AE-6: A Model Environment
  of Trapped Electrons for Solar Maximum, NSSDC/WDC-A-R&S 76-04, 1976

D.W. Sawyer, and J.I. Vette, AP-8 Trapped Proton Environment for
  Solar Maximum and Minimum, NSSDC/WDC-A-R&S 76-06, 1976.

J.I. Vette, K.W. Chan, and M.J. Teague, Problems in Modeling the
  Earth's Trapped Radiation Environment, AFGL-TR-78-0130, 1978.

K.W. Chan, M.J. Teague, N.J. Schofield, and J.I. Vette, Modeling of
  Electron Time Variation in the Radiation Belts, p. 121-149, in:
  Quantitative Modeling of Magnetospheric Processes, W.P. Olson
  (ed.), geophysical monograph 21, American Geophysical Union, 1979.

M.T. Teague, N.J. Schofield, K.W. Chan, and J.I. Vette, A Study of
  Inner Zone Electron Data and their Comparison with Trapped
  Radiation Models, NSSDC/WDC-A-R&S 79-06, 1979.

J.I. Vette, The AE-8 Trapped Electron Model Environment, 
  NSSDC/WDC-A-R&S 91-24, 1991.

J.I. Vette, The NASA/National Space Science Data Center Trapped 
  Radiation Environment Model Program (1964-1991), NSSDC/WDC-A-R&S 
  91-29, 1991.

(most of these references are available from NSSDC)

D. Heynderickx and A. Beliaev, J. Spacecraft and Rockets 32, 190-192, 1995.

National Space Science Data Center      Data set  PT-11B         Mar 1996 
=========================================================================
