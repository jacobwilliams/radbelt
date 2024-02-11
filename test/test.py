#
# Test of the Python version from: https://github.com/nasa/radbelt
#

from radbelt import get_flux
from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
import time

t = Time('2021-03-01')
energy = 20 * u.MeV

tstart = time.time()
n_cases = 0

# fortran loop indices:
# do ilat = -90, 90, 5
#     do ilon = -180, 180, 45
#         do ialt = 500, 1000, 100

for lat in range(-90, 91, 5):
    for lon in range(-180, 181, 45):
        for alt in range(500, 1001, 100):
            n_cases = n_cases + 1
            coords = EarthLocation(lon * u.deg, lat * u.deg, alt * u.km)
            f = get_flux(coords, t, energy, 'p', 'max')
            #print(t.utc.decimalyear, lat, lon, alt, f.value)
tend = time.time()
print(f'Python version runtime: {tend-tstart} sec. {int(n_cases/(tend-tstart))} (cases/sec)')

###############################################################
# read the results from the fortran test program and compare:
#
# this assumes the results.txt file has been generated by the radbelt_test test program.

print('')
print('python - fortran comparison:')

# example file:
# year, ilat, ilon, ialt, flux
# 2021.1616438356164      ,         -89 ,        -180 ,         500 ,   0.0000000000000000

n_cases = 0
n_failed_cases = 0
with open('results.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        n_cases = n_cases + 1
        line = line.strip().split(',')
        #year = line[0].strip()   # <----- assuming the same year as above
        lat  = float(line[1].strip())
        lon  = float(line[2].strip())
        alt  = float(line[3].strip())
        flux = float(line[4].strip())

        coords = EarthLocation(lon * u.deg, lat * u.deg, alt * u.km)
        f = get_flux(coords, t, energy, 'p', 'max')

        if flux == 0.0:
            error = abs(f.value-flux)
        else:
            error = abs(f.value-flux) / abs(flux)

        if (error > 1.0e-5):
            n_failed_cases = n_failed_cases + 1
            print(f'Error: {error}')

print(f'{n_failed_cases} failed cases out of {n_cases}')
if n_failed_cases>0:
    raise Exception('Test Failed')