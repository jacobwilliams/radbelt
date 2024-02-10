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