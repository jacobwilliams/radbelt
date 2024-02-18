#
# Test of the f2py-created Python interface.
#

import sys
import time
from pathlib import Path

dir = Path(__file__).resolve().parents[1] # root directory
sys.path.insert(0,str(dir / 'build-python'))   # assuming the radbelt lib is in python directory
from radbeltpy import RadbeltClass  # this is the module being tested

# location of the data files:
aep8_dir = str(dir / 'data' / 'aep8')
igrf_dir = str(dir / 'data' / 'igrf')

# create the class:
model = RadbeltClass(aep8_dir = aep8_dir, igrf_dir = igrf_dir)

EPS = sys.float_info.epsilon # machine precision for error checking
lon = -45.0
lat = -30.0
height = 500.0
year = 2021.1616438356164  # decimal year
particle = 'p'
solar = 'max'
e = 20.0

flux = model.get_flux(lon,lat,height,year,e,particle,solar)

print(f'flux = {flux}')

error = flux - 2642.50370051985726336559603128948869  #difference from real128 version
relerror = abs(error/flux)

print(f'Flux      = {flux}')
print(f'Error     = {error}')
print(f'Rel Error = {relerror}')
if relerror > 10.0*EPS:
    raise Exception('error')

print('')

n_cases = 0
tstart = time.time()
for lat in range(-90, 91, 5):
    for lon in range(-180, 181, 45):
        for alt in range(500, 1001, 100):
            n_cases = n_cases + 1
            f = model.get_flux(lon,lat,height,year, e, particle,solar)
tend = time.time()
print(f'Python version runtime: {tend-tstart} sec. {int(n_cases/(tend-tstart))} (cases/sec)')
