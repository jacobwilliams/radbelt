#
# Test of the f2py-created Python interface.
#

import sys
import time
from pathlib import Path

dir = Path(__file__).resolve().parents[1] # root directory
sys.path.insert(0,str(dir / 'python'))   # assuming the radbelt lib is in python directory
import radbelt  # this is the module being tested

class radbelt_class:
    """
    Class for using the radbelt model.
    """

    def __init__(self) -> None:
        """constructor for the fortran class"""

        #`ip` is an integer that represents a c pointer
        # to a `radbelt_type` in the Fortran library.
        self.ip = radbelt.radbelt_c_module.initialize_c()

    def __del__(self) -> None:
        """destructor for the fortran class"""

        radbelt.radbelt_c_module.destroy_c(self.ip)

    def set_data_files_paths(self, aep8_dir : str, igrf_dir : str) -> None:
        """Set the file paths"""

        radbelt.radbelt_c_module.set_data_files_paths_c(self.ip, aep8_dir, igrf_dir, len(aep8_dir), len(igrf_dir))

    def get_flux(self, lon : float, lat : float , height : float, year : float, e : float, imname : int) -> float:
        """Get the flux"""

        return radbelt.radbelt_c_module.get_flux_g_c(self.ip, lon,lat,height,year,e,imname)

#########################################################################################################
def set_data_files_paths(aep8_dir : str, igrf_dir : str) -> None:
    """Python function to set the file paths"""
    radbelt.radbelt_c_module.set_data_files_paths_c(aep8_dir, igrf_dir)

#########################################################################################################
#@np.vectorize   <-- makes it run slower for scalar inputs??
def get_flux(lon : float, lat : float , height : float, year : float, e : float, imname :int) -> float:
    """Python function to get the flux"""
    return radbelt.radbelt_c_module.get_flux_g_c(lon,lat,height,year,e,imname)

model = radbelt_class()

# set location of the data files:
aep8_dir = str(dir / 'data' / 'aep8')
igrf_dir = str(dir / 'data' / 'igrf')
model.set_data_files_paths(aep8_dir, igrf_dir)

EPS = sys.float_info.epsilon # machine precision for error checking
lon = -45.0
lat = -30.0
height = 500.0
year = 2021.1616438356164  # decimal year
imname = 4 # 'p', 'max'
e = 20.0

flux = model.get_flux(lon,lat,height,year,e,imname)

print(f'flux = {flux}')

error = flux - 2642.50370051985726336559603128948869  #difference from real128 version
relerror = abs(error/flux)

print(f'Flux      = {flux}')
print(f'Error     = {error}')
print(f'Rel Error = {relerror}')
if relerror>10*EPS:
    raise Exception('error')

print('')

n_cases = 0
tstart = time.time()
for lat in range(-90, 91, 5):
    for lon in range(-180, 181, 45):
        for alt in range(500, 1001, 100):
            n_cases = n_cases + 1
            f = model.get_flux(lon,lat,height,year, e, 4)
tend = time.time()
print(f'Python version runtime: {tend-tstart} sec. {int(n_cases/(tend-tstart))} (cases/sec)')

# print('... VECTORIZE ...')
# n_cases = 0
# tstart = time.time()

# lon_vec = []
# lat_vec = []
# height_vec = []
# year_vec = []
# e_vec = []
# ifile_vec = []

# for lat in range(-90, 91, 5):
#     for lon in range(-180, 181, 45):
#         for alt in range(500, 1001, 100):
#             n_cases = n_cases + 1
#             lon_vec.append(lon)
#             lat_vec.append(lat)
#             height_vec.append(height)
#             year_vec.append(year)
#             e_vec.append(e)
#             ifile_vec.append(4)
#             #f = get_flux(lon,lat,height,year, e, 4)

# f_vec = get_flux(lon_vec,lat_vec,height_vec,year_vec,e_vec,ifile_vec)
# tend = time.time()
# print(f'Python version runtime: {tend-tstart} sec. {int(n_cases/(tend-tstart))} (cases/sec)')

