#
# Python interface to the `radbelt_fortran` shared library.
#

from . import radbelt_fortran

class RadbeltClass:
    """
    Class for using the radbelt model.
    """

    def __init__(self, aep8_dir : str = None, igrf_dir: str = None) -> None:
        """
        Constructor for the fortran class

        Parameters
        ----------
        aep8_dir : str
            The directory containing the aep8 files. If None, then use the default directory ('data/aep8/')
        igrf_dir : str
            The directory containing the igrf files. If None, then use the default directory ('data/igrf/')
        """

        #`ip` is an integer that represents a c pointer
        # to a `radbelt_type` in the Fortran library.
        self.ip = radbelt_fortran.radbelt_c_module.initialize_c()

        # note that None means use defaults,
        # but '' means current directory
        if aep8_dir is not None:
            self.set_trm_file_path(aep8_dir)
        if igrf_dir is not None:
            self.set_igrf_file_path(igrf_dir)

    def __del__(self) -> None:
        """destructor for the fortran class"""

        radbelt_fortran.radbelt_c_module.destroy_c(self.ip)

    def set_trm_file_path(self, aep8_dir : str) -> None:
        """Set just the aep8 file path"""

        radbelt_fortran.radbelt_c_module.set_trm_file_path_c(self.ip, aep8_dir, len(aep8_dir))

    def set_igrf_file_path(self, igrf_dir : str) -> None:
        """Set just the igrf file path"""

        radbelt_fortran.radbelt_c_module.set_igrf_file_path_c(self.ip, igrf_dir, len(igrf_dir))

    def set_data_files_paths(self, aep8_dir : str, igrf_dir : str) -> None:
        """Set both the file paths"""

        radbelt_fortran.radbelt_c_module.set_data_files_paths_c(self.ip, aep8_dir, igrf_dir,
                                                                len(aep8_dir), len(igrf_dir))

    def get_flux(self, lon : float, lat : float , height : float, year : float,
                 e : float, particle : str, solar : str) -> float:
        """
        Compute the flux.

        Parameters
        ----------
        lon : float
            Longitude (deg)
        lat : float
            Latitude (deg)
        height : float
            Altitude (deg)
        year : float
            Decimal year (needed to account for drift of the Earth's magnetic field).
        energy : float
            Minimum energy.
        particle : {'e', 'p'}
            Particle species: 'e' for electrons, 'p' for protons.
        solar : {'min', 'max'}
            Solar activity: solar minimum or solar maximum.

        Returns
        -------
        flux : float
            The flux of particles above the given energy, (cm^-2 s^-1).
        """

        model = _model_index(particle, solar)
        return radbelt_fortran.radbelt_c_module.get_flux_g_c(self.ip,lon,lat,height,year,e,model)

def _model_index(particle : str, solar : str) -> int:
    """convert the particle and solar options to model number for the low-level routine"""

    p = particle.lower()
    s = solar.lower()
    if p not in ('e', 'p'): raise ValueError('particle must be "e" or "p"')
    if s not in ('min', 'max'): raise ValueError('solar must be "min" or "max"')
    return { ('e', 'min'): 1,
             ('e', 'max'): 2,
             ('p', 'min'): 3,
             ('p', 'max'): 4 }[(p, s)]