!*****************************************************************************************
!>
!  Main module.
!
!### See also
!   * https://ccmc.gsfc.nasa.gov/pub/modelweb/geomagnetic/igrf/fortran_code/bilcal.for
!   * https://ccmc.gsfc.nasa.gov/pub/modelweb/radiation_belt/radbelt/fortran_code/radbelt.for

module core

   use radbelt_kinds_module
   use trmfun_module
   use shellig_module

   implicit none

   type,public :: radbelt_type
      !! the main class that can be used to get the flux.
      private
      type(trm_type) :: trm
      type(shellig_type) :: igrf
      contains
      private
      procedure,public :: get_flux => get_flux_
      procedure,public :: set_data_files_paths
   end type radbelt_type

   public :: get_flux !! simple function version for testing

   contains

!*****************************************************************************************
!>
!  Set the paths to the data files.
!  If not used or blank, the folder `data/aep8` and `data/igrf` in the
!  current working directory is assumed

   subroutine set_data_files_paths(me, aep8_dir, igrf_dir)

      class(radbelt_type),intent(inout) :: me
      character(len=*),intent(in) :: aep8_dir
      character(len=*),intent(in) :: igrf_dir

      call me%trm%set_data_file_dir(trim(aep8_dir))
      call me%igrf%set_data_file_dir(trim(igrf_dir))

   end subroutine set_data_files_paths
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calculate the flux of trapped particles at a specific location and time.

   function get_flux_(me,lon,lat,height,year,e,imname) result(flux)

      class(radbelt_type),intent(inout) :: me
      real(wp),intent(in) :: lon !! geodetic longitude in degrees (east)
      real(wp),intent(in) :: lat !! geodetic latitude in degrees (north)
      real(wp),intent(in) :: height !! altitude in km above sea level
      real(wp),intent(in) :: year !! decimal year for which geomagnetic field is to
                                  !! be calculated (e.g.:1995.5 for day 185 of 1995)
      real(wp),intent(in) :: e !! minimum energy
      integer,intent(in) :: imname !! which method to use:
                                   !!
                                   !! * 1 -- particle species: electrons, solar activity: min
                                   !! * 2 -- particle species: electrons, solar activity: max
                                   !! * 3 -- particle species: protons, solar activity: min
                                   !! * 4 -- particle species: protons, solar activity: max
      real(wp) :: flux !! The flux of particles above the given energy, in units of cm^-2 s^-1.

      real(wp) :: xl !! l value
      real(wp) :: bbx

      call me%igrf%igrf(lon,lat,height,year,xl,bbx)
      call me%trm%aep8(e,xl,bbx,imname,flux)

   end function get_flux_
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calculate the flux of trapped particles at a specific location and time.
!  This is just a function version of the class method from [[radbelt_type]].
!
!@note This routine is not efficient at all since it will reload all the
!      files every time it is called.

   function get_flux(lon,lat,height,year,e,imname) result(flux)

      real(wp),intent(in) :: lon !! geodetic longitude in degrees (east)
      real(wp),intent(in) :: lat !! geodetic latitude in degrees (north)
      real(wp),intent(in) :: height !! altitude in km above sea level
      real(wp),intent(in) :: year !! decimal year for which geomagnetic field is to
                                  !! be calculated (e.g.:1995.5 for day 185 of 1995)
      real(wp),intent(in) :: e !! minimum energy
      integer,intent(in) :: imname !! which method to use:
                                   !!
                                   !! * 1 -- particle species: electrons, solar activity: min
                                   !! * 2 -- particle species: electrons, solar activity: max
                                   !! * 3 -- particle species: protons, solar activity: min
                                   !! * 4 -- particle species: protons, solar activity: max
      real(wp) :: flux !! The flux of particles above the given energy, in units of cm^-2 s^-1.

      type(radbelt_type) :: radbelt

      flux = radbelt%get_flux(lon,lat,height,year,e,imname)

   end function get_flux
!*****************************************************************************************

end module core
