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

   public :: get_flux

   contains

!*****************************************************************************************
!>
!  Calculate the flux of trapped particles at a specific location and time.
!
!@note This routine is not efficient since it will reload all the
!      files every time it is called.

   function get_flux(Lon,Lat,Height,Year,E,Imname) result(flux)

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
      type(trm_type) :: trm
      type(shellig_type) :: igrf

      call igrf%igrf(lon,lat,height,year,xl,bbx)
      call trm%aep8(e,xl,bbx,imname,flux)

   end function get_flux
!*****************************************************************************************

end module core
