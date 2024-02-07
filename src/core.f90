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
!  Main routine.

   function get_flux(Lon,Lat,Height,Year,E,Imname) result(flux)

      real(wp),intent(in) :: lon
      real(wp),intent(in) :: lat
      real(wp),intent(in) :: height
      real(wp),intent(in) :: year
      real(wp),intent(in) :: e
      integer,intent(in) :: imname
      real(wp) :: flux

      real(wp) :: xl, bbx
      type(trm_type) :: trm
      type(shellig_type) :: igrf

      call igrf%igrf(lon,lat,height,year,xl,bbx)
      call trm%aep8(e,xl,bbx,imname,flux)

   end function get_flux
!*****************************************************************************************

end module core
