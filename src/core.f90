!*****************************************************************************************
!>
!  Main module.
!
!  Adapted from
!   * https://ccmc.gsfc.nasa.gov/pub/modelweb/geomagnetic/igrf/fortran_code/bilcal.for
!   * https://ccmc.gsfc.nasa.gov/pub/modelweb/radiation_belt/radbelt/fortran_code/radbelt.for

module core

   use radbelt_kinds_module
   use trmfun_module
   use shellig_module

   implicit none

   public :: igrf
   public :: get_flux

   contains

   !*****************************************************************************************
   !>
   !  Main routine.
   !
   !@todo we need to read in the coefficients only once and keep them in memory,
   !      rather than everytime these functions are called !

   function get_flux(Lon,Lat,Height,Year,E,Imname) result(flux)

      real(wp) :: lon, lat, height, year, e
      integer :: imname

      real(wp) :: flux,xl,bbx
      type(trm_type) :: trm

      call igrf(Lon,Lat,Height,Year,Xl,Bbx)
      call trm%aep8(e,Xl,Bbx,Imname,flux)

   end function get_flux

   !*****************************************************************************************
   !>
   ! Wrapper for IGRF functions.

   subroutine igrf(lon,lat,height,year,xl,bbx)

      real(wp) :: bab1 , babs , bdel , bdown , beast , beq , bequ , bnorth , dimo , rr0
      integer :: icode
      real(wp) :: lon , lat , height , year , xl , bbx
      logical :: val

      CALL feldcof(Year,dimo)
      CALL feldg(Lat,Lon,Height,bnorth,beast,bdown,babs)
      CALL shellg(Lat,Lon,Height,dimo,Xl,icode,bab1)
      bequ = dimo/(Xl*Xl*Xl)
      IF ( icode==1 ) THEN
         bdel = 1.0e-3_wp
         CALL findb0(0.05_wp,bdel,val,beq,rr0)
         IF ( val ) bequ = beq
      ENDIF
      Bbx = babs/bequ

   end subroutine igrf

end module core
