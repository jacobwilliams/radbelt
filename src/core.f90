
!>
! Adapted from
!  * https://ccmc.gsfc.nasa.gov/pub/modelweb/geomagnetic/igrf/fortran_code/bilcal.for
!  * https://ccmc.gsfc.nasa.gov/pub/modelweb/radiation_belt/radbelt/fortran_code/radbelt.for

module core 

   use radbelt_kinds_module
   use radbelt_module
   use shellig_module

   implicit none 

   public :: igrf 
   public :: aep8
   public :: get_flux

   contains

   !TODO: we need to read in the coefficients only once and keep them in memory,
   !      rather than everytime these functions are called !

   function get_flux(Lon,Lat,Height,Year,E,Imname) result(flux)

      real(wp) :: lon, lat, height, year, e
      integer :: imname

      real(wp) :: flux,xl,bbx 
      real(wp), dimension(1) :: flux_, e_
      
      e_(1) = e
      
      call igrf(Lon,Lat,Height,Year,Xl,Bbx)
      call aep8(e_,Xl,Bbx,Imname,flux_)
      
      flux = flux_(1)
      
   end function get_flux
  
subroutine igrf(lon,lat,height,year,xl,bbx)

   real(wp) bab1 , babs , bdel , bdown , beast , beq , bequ , bnorth , dimo , rr0
   integer icode
   real(wp) lon , lat , height , year , xl , bbx
   logical val
 
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
END SUBROUTINE igrf
  
subroutine aep8(e,l,bb0,imname,flux)

   real(wp) e(1) , ee(1) , flux(1)
   integer i , ierr , ihead(8) , imname , iuaeap , nmap
   integer,dimension(:),allocatable :: map
   real(wp) l , bb0
   character(len=10) :: name 

   character(len=10),dimension(4),parameter :: mname = ['ae8min.asc' , &
                                                        'ae8max.asc' , &
                                                        'ap8min.asc' , &
                                                        'ap8max.asc']
 
   iuaeap = 15
   name = mname(Imname)
 
   OPEN (iuaeap,FILE=name,STATUS='OLD',IOSTAT=ierr,FORM='FORMATTED')
   IF ( ierr/=0 ) then
      error stop 'error reading '//trim(name)
   end if
   READ (iuaeap,'(1X,12I6)') ihead
   nmap = ihead(8)
   allocate(map(nmap))
   READ (iuaeap,'(1X,12I6)') (map(i),i=1,nmap)
   CLOSE (iuaeap)
 
   ee(1) = E(1)
   CALL trara1(ihead,map,L,Bb0,E,Flux,1)
   IF ( Flux(1)>0.0_wp ) Flux(1) = 10.0_wp**Flux(1)

END SUBROUTINE aep8

end module core 
