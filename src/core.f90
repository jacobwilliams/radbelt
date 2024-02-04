
!>
! Adapted from
!  * https://ccmc.gsfc.nasa.gov/pub/modelweb/geomagnetic/igrf/fortran_code/bilcal.for
!  * https://ccmc.gsfc.nasa.gov/pub/modelweb/radiation_belt/radbelt/fortran_code/radbelt.for

module core 

   use radbelt_module
   use shellig_module

   implicit none 

   contains
  
SUBROUTINE igrf(Lon,Lat,Height,Year,Xl,Bbx)
   IMPLICIT NONE
   REAL bab1 , babs , bdel , bdown , beast , beq , bequ , bnorth , dimo , rr0
   INTEGER icode
   REAL Lon , Lat , Height , Year , Xl , Bbx
   LOGICAL val
 
   CALL initize()
   CALL feldcof(Year,dimo)
   CALL feldg(Lat,Lon,Height,bnorth,beast,bdown,babs)
   CALL shellg(Lat,Lon,Height,dimo,Xl,icode,bab1)
   bequ = dimo/(Xl*Xl*Xl)
   IF ( icode==1 ) THEN
      bdel = 1.E-3
      CALL findb0(0.05,bdel,val,beq,rr0)
      IF ( val ) bequ = beq
   ENDIF
   Bbx = babs/bequ
END SUBROUTINE igrf
  
SUBROUTINE aep8(E,L,Bb0,Imname,Flux)
   IMPLICIT NONE
   REAL E(1) , ee , Flux(1)
   INTEGER i , ierr , ihead , Imname , iuaeap , nmap
   integer,dimension(:),allocatable :: map
   REAL L , Bb0
   DIMENSION ihead(8) , ee(1)
   CHARACTER*10 name , mname(4)
   DATA mname/'ae8min.asc' , 'ae8max.asc' , 'ap8min.asc' , 'ap8max.asc'/
 
   iuaeap = 15
   name = mname(Imname)
 
   OPEN (iuaeap,FILE=name,STATUS='OLD',IOSTAT=ierr,ERR=100,FORM='FORMATTED')
   READ (iuaeap,99001) ihead
   nmap = ihead(8)
   allocate(map(nmap))
   READ (iuaeap,99001) (map(i),i=1,nmap)
 
 100  CLOSE (iuaeap)
   IF ( ierr/=0 ) then
      error stop 'error reading '//trim(name)
   end if
 
   ee(1) = E(1)
   CALL trara1(ihead,map,L,Bb0,E,Flux,1)
   IF ( Flux(1)>0.0 ) Flux(1) = 10.**Flux(1)
99001 FORMAT (1X,12I6)
END SUBROUTINE aep8

end module core 
