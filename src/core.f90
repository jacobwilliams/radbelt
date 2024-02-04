! Copyright (C) 2021 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. No copyright is claimed
! in the United States under Title 17, U.S. Code. All Other Rights Reserved.
!
! SPDX-License-Identifier: NASA-1.3

module core 

   use radbelt_module
   use shellig_module

   implicit none 

   contains
 
! Adapted from
! https://ccmc.gsfc.nasa.gov/pub/modelweb/geomagnetic/igrf/fortran_code/bilcal.for
 
SUBROUTINE igrf(Lon,Lat,Height,Year,Xl,Bbx)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL bab1 , babs , bdel , bdown , beast , beq , bequ , bnorth , dimo , rr0
   INTEGER icode
!*** End of declarations inserted by SPAG
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
 
! Adapted from
! https://ccmc.gsfc.nasa.gov/pub/modelweb/radiation_belt/radbelt/fortran_code/radbelt.for
 
SUBROUTINE aep8(E,L,Bb0,Imname,Flux)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL E(1) , ee , Flux(1)
   INTEGER i , ier , ierr , ihead , Imname , iuaeap , map , nmap
!*** End of declarations inserted by SPAG
   REAL L , Bb0
   DIMENSION map(20000) , ihead(8) , ee(1)     !JW WARNING: map doesn't have to be this big! should allocate it to size nmap !
   CHARACTER*10 name , mname(4)
   DATA mname/'ae8min.asc' , 'ae8max.asc' , 'ap8min.asc' , 'ap8max.asc'/
 
   iuaeap = 15
   name = mname(Imname)
 
   OPEN (iuaeap,FILE=name,STATUS='OLD',IOSTAT=ierr,ERR=100,FORM='FORMATTED')
   READ (iuaeap,99001) ihead
   nmap = ihead(8)
   READ (iuaeap,99001) (map(i),i=1,nmap)
 
 100  CLOSE (iuaeap)
   IF ( ier/=0 ) STOP
 
   ee(1) = E(1)
   CALL trara1(ihead,map,L,Bb0,E,Flux,1)
   IF ( Flux(1)>0.0 ) Flux(1) = 10.**Flux(1)
99001 FORMAT (1X,12I6)
END SUBROUTINE aep8

end module core 
