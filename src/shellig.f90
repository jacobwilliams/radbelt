
   module shellig_module 

      implicit none 

      contains 

! SHELLIG.FOR, Version 2.0, January 1992
!
! 11/01/91-DKB- SHELLG: lowest starting point for B0 search is 2
!  1/27/92-DKB- Adopted to IGRF-91 coeffcients model
!  2/05/92-DKB- Reduce variable-names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE
!  8/08/95-DKB- Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S
!  5/31/00-DKB- Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s
!  3/24/05-DKB- Updated to IGRF-45-10; new coeff.: IGRF05, IGRF05s
!

!
! SUBROUTINE findb0(Stps,Bdel,Value,Bequ,Rr0)
!    IMPLICIT NONE
! !*** Start of declarations inserted by SPAG
!    REAL b , Bdel , bdelta , Bequ , bmin , bold , bq1 , bq2 , bq3 , p , r1 , r2 , r3 , rold , Rr0 , Sp , step , step12 , Stps , zz
!    INTEGER i , irun , j , n
! !*** End of declarations inserted by SPAG
! !--------------------------------------------------------------------
! ! FINDS SMALLEST MAGNETIC FIELD STRENGTH ON FIELD LINE
! !
! ! INPUT:   STPS   STEP SIZE FOR FIELD LINE TRACING
! !  	COMMON/FIDB0/
! !	   SP     DIPOLE ORIENTED COORDINATES FORM SHELLG; P(1,*),
! ! 		  P(2,*), P(3,*) CLOSEST TO MAGNETIC EQUATOR
! !	   BDEL   REQUIRED ACCURACY  = [ B(LAST) - BEQU ] / BEQU
! !		  B(LAST)  IS FIELD STRENGTH BEFORE BEQU
! !
! ! OUTPUT:  VALUE  =.FALSE., IF BEQU IS NOT MINIMAL VALUE ON FIELD LINE
! !	   BEQU	  MAGNETIC FIELD STRENGTH AT MAGNETIC EQUATOR
! !	   RR0	  EQUATORIAL RADIUS NORMALIZED TO EARTH RADIUS
! !	   BDEL	  FINAL ACHIEVED ACCURACY
! !--------------------------------------------------------------------
!    DIMENSION p(8,4) , Sp(3)
!    LOGICAL Value
!    COMMON /fidb0 / Sp
!    INTEGER :: spag_nextblock_1
!    spag_nextblock_1 = 1
!    SPAG_DispatchLoop_1: DO
!       SELECT CASE (spag_nextblock_1)
!       CASE (1)
! !
!          step = Stps
!          irun = 0
!          spag_nextblock_1 = 2
!       CASE (2)
!          irun = irun + 1
!          IF ( irun>5 ) THEN
!             Value = .FALSE.
!             spag_nextblock_1 = 3
!             CYCLE SPAG_DispatchLoop_1
!          ENDIF
! !*********************FIRST THREE POINTS
!          p(1,2) = Sp(1)
!          p(2,2) = Sp(2)
!          p(3,2) = Sp(3)
!          step = -sign(step,p(3,2))
!          CALL stoer(p(1,2),bq2,r2)
!          p(1,3) = p(1,2) + 0.5*step*p(4,2)
!          p(2,3) = p(2,2) + 0.5*step*p(5,2)
!          p(3,3) = p(3,2) + 0.5*step
!          CALL stoer(p(1,3),bq3,r3)
!          p(1,1) = p(1,2) - step*(2.*p(4,2)-p(4,3))
!          p(2,1) = p(2,2) - step*(2.*p(5,2)-p(5,3))
!          p(3,1) = p(3,2) - step
!          CALL stoer(p(1,1),bq1,r1)
!          p(1,3) = p(1,2) + step*(20.*p(4,3)-3.*p(4,2)+p(4,1))/18.
!          p(2,3) = p(2,2) + step*(20.*p(5,3)-3.*p(5,2)+p(5,1))/18.
!          p(3,3) = p(3,2) + step
!          CALL stoer(p(1,3),bq3,r3)
! !******************INVERT SENSE IF REQUIRED
!          IF ( bq3>bq1 ) THEN
!             step = -step
!             r3 = r1
!             bq3 = bq1
!             DO i = 1 , 5
!                zz = p(i,1)
!                p(i,1) = p(i,3)
!                p(i,3) = zz
!             ENDDO
!          ENDIF
! !******************INITIALIZATION
!          step12 = step/12.
!          Value = .TRUE.
!          bmin = 1.E4
!          bold = 1.E4
! !******************CORRECTOR (FIELD LINE TRACING)
!          n = 0
!          SPAG_Loop_1_1: DO
!             p(1,3) = p(1,2) + step12*(5.*p(4,3)+8.*p(4,2)-p(4,1))
!             n = n + 1
!             p(2,3) = p(2,2) + step12*(5.*p(5,3)+8.*p(5,2)-p(5,1))
! !******************PREDICTOR (FIELD LINE TRACING)
!             p(1,4) = p(1,3) + step12*(23.*p(4,3)-16.*p(4,2)+5.*p(4,1))
!             p(2,4) = p(2,3) + step12*(23.*p(5,3)-16.*p(5,2)+5.*p(5,1))
!             p(3,4) = p(3,3) + step
!             CALL stoer(p(1,4),bq3,r3)
!             DO j = 1 , 3
!                DO i = 1 , 8
!                   p(i,j) = p(i,j+1)
!                ENDDO
!             ENDDO
!             b = sqrt(bq3)
!             IF ( b<bmin ) bmin = b
!             IF ( b<=bold ) THEN
!                bold = b
!                rold = 1./r3
!                Sp(1) = p(1,4)
!                Sp(2) = p(2,4)
!                Sp(3) = p(3,4)
!                CYCLE
!             ENDIF
!             IF ( bold/=bmin ) Value = .FALSE.
!             bdelta = (b-bold)/bold
!             IF ( bdelta>Bdel ) THEN
!                step = step/10.
!                spag_nextblock_1 = 2
!                CYCLE SPAG_DispatchLoop_1
!             ENDIF
!             EXIT SPAG_Loop_1_1
!          ENDDO SPAG_Loop_1_1
!          spag_nextblock_1 = 3
!       CASE (3)
!          Rr0 = rold
!          Bequ = bold
!          Bdel = bdelta
!          EXIT SPAG_DispatchLoop_1
!       END SELECT
!    ENDDO SPAG_DispatchLoop_1
! END SUBROUTINE findb0

    subroutine findb0(stps,bdel,value,bequ,rr0)

      REAL b , Bdel , bdelta , Bequ , bmin , bold , bq1 , bq2 , bq3 , p , r1 , r2 , r3 , &
           rold , Rr0 , Sp , step , step12 , Stps , zz
      INTEGER i , irun , j , n
      dimension p(8,4),sp(3)
      logical   value

      common/fidb0/ sp

      step=stps
      irun=0

      main : do
        irun=irun+1
        if (irun>5) then
            value=.false.
            exit main
        endif
        !*********************first three points 
        p(1,2)=sp(1)
        p(2,2)=sp(2)
        p(3,2)=sp(3)
        step=-sign(step,p(3,2))
        call stoer(p(1,2),bq2,r2)
        p(1,3)=p(1,2)+0.5*step*p(4,2)
        p(2,3)=p(2,2)+0.5*step*p(5,2)
        p(3,3)=p(3,2)+0.5*step
        call stoer(p(1,3),bq3,r3)
        p(1,1)=p(1,2)-step*(2.*p(4,2)-p(4,3))
        p(2,1)=p(2,2)-step*(2.*p(5,2)-p(5,3))
        p(3,1)=p(3,2)-step
        call stoer(p(1,1),bq1,r1)
        p(1,3)=p(1,2)+step*(20.*p(4,3)-3.*p(4,2)+p(4,1))/18.
        p(2,3)=p(2,2)+step*(20.*p(5,3)-3.*p(5,2)+p(5,1))/18.
        p(3,3)=p(3,2)+step
        call stoer(p(1,3),bq3,r3)
        !******************invert sense if required
        if (bq3>bq1) then
            step=-step
            r3=r1
            bq3=bq1
            do i=1,5
                zz=p(i,1)
                p(i,1)=p(i,3)
                p(i,3)=zz
            end do
        end if
        !******************initialization 
        step12=step/12.
        value=.true.
        bmin=1.e4
        bold=1.e4
        !******************corrector (field line tracing)
        n=0
        corrector : do
            p(1,3)=p(1,2)+step12*(5.*p(4,3)+8.*p(4,2)-p(4,1))
            n=n+1
            p(2,3)=p(2,2)+step12*(5.*p(5,3)+8.*p(5,2)-p(5,1))
            !******************predictor (field line tracing)
            p(1,4)=p(1,3)+step12*(23.*p(4,3)-16.*p(4,2)+5.*p(4,1))
            p(2,4)=p(2,3)+step12*(23.*p(5,3)-16.*p(5,2)+5.*p(5,1))
            p(3,4)=p(3,3)+step
            call stoer(p(1,4),bq3,r3)
            do j=1,3
                do i=1,8
                    p(i,j)=p(i,j+1)
                end do 
            end do
            b=sqrt(bq3)
            if (b<bmin) bmin=b
            if (b>bold) exit corrector
            bold=b
            rold=1./r3
            sp(1)=p(1,4)
            sp(2)=p(2,4)
            sp(3)=p(3,4)
        end do corrector
        if (bold/=bmin) then
            value=.false.
        endif
        bdelta=(b-bold)/bold
        if (bdelta<=bdel) exit main
        step=step/10.
    end do main

    rr0=rold
    bequ=bold
    bdel=bdelta
     
  end subroutine findb0

!--------------------------------------------------------------------
! CALCULATES L-VALUE FOR SPECIFIED GEODAETIC COORDINATES, ALTITUDE
! AND GEMAGNETIC FIELD MODEL.
! REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE
!      NO. 67, 1970.
!      G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972
!--------------------------------------------------------------------
! CHANGES (D. BILITZA, NOV 87):
!   - USING CORRECT DIPOL MOMENT I.E.,DIFFERENT COMMON/MODEL/
!   - USING IGRF EARTH MAGNETIC FIELD MODELS FROM 1945 TO 1990
!--------------------------------------------------------------------
!  INPUT:  ENTRY POINT SHELLG
!	 	GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
!         	GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
!         	ALT   ALTITUDE IN KM ABOVE SEA LEVEL
!
!	   ENTRY POINT SHELLC
!		V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
!			X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
!			Y-AXIS POINTING TO EQUATOR AT 90 LONG.
!			Z-AXIS POINTING TO NORTH POLE
!
!	   DIMO	      DIPOL MOMENT IN GAUSS (NORMALIZED TO EARTH RADIUS)
!
!	   COMMON
!		X(3)	NOT USED
!		H(144)	FIELD MODEL COEFFICIENTS ADJUSTED FOR SHELLG
!-----------------------------------------------------------------------
!  OUTPUT: FL   	L-VALUE
!	   ICODE  	=1 NORMAL COMPLETION
!			=2 UNPHYSICAL CONJUGATE POINT (FL MEANINGLESS)
!			=3 SHELL PARAMETER GREATER THAN LIMIT UP TO
!			   WHICH ACCURATE CALCULATION IS REQUIRED;
!			   APPROXIMATION IS USED.
! 	   B0   	MAGNETIC FIELD STRENGTH IN GAUSS
!-----------------------------------------------------------------------
  
  SUBROUTINE shellg(Glat,Glon,Alt,Dimo,Fl,Icode,B0)
   IMPLICIT NONE
   REAL Alt , Aquad , arg1 , arg2 , B0 , bequ , bq1 , bq2 , bq3 , Bquad , c0 , c1 , c2 , c3 , ct , d , d0 , d1 , d2 , Dimo
   REAL dimob0 , e0 , e1 , e2 , Era , ff , fi , Fl , gg , Glat , Glon , H , hli , oradik , oterm , p , r , r1 , r2 , r3
   REAL r3h , radik , rlat , rlon , rmax , rmin , rq , Sp , st , step , step12 , step2 , steq , stp , t , term , Umr , V , X
   REAL xx , z , zq , zz
   INTEGER i , Icode , iequ , n

   DIMENSION V(3) , p(8,100) , Sp(3)

   real,dimension(3,3),parameter ::  u = reshape([ +0.3511737 , -0.9148385 , -0.1993679 , &
                                                   +0.9335804 , +0.3583680 , +0.0000000 , &
                                                   +0.0714471 , -0.1861260 , +0.9799247], [3,3])

   COMMON X(3) , H(144)
   COMMON /fidb0 / Sp
   COMMON /gener / Umr , Era , Aquad , Bquad

   !-- RMIN, RMAX ARE BOUNDARIES FOR IDENTIFICATION OF ICODE=2 AND 3
   !-- STEP IS STEP SIZE FOR FIELD LINE TRACING
   !-- STEQ IS STEP SIZE FOR INTEGRATION

   DATA rmin , rmax/0.05 , 1.01/
   DATA step , steq/0.20 , 0.03/
   bequ = 1.E10

   !*****ENTRY POINT  SHELLG  TO BE USED WITH GEODETIC CO-ORDINATES
   rlat = Glat*Umr
   ct = sin(rlat)
   st = cos(rlat)
   d = sqrt(Aquad-(Aquad-Bquad)*ct*ct)
   X(1) = (Alt+Aquad/d)*st/Era
   X(3) = (Alt+Bquad/d)*ct/Era
   rlon = Glon*Umr
   X(2) = X(1)*sin(rlon)
   X(1) = X(1)*cos(rlon)
   CALL spag_block_1()
   RETURN

!*****ENTRY POINT  SHELLC  TO BE USED WITH CARTESIAN CO-ORDINATES
   ENTRY shellc(V,Fl,B0)
   X(1) = V(1)
   X(2) = V(2)
   X(3) = V(3)
   CALL spag_block_1()

CONTAINS

   SUBROUTINE spag_block_1

      integer,parameter :: max_loop_index = 100  ! 3333   <--- original code had 3333 ... was this a bug ????

      !*****CONVERT TO DIPOL-ORIENTED CO-ORDINATES
      rq = 1./(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      r3h = sqrt(rq*sqrt(rq))
      p(1,2) = (X(1)*u(1,1)+X(2)*u(2,1)+X(3)*u(3,1))*r3h
      p(2,2) = (X(1)*u(1,2)+X(2)*u(2,2))*r3h
      p(3,2) = (X(1)*u(1,3)+X(2)*u(2,3)+X(3)*u(3,3))*rq
      !     *****FIRST THREE POINTS OF FIELD LINE
      step = -sign(step,p(3,2))
      CALL stoer(p(1,2),bq2,r2)
      B0 = sqrt(bq2)
      p(1,3) = p(1,2) + 0.5*step*p(4,2)
      p(2,3) = p(2,2) + 0.5*step*p(5,2)
      p(3,3) = p(3,2) + 0.5*step
      CALL stoer(p(1,3),bq3,r3)
      p(1,1) = p(1,2) - step*(2.*p(4,2)-p(4,3))
      p(2,1) = p(2,2) - step*(2.*p(5,2)-p(5,3))
      p(3,1) = p(3,2) - step
      CALL stoer(p(1,1),bq1,r1)
      p(1,3) = p(1,2) + step*(20.*p(4,3)-3.*p(4,2)+p(4,1))/18.
      p(2,3) = p(2,2) + step*(20.*p(5,3)-3.*p(5,2)+p(5,1))/18.
      p(3,3) = p(3,2) + step
      CALL stoer(p(1,3),bq3,r3)
      !*****INVERT SENSE IF REQUIRED
      IF ( bq3>bq1 ) THEN
         step = -step
         r3 = r1
         bq3 = bq1
         DO i = 1 , 7
            zz = p(i,1)
            p(i,1) = p(i,3)
            p(i,3) = zz
         ENDDO
      ENDIF
      !*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
      IF ( bq1<bequ ) THEN
         bequ = bq1
         iequ = 1
      ENDIF
      IF ( bq2<bequ ) THEN
         bequ = bq2
         iequ = 2
      ENDIF
      IF ( bq3<bequ ) THEN
         bequ = bq3
         iequ = 3
      ENDIF
      !*****INITIALIZATION OF INTEGRATION LOOPS
      step12 = step/12.
      step2 = step + step
      steq = sign(steq,step)
      fi = 0.
      Icode = 1
      oradik = 0.
      oterm = 0.
      stp = r2*steq
      z = p(3,2) + stp
      stp = stp/0.75
      p(8,1) = step2*(p(1,1)*p(4,1)+p(2,1)*p(5,1))
      p(8,2) = step2*(p(1,2)*p(4,2)+p(2,2)*p(5,2))
      !*****MAIN LOOP (FIELD LINE TRACING)
      main: DO n = 3 , max_loop_index
         !*****CORRECTOR (FIELD LINE TRACING)
         p(1,n) = p(1,n-1) + step12*(5.*p(4,n)+8.*p(4,n-1)-p(4,n-2))
         p(2,n) = p(2,n-1) + step12*(5.*p(5,n)+8.*p(5,n-1)-p(5,n-2))
         !*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION
         !*****OF SLOWLY VARYING QUANTITIES
         p(8,n) = step2*(p(1,n)*p(4,n)+p(2,n)*p(5,n))
         c0 = p(1,n-1)**2 + p(2,n-1)**2
         c1 = p(8,n-1)
         c2 = (p(8,n)-p(8,n-2))*0.25
         c3 = (p(8,n)+p(8,n-2)-c1-c1)/6.0
         d0 = p(6,n-1)
         d1 = (p(6,n)-p(6,n-2))*0.5
         d2 = (p(6,n)+p(6,n-2)-d0-d0)*0.5
         e0 = p(7,n-1)
         e1 = (p(7,n)-p(7,n-2))*0.5
         e2 = (p(7,n)+p(7,n-2)-e0-e0)*0.5
         inner: DO
            !*****INNER LOOP (FOR QUADRATURE)
            t = (z-p(3,n-1))/step
            IF ( t>1. ) THEN
               !*****PREDICTOR (FIELD LINE TRACING)
               p(1,n+1) = p(1,n) + step12*(23.*p(4,n)-16.*p(4,n-1)+5.*p(4,n-2))
               p(2,n+1) = p(2,n) + step12*(23.*p(5,n)-16.*p(5,n-1)+5.*p(5,n-2))
               p(3,n+1) = p(3,n) + step
               CALL stoer(p(1,n+1),bq3,r3)
               !*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
               IF ( bq3<bequ ) THEN
                  iequ = n + 1
                  bequ = bq3
               ENDIF
               EXIT inner
            ELSE
               hli = 0.5*(((c3*t+c2)*t+c1)*t+c0)
               zq = z*z
               r = hli + sqrt(hli*hli+zq)
               IF ( r<=rmin ) THEN
                  !*****APPROXIMATION FOR HIGH VALUES OF L.
                  Icode = 3
                  t = -p(3,n-1)/step
                  Fl = 1./(abs(((c3*t+c2)*t+c1)*t+c0)+1E-15)
                  RETURN
               ENDIF
               rq = r*r
               ff = sqrt(1.+3.*zq/rq)
               radik = B0 - ((d2*t+d1)*t+d0)*r*rq*ff
               IF ( r>rmax ) THEN
                  Icode = 2
                  radik = radik - 12.*(r-rmax)**2
               ENDIF
               IF ( radik+radik<=oradik ) EXIT main
               term = sqrt(radik)*ff*((e2*t+e1)*t+e0)/(rq+zq)
               fi = fi + stp*(oterm+term)
               oradik = radik
               oterm = term
               stp = r*steq
               z = z + stp
            ENDIF
         ENDDO inner
      ENDDO main
      IF ( iequ<2 ) iequ = 2
      Sp(1) = p(1,iequ-1)
      Sp(2) = p(2,iequ-1)
      Sp(3) = p(3,iequ-1)
      IF ( oradik>=1E-15 ) fi = fi + stp/0.75*oterm*oradik/(oradik-radik)
      !
      !-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
      !-- because 1E-38 is the minimal allowable arg. for ALOG in our envir.
      !-- D. Bilitza, Nov 87.
      !
      fi = 0.5*abs(fi)/sqrt(B0) + 1E-12
      !*****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR.
      !
      !-- Correct dipole moment is used here. D. Bilitza, Nov 87.
      !
      dimob0 = Dimo/B0
      arg1 = alog(fi)
      arg2 = alog(dimob0)
!       arg = FI*FI*FI/DIMOB0
!       if(abs(arg).gt.88.0) arg=88.0
      xx = 3*arg1 - arg2
      IF ( xx>23.0 ) THEN
         gg = xx - 3.0460681E0
      ELSEIF ( xx>11.7 ) THEN
         gg = (((((2.8212095E-8*xx-3.8049276E-6)*xx+2.170224E-4)*xx-6.7310339E-3)*xx+1.2038224E-1)*xx-1.8461796E-1)                &
            & *xx + 2.0007187E0
      ELSEIF ( xx>+3.0 ) THEN
         gg = ((((((((6.3271665E-10*xx-3.958306E-8)*xx+9.9766148E-07)*xx-1.2531932E-5)*xx+7.9451313E-5)*xx-3.2077032E-4)           &
            & *xx+2.1680398E-3)*xx+1.2817956E-2)*xx+4.3510529E-1)*xx + 6.222355E-1
      ELSEIF ( xx>-3.0 ) THEN
         gg = ((((((((2.6047023E-10*xx+2.3028767E-9)*xx-2.1997983E-8)*xx-5.3977642E-7)*xx-3.3408822E-6)*xx+3.8379917E-5)           &
            & *xx+1.1784234E-3)*xx+1.4492441E-2)*xx+4.3352788E-1)*xx + 6.228644E-1
      ELSEIF ( xx>-22. ) THEN
         gg = ((((((((-8.1537735E-14*xx+8.3232531E-13)*xx+1.0066362E-9)*xx+8.1048663E-8)*xx+3.2916354E-6)*xx+8.2711096E-5)         &
            & *xx+1.3714667E-3)*xx+1.5017245E-2)*xx+4.3432642E-1)*xx + 6.2337691E-1
      ELSE
         gg = 3.33338E-1*xx + 3.0062102E-1
      ENDIF
      Fl = exp(alog((1.+exp(gg))*dimob0)/3.0)
      RETURN
   END SUBROUTINE spag_block_1

END SUBROUTINE shellg

!*==stoer.f90 processed by SPAG 8.01MH 09:18  3 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
SUBROUTINE stoer(P,Bq,R)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL Bq , dr , dsq , dx , dxm , dy , dym , dz , dzm , fli , H , P , q , R , rq , u , wr , Xi , xm , ym
   REAL zm
!*** End of declarations inserted by SPAG
!*******************************************************************
!* SUBROUTINE USED FOR FIELD LINE TRACING IN SHELLG                *
!* CALLS ENTRY POINT FELDI IN GEOMAGNETIC FIELD SUBROUTINE FELDG   *
!*******************************************************************
   DIMENSION P(7) , u(3,3)
   COMMON Xi(3) , H(144)
!*****XM,YM,ZM ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES
   zm = P(3)
   fli = P(1)*P(1) + P(2)*P(2) + 1E-15
   R = 0.5*(fli+sqrt(fli*fli+(zm+zm)**2))
   rq = R*R
   wr = sqrt(R)
   xm = P(1)*wr
   ym = P(2)*wr
!*****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM
   DATA u/ + 0.3511737 , -0.9148385 , -0.1993679 , +0.9335804 , +0.3583680 , +0.0000000 , +0.0714471 , -0.1861260 , +0.9799247/
   Xi(1) = xm*u(1,1) + ym*u(1,2) + zm*u(1,3)
   Xi(2) = xm*u(2,1) + ym*u(2,2) + zm*u(2,3)
   Xi(3) = xm*u(3,1) + zm*u(3,3)
!*****COMPUTE DERIVATIVES
! Changed from CALL FELDI(XI,H); XI, H are in COMMON block; results
! are the same; dkb Feb 1998
   CALL feldi
   q = H(1)/rq
   dx = H(3) + H(3) + q*Xi(1)
   dy = H(4) + H(4) + q*Xi(2)
   dz = H(2) + H(2) + q*Xi(3)
!*****TRANSFORM BACK TO GEOMAGNETIC CO-ORDINATE SYSTEM
   dxm = u(1,1)*dx + u(2,1)*dy + u(3,1)*dz
   dym = u(1,2)*dx + u(2,2)*dy
   dzm = u(1,3)*dx + u(2,3)*dy + u(3,3)*dz
   dr = (xm*dxm+ym*dym+zm*dzm)/R
!*****FORM SLOWLY VARYING EXPRESSIONS
   P(4) = (wr*dxm-0.5*P(1)*dr)/(R*dzm)
   P(5) = (wr*dym-0.5*P(2)*dr)/(R*dzm)
   dsq = rq*(dxm*dxm+dym*dym+dzm*dzm)
   Bq = dsq*rq*rq
   P(6) = sqrt(dsq/(rq+3.*zm*zm))
   P(7) = P(6)*(rq+zm*zm)/(rq*dzm)
END SUBROUTINE stoer

!*==feldg.f90 processed by SPAG 8.01MH 09:18  3 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
! SUBROUTINE feldg(Glat,Glon,Alt,Bnorth,Beast,Bdown,Babs)
!    IMPLICIT NONE
! !*** Start of declarations inserted by SPAG
!    REAL Alt , Aquad , B , Babs , Bdown , Beast , Bnorth , Bquad , brho , bxxx , byyy , bzzz , cp , ct , d , Era , f , G , Glat ,   &
!       & Glon
!    REAL H , rho , rlat , rlon , rq , s , sp , st , t , Time , Umr , V , x , Xi , xxx , y , yyy , z , zzz
!    INTEGER i , ih , ihmax , il , imax , is , k , last , m , Nmax
! !*** End of declarations inserted by SPAG
! !-------------------------------------------------------------------
! ! CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
! ! REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61,
! !      1970.
! !--------------------------------------------------------------------
! ! CHANGES (D. BILITZA, NOV 87):
! !   - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
! !   - CALCULATES DIPOL MOMENT
! !--------------------------------------------------------------------
! !  INPUT:  ENTRY POINT FELDG
! !	 	GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
! !         	GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
! !         	ALT   ALTITUDE IN KM ABOVE SEA LEVEL
! !
! !	   ENTRY POINT FELDC
! !		V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
! !			X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
! !			Y-AXIS POINTING TO EQUATOR AT 90 LONG.
! !			Z-AXIS POINTING TO NORTH POLE
! !
! !	   COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED
! !	     IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG.
! !
! !	   COMMON /MODEL/ AND /GENER/
! !		UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
! !		ERA	EARTH RADIUS FOR NORMALIZATION OF CARTESIAN
! !			COORDINATES (6371.2 KM)
! !		AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS FOR
! !			EARTH ELLIPSOID AS RECOMMENDED BY INTERNATIONAL
! !			ASTRONOMICAL UNION (6378.160, 6356.775 KM).
! !		NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
! !		TIME	YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC
! !			FIELD IS TO BE CALCULATED
! !		G(M)	NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
! !			M=NMAX*(NMAX+2)
! !------------------------------------------------------------------------
! !  OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS
! !	   BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
! !		  TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
! !		  POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
! !		  AND DOWNWARD.
! !-----------------------------------------------------------------------
!    DIMENSION V(3) , B(3)
!    CHARACTER*14 Name
!    COMMON Xi(3) , H(144)
!    COMMON /model / Name , Nmax , Time , G(144)
!    COMMON /gener / Umr , Era , Aquad , Bquad
!    INTEGER :: spag_nextblock_1
!    INTEGER :: spag_nextblock_2
!    spag_nextblock_1 = 1
!    SPAG_DispatchLoop_1: DO
!       SELECT CASE (spag_nextblock_1)
!       CASE (1)
! !
! !-- IS RECORDS ENTRY POINT
! !
! !*****ENTRY POINT  FELDG  TO BE USED WITH GEODETIC CO-ORDINATES
!          is = 1
!          rlat = Glat*Umr
!          ct = sin(rlat)
!          st = cos(rlat)
!          d = sqrt(Aquad-(Aquad-Bquad)*ct*ct)
!          rlon = Glon*Umr
!          cp = cos(rlon)
!          sp = sin(rlon)
!          zzz = (Alt+Bquad/d)*ct/Era
!          rho = (Alt+Aquad/d)*st/Era
!          xxx = rho*cp
!          yyy = rho*sp
!          spag_nextblock_1 = 2
!          CYCLE SPAG_DispatchLoop_1
! !
! !*****ENTRY POINT  FELDC  TO BE USED WITH CARTESIAN CO-ORDINATES
!          ENTRY feldc(V,B)
!          is = 2
!          xxx = V(1)
!          yyy = V(2)
!          zzz = V(3)
!          spag_nextblock_1 = 2
!       CASE (2)
!          rq = 1./(xxx*xxx+yyy*yyy+zzz*zzz)
!          Xi(1) = xxx*rq
!          Xi(2) = yyy*rq
!          Xi(3) = zzz*rq
!          spag_nextblock_1 = 3
!          CYCLE SPAG_DispatchLoop_1
! !
! !*****ENTRY POINT  FELDI  USED FOR L COMPUTATION
!          ENTRY feldi
!          is = 3
!          spag_nextblock_1 = 3
!       CASE (3)
!          ihmax = Nmax*Nmax + 1
!          last = ihmax + Nmax + Nmax
!          imax = Nmax + Nmax - 1
!          DO i = ihmax , last
!             H(i) = G(i)
!          ENDDO
!          DO k = 1 , 3 , 2
!             spag_nextblock_2 = 1
!             SPAG_DispatchLoop_2: DO
!                SELECT CASE (spag_nextblock_2)
!                CASE (1)
!                   i = imax
!                   ih = ihmax
!                   spag_nextblock_2 = 2
!                CASE (2)
!                   il = ih - i
!                   f = 2./float(i-k+2)
!                   x = Xi(1)*f
!                   y = Xi(2)*f
!                   z = Xi(3)*(f+f)
!                   i = i - 2
!                   IF ( i<1 ) THEN
!                      spag_nextblock_2 = 3
!                      CYCLE SPAG_DispatchLoop_2
!                   ENDIF
!                   IF ( i/=1 ) THEN
!                      DO m = 3 , i , 2
!                         H(il+m+1) = G(il+m+1) + z*H(ih+m+1) + x*(H(ih+m+3)-H(ih+m-1)) - y*(H(ih+m+2)+H(ih+m-2))
!                         H(il+m) = G(il+m) + z*H(ih+m) + x*(H(ih+m+2)-H(ih+m-2)) + y*(H(ih+m+3)+H(ih+m-1))
!                      ENDDO
!                   ENDIF
!                   H(il+2) = G(il+2) + z*H(ih+2) + x*H(ih+4) - y*(H(ih+3)+H(ih))
!                   H(il+1) = G(il+1) + z*H(ih+1) + y*H(ih+4) + x*(H(ih+3)-H(ih))
!                   spag_nextblock_2 = 3
!                CASE (3)
!                   H(il) = G(il) + z*H(ih) + 2.*(x*H(ih+1)+y*H(ih+2))
!                   ih = il
!                   IF ( i>=k ) THEN
!                      spag_nextblock_2 = 2
!                      CYCLE SPAG_DispatchLoop_2
!                   ENDIF
!                   EXIT SPAG_DispatchLoop_2
!                END SELECT
!             ENDDO SPAG_DispatchLoop_2
!          ENDDO
!          IF ( is==3 ) RETURN
!          s = .5*H(1) + 2.*(H(2)*Xi(3)+H(3)*Xi(1)+H(4)*Xi(2))
!          t = (rq+rq)*sqrt(rq)
!          bxxx = t*(H(3)-s*xxx)
!          byyy = t*(H(4)-s*yyy)
!          bzzz = t*(H(2)-s*zzz)
!          IF ( is==2 ) THEN
!             B(1) = bxxx
!             B(2) = byyy
!             B(3) = bzzz
!             RETURN
!          ENDIF
!          Babs = sqrt(bxxx*bxxx+byyy*byyy+bzzz*bzzz)
!          Beast = byyy*cp - bxxx*sp
!          brho = byyy*sp + bxxx*cp
!          Bnorth = bzzz*st - brho*ct
!          Bdown = -bzzz*ct - brho*st
!          RETURN
!       END SELECT
!    ENDDO SPAG_DispatchLoop_1
! END SUBROUTINE feldg

!.... NOTE: this one was broken by SPAG (above)... so revert to original and see
!           if we can manually refactor it.
subroutine feldg(glat,glon,alt,bnorth,beast,bdown,babs)           
   !-------------------------------------------------------------------
   ! calculates earth magnetic field from spherical harmonics model
   ! ref: g. kluge, european space operations centre, internal note 61, 
   !      1970.
   !--------------------------------------------------------------------
   ! changes (d. bilitza, nov 87):
   !   - field coefficients in binary data files instead of block data
   !   - calculates dipol moment
   !--------------------------------------------------------------------
   !  input:  entry point feldg
   !         glat  geodetic latitude in degrees (north)
   !             glon  geodetic longitude in degrees (east)
   !             alt   altitude in km above sea level
   !
   !       entry point feldc
   !        v(3)  cartesian coordinates in earth radii (6371.2 km)
   !            x-axis pointing to equator at 0 longitude
   !            y-axis pointing to equator at 90 long.
   !            z-axis pointing to north pole
   !
   !       common blank and entry point feldi are needed when used
   !         in connection with l-calculation program shellg.
   !    
   !       common /model/ and /gener/
   !        umr     = atan(1.0)*4./180.   <degree>*umr=<radiant>
   !        era    earth radius for normalization of cartesian 
   !            coordinates (6371.2 km)
   !        aquad, bquad   square of major and minor half axis for 
   !            earth ellipsoid as recommended by international 
   !            astronomical union (6378.160, 6356.775 km).
   !        nmax    maximum order of spherical harmonics
   !        time    year (decimal: 1973.5) for which magnetic 
   !            field is to be calculated
   !        g(m)    normalized field coefficients (see feldcof)
   !            m=nmax*(nmax+2)
   !------------------------------------------------------------------------
   !  output: babs   magnetic field strength in gauss
   !       bnorth, beast, bdown   components of the field with respect
   !          to the local geodetic coordinate system, with axis
   !          pointing in the tangential plane to the north, east
   !          and downward.   
   !-----------------------------------------------------------------------
   
   !*** Start of declarations inserted by SPAG
       REAL Alt , Aquad , B , Babs , Bdown , Beast , Bnorth , Bquad , brho , bxxx , &
            byyy , bzzz , cp , ct , d , Era , f , G , Glat , Glon
        REAL H , rho , rlat , rlon , rq , s , sp , st , t , Time , Umr , V , x , Xi , xxx , y , yyy , z , zzz
        INTEGER i , ih , ihmax , il , imax , is , k , last , m , Nmax
     !*** End of declarations inserted by SPAG
     
         dimension     v(3),b(3)   
         character*14     name

         common         xi(3),h(144)
         common/model/    name,nmax,time,g(144)  
         common/gener/    umr,era,aquad,bquad
   !
   !-- is records entry point
   !
   !*****entry point  feldg  to be used with geodetic co-ordinates         
         is=1                                                              
         rlat=glat*umr
         ct=sin(rlat)                                                      
         st=cos(rlat)                                                      
         d=sqrt(aquad-(aquad-bquad)*ct*ct)                                 
         rlon=glon*umr
         cp=cos(rlon)                                                      
         sp=sin(rlon)                                                      
          zzz=(alt+bquad/d)*ct/era
          rho=(alt+aquad/d)*st/era
          xxx=rho*cp                                                       
          yyy=rho*sp                                                       
          goto 10
   !
   !*****entry point  feldc  to be used with cartesian co-ordinates        
         entry feldc(v,b)                                                  
         is=2                                                              
         xxx=v(1)                                                          
         yyy=v(2)                                                          
         zzz=v(3)                                                          
   10    rq=1./(xxx*xxx+yyy*yyy+zzz*zzz) 
         xi(1)=xxx*rq                                                      
         xi(2)=yyy*rq                                                      
         xi(3)=zzz*rq                                                      
          goto 20                                                            
   !
   !*****entry point  feldi  used for l computation                        
         entry feldi                                                       
         is=3                                                              
   20    ihmax=nmax*nmax+1                                                 
         last=ihmax+nmax+nmax                                              
         imax=nmax+nmax-1                                                  
         do i=ihmax,last                                                 
             h(i)=g(i)    
         end do                                                     
         do k=1,3,2                                                      
            i=imax                                                            
            ih=ihmax                                                          
      1     il=ih-i                                                           
            f=2./float(i-k+2)                                                 
            x=xi(1)*f                                                         
            y=xi(2)*f                                                         
            z=xi(3)*(f+f)                                                     
            i=i-2         
            if ((i-1)>=0) then 
               if ((i-1)>0) then
                  do m=3,i,2                                                      
                     h(il+m+1)=g(il+m+1)+z*h(ih+m+1)+x*(h(ih+m+3)-h(ih+m-1))-y*(h(ih+m+2)+h(ih+m-2))           
                     h(il+m)=g(il+m)+z*h(ih+m)+x*(h(ih+m+2)-h(ih+m-2))+y*(h(ih+m+3)+h(ih+m-1))    
                  end do  
               end if 
               h(il+2)=g(il+2)+z*h(ih+2)+x*h(ih+4)-y*(h(ih+3)+h(ih))             
               h(il+1)=g(il+1)+z*h(ih+1)+y*h(ih+4)+x*(h(ih+3)-h(ih))   
            end if
            h(il)=g(il)+z*h(ih)+2.*(x*h(ih+1)+y*h(ih+2))
            ih=il                                                             
            if (i>=k) goto 1                                                   
         end do            

         if (is==3) return                                                 
         s=.5*h(1)+2.*(h(2)*xi(3)+h(3)*xi(1)+h(4)*xi(2))                   
         t=(rq+rq)*sqrt(rq)                                                
         bxxx=t*(h(3)-s*xxx)                                               
         byyy=t*(h(4)-s*yyy)                                               
         bzzz=t*(h(2)-s*zzz)                                               
         if (is==2) then
            b(1)=bxxx                                                         
            b(2)=byyy                                                         
            b(3)=bzzz                                                         
         else
            babs=sqrt(bxxx*bxxx+byyy*byyy+bzzz*bzzz)
            beast=byyy*cp-bxxx*sp                                             
            brho=byyy*sp+bxxx*cp                                              
            bnorth=bzzz*st-brho*ct                                            
            bdown=-bzzz*ct-brho*st                                            
         end if

   end subroutine feldg        

!*==feldcof.f90 processed by SPAG 8.01MH 09:18  3 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
SUBROUTINE feldcof(Year,Dimo)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL Aquad , Bquad , Dimo , dte1 , dte2 , dtemod , Erad , Gh1 , gh2 , gha , sqrt2 , Time , Umr , Year
   INTEGER i , ier , is , iu , iyea , j , l , m , n , Nmax , nmax1 , nmax2 , numye
!*** End of declarations inserted by SPAG
!------------------------------------------------------------------------
!  DETERMINES COEFFICIENTS AND DIPOL MOMENT FROM IGRF MODELS
!
!	INPUT:  YEAR	DECIMAL YEAR FOR WHICH GEOMAGNETIC FIELD IS TO
!			BE CALCULATED (e.g.:1995.5 for day 185 of 1995)
!	OUTPUT:	DIMO	GEOMAGNETIC DIPOL MOMENT IN GAUSS (NORMALIZED
!			TO EARTH'S RADIUS) AT THE TIME (YEAR)
!  D. BILITZA, NSSDC, GSFC, CODE 633, GREENBELT, MD 20771,
!	(301)286-9536   NOV 1987.
!  -corrected for 2000 update - dkb- 5/31/2000
!  ### updated to IGRF-2000 version -dkb- 5/31/2000
!  ### updated to IGRF-2005 version -dkb- 3/24/2000
!-----------------------------------------------------------------------
   CHARACTER*14 filmod , Fil1 , fil2
! ### FILMOD, DTEMOD arrays +1
   DIMENSION Gh1(144) , gh2(120) , gha(144) , filmod(17) , dtemod(17)
   DOUBLE PRECISION x , f0 , f
   COMMON /model/ Fil1 , Nmax , Time , Gh1
   COMMON /gener/ Umr , Erad , Aquad , Bquad
! ### changed to conform with IGRF 45-95, also FILMOD, DTEMOD arrays +1
   DATA filmod/'dgrf1945.dat' , 'dgrf1950.dat' , 'dgrf1955.dat' , 'dgrf1960.dat' , 'dgrf1965.dat' , 'dgrf1970.dat' ,               &
      & 'dgrf1975.dat' , 'dgrf1980.dat' , 'dgrf1985.dat' , 'dgrf1990.dat' , 'dgrf1995.dat' , 'dgrf2000.dat' , 'dgrf2005.dat' ,     &
       &'dgrf2010.dat' , 'dgrf2015.dat' , 'igrf2020.dat' , 'igrf2020s.dat'/
   DATA dtemod/1945. , 1950. , 1955. , 1960. , 1965. , 1970. , 1975. , 1980. , 1985. , 1990. , 1995. , 2000. , 2005. , 2010. ,     &
      & 2015. , 2020. , 2025./
!
! ### numye is number of 5-year priods represented by IGRF
!
   numye = 16
!
!  IS=0 FOR SCHMIDT NORMALIZATION   IS=1 GAUSS NORMALIZATION
!  IU  IS INPUT UNIT NUMBER FOR IGRF COEFFICIENT SETS
!
   iu = 10
   is = 0
!-- DETERMINE IGRF-YEARS FOR INPUT-YEAR
   Time = Year
   iyea = int(Year/5.)*5
   l = (iyea-1945)/5 + 1
   IF ( l<1 ) l = 1
   IF ( l>numye ) l = numye
   dte1 = dtemod(l)
   Fil1 = filmod(l)
   dte2 = dtemod(l+1)
   fil2 = filmod(l+1)
!-- GET IGRF COEFFICIENTS FOR THE BOUNDARY YEARS
   CALL getshc(iu,Fil1,nmax1,Erad,Gh1,ier)
   IF ( ier/=0 ) STOP
   CALL getshc(iu,fil2,nmax2,Erad,gh2,ier)
   IF ( ier/=0 ) STOP
!-- DETERMINE IGRF COEFFICIENTS FOR YEAR
   IF ( l<=numye-1 ) THEN
      CALL intershc(Year,dte1,nmax1,Gh1,dte2,nmax2,gh2,Nmax,gha)
   ELSE
      CALL extrashc(Year,dte1,nmax1,Gh1,nmax2,gh2,Nmax,gha)
   ENDIF
!-- DETERMINE MAGNETIC DIPOL MOMENT AND COEFFIECIENTS G
   f0 = 0.D0
   DO j = 1 , 3
      f = gha(j)*1.D-5
      f0 = f0 + f*f
   ENDDO
   Dimo = sqrt(f0)
 
   Gh1(1) = 0.0
   i = 2
   f0 = 1.D-5
   IF ( is==0 ) f0 = -f0
   sqrt2 = sqrt(2.)
 
   DO n = 1 , Nmax
      x = n
      f0 = f0*x*x/(4.D0*x-2.D0)
      IF ( is==0 ) f0 = f0*(2.D0*x-1.D0)/x
      f = f0*0.5D0
      IF ( is==0 ) f = f*sqrt2
      Gh1(i) = gha(i-1)*f0
      i = i + 1
      DO m = 1 , n
         f = f*(x+m)/(x-m+1.D0)
         IF ( is==0 ) f = f*sqrt((x-m+1.D0)/(x+m))
         Gh1(i) = gha(i-1)*f
         Gh1(i+1) = gha(i)*f
         i = i + 2
      ENDDO
   ENDDO
END SUBROUTINE feldcof

!*==getshc.f90 processed by SPAG 8.01MH 09:18  3 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
SUBROUTINE getshc(Iu,Fspec,Nmax,Erad,Gh,Ier)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL Erad , g , Gh , h
   INTEGER i , Ier , Iu , m , mm , n , Nmax , nn
!*** End of declarations inserted by SPAG
 
! ===============================================================
!
!	Version 1.01
!
!	Reads spherical harmonic coefficients from the specified
!	file into an array.
!
!	Input:
!	    IU    - Logical unit number
!	    FSPEC - File specification
!
!	Output:
!	    NMAX  - Maximum degree and order of model
!	    ERAD  - Earth's radius associated with the spherical
!		    harmonic coefficients, in the same units as
!		    elevation
!	    GH    - Schmidt quasi-normal internal spherical
!		    harmonic coefficients
!	    IER   - Error number: =  0, no error
!				  = -2, records out of order
!			     	  = FORTRAN run-time error number
!
!	A. Zunde
!	USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
!
! ===============================================================
 
   CHARACTER Fspec*(*)
   DIMENSION Gh(*)
 
! ---------------------------------------------------------------
!	Open coefficient file. Read past first header record.
!	Read degree and order of model and Earth's radius.
! ---------------------------------------------------------------
   OPEN (Iu,FILE=Fspec,STATUS='OLD',IOSTAT=Ier,ERR=100)
   READ (Iu,*,IOSTAT=Ier,ERR=100)
   READ (Iu,*,IOSTAT=Ier,ERR=100) Nmax , Erad
! ---------------------------------------------------------------
!	Read the coefficient file, arranged as follows:
!
!					N     M     G     H
!					----------------------
!				    /   1     0    GH(1)  -
!				   /	1     1    GH(2) GH(3)
!				  /	2     0    GH(4)  -
!				 /	2     1    GH(5) GH(6)
!	    NMAX*(NMAX+3)/2 	/	2     2    GH(7) GH(8)
!	       records		\	3     0    GH(9)  -
!				 \      .     .     .     .
!				  \	.     .     .     .
!	    NMAX*(NMAX+2)	   \	.     .     .     .
!	    elements in GH	    \  NMAX  NMAX   .     .
!
!	N and M are, respectively, the degree and order of the
!	coefficient.
! ---------------------------------------------------------------
 
   i = 0
   main: DO nn = 1 , Nmax
      DO mm = 0 , nn
         READ (Iu,*,IOSTAT=Ier,ERR=100) n , m , g , h
         IF ( nn/=n .OR. mm/=m ) THEN
            Ier = -2
            EXIT main
         ENDIF
         i = i + 1
         Gh(i) = g
         IF ( m/=0 ) THEN
            i = i + 1
            Gh(i) = h
         ENDIF
      ENDDO
   ENDDO main
 
 100  CLOSE (Iu)
 
END SUBROUTINE getshc

!*==intershc.f90 processed by SPAG 8.01MH 09:18  3 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
SUBROUTINE intershc(Date,Dte1,Nmax1,Gh1,Dte2,Nmax2,Gh2,Nmax,Gh)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL Date , Dte1 , Dte2 , factor , Gh , Gh1 , Gh2
   INTEGER i , k , l , Nmax , Nmax1 , Nmax2
!*** End of declarations inserted by SPAG
 
! ===============================================================
!
!	Version 1.01
!
!	Interpolates linearly, in time, between two spherical
!	harmonic models.
!
!	Input:
!	    DATE  - Date of resulting model (in decimal year)
!	    DTE1  - Date of earlier model
!	    NMAX1 - Maximum degree and order of earlier model
!	    GH1   - Schmidt quasi-normal internal spherical
!		    harmonic coefficients of earlier model
!	    DTE2  - Date of later model
!	    NMAX2 - Maximum degree and order of later model
!	    GH2   - Schmidt quasi-normal internal spherical
!		    harmonic coefficients of later model
!
!	Output:
!	    GH    - Coefficients of resulting model
!	    NMAX  - Maximum degree and order of resulting model
!
!	A. Zunde
!	USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
!
! ===============================================================
 
   DIMENSION Gh1(*) , Gh2(*) , Gh(*)
 
! ---------------------------------------------------------------
!	The coefficients (GH) of the resulting model, at date
!	DATE, are computed by linearly interpolating between the
!	coefficients of the earlier model (GH1), at date DTE1,
!	and those of the later model (GH2), at date DTE2. If one
!	model is smaller than the other, the interpolation is
!	performed with the missing coefficients assumed to be 0.
! ---------------------------------------------------------------
 
   factor = (Date-Dte1)/(Dte2-Dte1)
 
   IF ( Nmax1==Nmax2 ) THEN
      k = Nmax1*(Nmax1+2)
      Nmax = Nmax1
   ELSEIF ( Nmax1>Nmax2 ) THEN
      k = Nmax2*(Nmax2+2)
      l = Nmax1*(Nmax1+2)
      DO i = k + 1 , l
         Gh(i) = Gh1(i) + factor*(-Gh1(i))
      ENDDO
      Nmax = Nmax1
   ELSE
      k = Nmax1*(Nmax1+2)
      l = Nmax2*(Nmax2+2)
      DO i = k + 1 , l
         Gh(i) = factor*Gh2(i)
      ENDDO
      Nmax = Nmax2
   ENDIF
 
   DO i = 1 , k
      Gh(i) = Gh1(i) + factor*(Gh2(i)-Gh1(i))
   ENDDO
 
END SUBROUTINE intershc

!*==extrashc.f90 processed by SPAG 8.01MH 09:18  3 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
SUBROUTINE extrashc(Date,Dte1,Nmax1,Gh1,Nmax2,Gh2,Nmax,Gh)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL Date , Dte1 , factor , Gh , Gh1 , Gh2
   INTEGER i , k , l , Nmax , Nmax1 , Nmax2
!*** End of declarations inserted by SPAG
 
! ===============================================================
!
!	Version 1.01
!
!	Extrapolates linearly a spherical harmonic model with a
!	rate-of-change model.
!
!	Input:
!	    DATE  - Date of resulting model (in decimal year)
!	    DTE1  - Date of base model
!	    NMAX1 - Maximum degree and order of base model
!	    GH1   - Schmidt quasi-normal internal spherical
!		    harmonic coefficients of base model
!	    NMAX2 - Maximum degree and order of rate-of-change
!		    model
!	    GH2   - Schmidt quasi-normal internal spherical
!		    harmonic coefficients of rate-of-change model
!
!	Output:
!	    GH    - Coefficients of resulting model
!	    NMAX  - Maximum degree and order of resulting model
!
!	A. Zunde
!	USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
!
! ===============================================================
 
   DIMENSION Gh1(*) , Gh2(*) , Gh(*)
 
! ---------------------------------------------------------------
!	The coefficients (GH) of the resulting model, at date
!	DATE, are computed by linearly extrapolating the coef-
!	ficients of the base model (GH1), at date DTE1, using
!	those of the rate-of-change model (GH2), at date DTE2. If
!	one model is smaller than the other, the extrapolation is
!	performed with the missing coefficients assumed to be 0.
! ---------------------------------------------------------------
 
   factor = (Date-Dte1)
 
   IF ( Nmax1==Nmax2 ) THEN
      k = Nmax1*(Nmax1+2)
      Nmax = Nmax1
   ELSEIF ( Nmax1>Nmax2 ) THEN
      k = Nmax2*(Nmax2+2)
      l = Nmax1*(Nmax1+2)
      DO i = k + 1 , l
         Gh(i) = Gh1(i)
      ENDDO
      Nmax = Nmax1
   ELSE
      k = Nmax1*(Nmax1+2)
      l = Nmax2*(Nmax2+2)
      DO i = k + 1 , l
         Gh(i) = factor*Gh2(i)
      ENDDO
      Nmax = Nmax2
   ENDIF
 
   DO i = 1 , k
      Gh(i) = Gh1(i) + factor*Gh2(i)
   ENDDO
 
END SUBROUTINE extrashc

!*==initize.f90 processed by SPAG 8.01MH 09:18  3 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!
SUBROUTINE initize()
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL Aquad , Bquad , Era , erequ , erpol , Umr
!*** End of declarations inserted by SPAG
!----------------------------------------------------------------
! Initializes the parameters in COMMON/GENER/
!
!	UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
!	ERA	EARTH RADIUS FOR NORMALIZATION OF CARTESIAN
!			COORDINATES (6371.2 KM)
!	EREQU	MAJOR HALF AXIS FOR EARTH ELLIPSOID (6378.160 KM)
!	ERPOL	MINOR HALF AXIS FOR EARTH ELLIPSOID (6356.775 KM)
!	AQUAD	SQUARE OF MAJOR HALF AXIS FOR EARTH ELLIPSOID
!	BQUAD   SQUARE OF MINOR HALF AXIS FOR EARTH ELLIPSOID
!
! ERA, EREQU and ERPOL as recommended by the INTERNATIONAL
! ASTRONOMICAL UNION .
!-----------------------------------------------------------------
   COMMON /gener / Umr , Era , Aquad , Bquad
   Era = 6371.2
   erequ = 6378.16
   erpol = 6356.775
   Aquad = erequ*erequ
   Bquad = erpol*erpol
   Umr = atan(1.0)*4./180.
END SUBROUTINE initize

end module SHELLIG_module