!*****************************************************************************************
!>
!  IGRF model
!
!### History
!  * SHELLIG.FOR, Version 2.0, January 1992
!  * 11/01/91-DKB- SHELLG: lowest starting point for B0 search is 2
!  * 1/27/92-DKB- Adopted to IGRF-91 coefficients model
!  * 2/05/92-DKB- Reduce variable-names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE
!  * 8/08/95-DKB- Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S
!  * 5/31/00-DKB- Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s
!  * 3/24/05-DKB- Updated to IGRF-45-10; new coeff.: IGRF05, IGRF05s

   module shellig_module

      use radbelt_kinds_module

      implicit none

      private

      integer,parameter :: filename_len = 14 !! length of the model data file names

      ! parameters formerly in `gener` common block
      real(wp),parameter :: Era = 6371.2_wp !! earth radius for normalization of cartesian coordinates (6371.2 km)
      real(wp),parameter :: erequ = 6378.16_wp
      real(wp),parameter :: erpol = 6356.775_wp
      real(wp),parameter :: Aquad = erequ*erequ !! square of major half axis for
                                                !! earth ellipsoid as recommended by international
                                                !! astronomical union
      real(wp),parameter :: Bquad = erpol*erpol !! square of minor half axis for
                                                !! earth ellipsoid as recommended by international
                                                !! astronomical union
      real(wp),parameter :: Umr = atan(1.0_wp)*4.0_wp/180.0_wp !! atan(1.0)*4./180.   <degree>*umr=<radiant>

      real(wp),dimension(3,3),parameter ::  u = reshape([ +0.3511737_wp , -0.9148385_wp , -0.1993679_wp , &
                                                          +0.9335804_wp , +0.3583680_wp , +0.0000000_wp , &
                                                          +0.0714471_wp , -0.1861260_wp , +0.9799247_wp], [3,3])

       type,public :: shellig_type
         private

         ! formerly in the `fidb0` common block
         real(wp),dimension(3) :: sp = 0.0_wp

         ! formerly in blank common
         real(wp),dimension(3) :: xi = 0.0_wp
         real(wp),dimension(144) :: h = 0.0_wp !! Field model coefficients adjusted for [[shellg]]

         ! formerly in `model` common block
         integer :: iyea = 0 !! the int year corresponding to the file `name` that has been read
         character(len=filename_len) :: name = '' !! file name
         integer :: nmax = 0 !! maximum order of spherical harmonics
         real(wp) :: Time = 0.0_wp !! year (decimal: 1973.5) for which magnetic field is to be calculated
         real(wp),dimension(144) :: g = 0.0_wp  !! `g(m)` -- normalized field coefficients (see [[feldcof]]) m=nmax*(nmax+2)

         ! formerly saved vars in shellg:
         real(wp) :: step = 0.20_wp !! step size for field line tracing
         real(wp) :: steq = 0.03_wp !! step size for integration

         ! from feldcof, so we can cache the coefficients
         real(wp),dimension(120) :: gh2 = 0.0_wp   ! JW : why is this 120 and g is 144 ???

         contains
         private

         procedure,public :: igrf

         procedure, public :: feldcof
         procedure, public :: feldg
         procedure, public :: shellg
         procedure, public :: findb0
         procedure :: stoer, getshc, intershc, extrashc

      end type shellig_type

   contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Wrapper for IGRF functions.

   subroutine igrf(me,lon,lat,height,year,xl,bbx)

      class(shellig_type),intent(inout) :: me
      real(wp),intent(in) :: lon !! geodetic longitude in degrees (east)
      real(wp),intent(in) :: lat !! geodetic latitude in degrees (north)
      real(wp),intent(in) :: height !! altitude in km above sea level
      real(wp),intent(in) :: year !! decimal year for which geomagnetic field is to
                                  !! be calculated (e.g.:1995.5 for day 185 of 1995)
      real(wp),intent(out) :: xl !! l-value
      real(wp),intent(out) :: bbx

      real(wp) :: bab1 , babs , bdel , bdown , beast , &
                  beq , bequ , bnorth , dimo , rr0
      integer :: icode
      logical :: val

      real(wp),parameter :: stps = 0.05_wp

      call me%feldcof(year,dimo)
      call me%feldg(lat,lon,height,bnorth,beast,bdown,babs)
      call me%shellg(lat,lon,height,dimo,xl,icode,bab1)

      bequ = dimo/(xl*xl*xl)
      if ( icode==1 ) then
         bdel = 1.0e-3_wp
         call me%findb0(stps,bdel,val,beq,rr0)
         if ( val ) bequ = beq
      endif
      bbx = babs/bequ

   end subroutine igrf
!*****************************************************************************************

!*****************************************************************************************
!>
    subroutine findb0(me,stps,bdel,value,bequ,rr0)

      class(shellig_type),intent(inout) :: me
      real(wp),intent(in) :: stps
      real(wp),intent(inout) :: bdel
      real(wp),intent(out) :: bequ
      logical,intent(out) :: value
      real(wp),intent(out) :: rr0

      real(wp) :: b , bdelta , bmin , bold , bq1 , &
                  bq2 , bq3 , p(8,4) , r1 , r2 , r3 , &
                  rold , step , step12 , zz
      integer :: i , irun , j , n

      step=stps
      irun=0

      main : do
        irun=irun+1
        if (irun>5) then
            value=.false.
            exit main
        endif
        !*********************first three points
        p(1,2)=me%sp(1)
        p(2,2)=me%sp(2)
        p(3,2)=me%sp(3)
        step=-sign(step,p(3,2))
        call me%stoer(p(1,2),bq2,r2)
        p(1,3)=p(1,2)+0.5_wp*step*p(4,2)
        p(2,3)=p(2,2)+0.5_wp*step*p(5,2)
        p(3,3)=p(3,2)+0.5_wp*step
        call me%stoer(p(1,3),bq3,r3)
        p(1,1)=p(1,2)-step*(2.0_wp*p(4,2)-p(4,3))
        p(2,1)=p(2,2)-step*(2.0_wp*p(5,2)-p(5,3))
        p(3,1)=p(3,2)-step
        call me%stoer(p(1,1),bq1,r1)
        p(1,3)=p(1,2)+step*(20.0_wp*p(4,3)-3.*p(4,2)+p(4,1))/18.0_wp
        p(2,3)=p(2,2)+step*(20.0_wp*p(5,3)-3.*p(5,2)+p(5,1))/18.0_wp
        p(3,3)=p(3,2)+step
        call me%stoer(p(1,3),bq3,r3)
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
        step12=step/12.0_wp
        value=.true.
        bmin=1.0e4_wp
        bold=1.0e4_wp
        !******************corrector (field line tracing)
        n=0
        corrector : do
            p(1,3)=p(1,2)+step12*(5.0_wp*p(4,3)+8.0_wp*p(4,2)-p(4,1))
            n=n+1
            p(2,3)=p(2,2)+step12*(5.0_wp*p(5,3)+8.0_wp*p(5,2)-p(5,1))
            !******************predictor (field line tracing)
            p(1,4)=p(1,3)+step12*(23.0_wp*p(4,3)-16.0_wp*p(4,2)+5.0_wp*p(4,1))
            p(2,4)=p(2,3)+step12*(23.0_wp*p(5,3)-16.0_wp*p(5,2)+5.0_wp*p(5,1))
            p(3,4)=p(3,3)+step
            call me%stoer(p(1,4),bq3,r3)
            do j=1,3
                do i=1,8
                    p(i,j)=p(i,j+1)
                end do
            end do
            b=sqrt(bq3)
            if (b<bmin) bmin=b
            if (b>bold) exit corrector
            bold=b
            rold=1.0_wp/r3
            me%sp(1)=p(1,4)
            me%sp(2)=p(2,4)
            me%sp(3)=p(3,4)
        end do corrector
        if (bold/=bmin) value=.false.
        bdelta=(b-bold)/bold
        if (bdelta<=bdel) exit main
        step=step/10.0_wp
    end do main

    rr0=rold
    bequ=bold
    bdel=bdelta

  end subroutine findb0

!*****************************************************************************************
!>
!  calculates l-value for specified geodaetic coordinates, altitude
!  and gemagnetic field model.
!
!### Reference
!  * G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE
!    NO. 67, 1970.
!  * G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972
!
!### History
! * CHANGES (D. BILITZA, NOV 87):
!   - USING CORRECT DIPOL MOMENT I.E.,DIFFERENT COMMON/MODEL/
!   - USING IGRF EARTH MAGNETIC FIELD MODELS FROM 1945 TO 1990

  subroutine shellg(me,glat,glon,alt,dimo,fl,icode,b0)

   class(shellig_type),intent(inout) :: me
   real(wp),intent(in) :: glat !! GEODETIC LATITUDE IN DEGREES (NORTH)
   real(wp),intent(in) :: glon !! GEODETIC LONGITUDE IN DEGREES (EAST)
   real(wp),intent(in) :: alt  !! ALTITUDE IN KM ABOVE SEA LEVEL
   real(wp),intent(in) :: dimo !! DIPOL MOMENT IN GAUSS (NORMALIZED TO EARTH RADIUS)
   real(wp),intent(out) :: fl  !! l-value
   integer,intent(out) :: icode  !! * =1 normal completion
                                 !! * =2 unphysical conjugate point (fl meaningless)
                                 !! * =3 shell parameter greater than limit up to
                                 !!   which accurate calculation is required;
                                 !!   approximation is used.
   real(wp),intent(out) :: b0 !! magnetic field strength in gauss

   real(wp) :: arg1 , arg2 , bequ , bq1 , bq2 , bq3 , c0 , c1 , c2 , c3 , &
               ct , d , d0 , d1 , d2, dimob0 , e0 , e1 , e2 , ff , fi , gg , &
               hli , oradik , oterm , p(8,100) , r , r1 , r2 , r3 , r3h , radik , &
               rlat , rlon , rq , st , step12 , step2 , &
               stp , t , term , v(3) , xx , z , zq , zz
   integer :: i , iequ , n

   real(wp),parameter :: rmin = 0.05_wp !! boundaries for identification of `icode=2 and 3`
   real(wp),parameter :: rmax = 1.01_wp !! boundaries for identification of `icode=2 and 3`

   bequ = 1.0e10_wp

   !*****ENTRY POINT SHELLG TO BE USED WITH GEODETIC CO-ORDINATES
   rlat = Glat*Umr
   ct = sin(rlat)
   st = cos(rlat)
   d = sqrt(Aquad-(Aquad-Bquad)*ct*ct)
   me%Xi(1) = (Alt+Aquad/d)*st/Era
   me%Xi(3) = (Alt+Bquad/d)*ct/Era
   rlon = Glon*Umr
   me%Xi(2) = me%Xi(1)*sin(rlon)
   me%Xi(1) = me%Xi(1)*cos(rlon)
   CALL spag_block_1()
   RETURN

!*****ENTRY POINT SHELLC TO BE USED WITH CARTESIAN CO-ORDINATES   ! JW : does not appear to be used anywhere.
!    V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
!      X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
!      Y-AXIS POINTING TO EQUATOR AT 90 LONG.
!      Z-AXIS POINTING TO NORTH POLE
   ENTRY shellc(me,V,Fl,B0)
   me%Xi(1) = V(1)
   me%Xi(2) = V(2)
   me%Xi(3) = V(3)
   CALL spag_block_1()

CONTAINS

   subroutine spag_block_1

      integer,parameter :: max_loop_index = 100  ! 3333   <--- JW : original code had 3333 ... was this a bug ????

      !*****CONVERT TO DIPOL-ORIENTED CO-ORDINATES
      rq = 1.0_wp/(me%Xi(1)*me%Xi(1)+me%Xi(2)*me%Xi(2)+me%Xi(3)*me%Xi(3))
      r3h = sqrt(rq*sqrt(rq))
      p(1,2) = (me%Xi(1)*u(1,1)+me%Xi(2)*u(2,1)+me%Xi(3)*u(3,1))*r3h
      p(2,2) = (me%Xi(1)*u(1,2)+me%Xi(2)*u(2,2))*r3h
      p(3,2) = (me%Xi(1)*u(1,3)+me%Xi(2)*u(2,3)+me%Xi(3)*u(3,3))*rq
      !     *****FIRST THREE POINTS OF FIELD LINE
      me%step = -sign(me%step,p(3,2))
      call me%stoer(p(1,2),bq2,r2)
      B0 = sqrt(bq2)
      p(1,3) = p(1,2) + 0.5_wp*me%step*p(4,2)
      p(2,3) = p(2,2) + 0.5_wp*me%step*p(5,2)
      p(3,3) = p(3,2) + 0.5_wp*me%step
      call me%stoer(p(1,3),bq3,r3)
      p(1,1) = p(1,2) - me%step*(2.0_wp*p(4,2)-p(4,3))
      p(2,1) = p(2,2) - me%step*(2.0_wp*p(5,2)-p(5,3))
      p(3,1) = p(3,2) - me%step
      call me%stoer(p(1,1),bq1,r1)
      p(1,3) = p(1,2) + me%step*(20.0_wp*p(4,3)-3.*p(4,2)+p(4,1))/18.0_wp
      p(2,3) = p(2,2) + me%step*(20.0_wp*p(5,3)-3.*p(5,2)+p(5,1))/18.0_wp
      p(3,3) = p(3,2) + me%step
      call me%stoer(p(1,3),bq3,r3)
      !*****INVERT SENSE IF REQUIRED
      IF ( bq3>bq1 ) THEN
         me%step = -me%step
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
      step12 = me%step/12.0_wp
      step2 = me%step + me%step
      me%steq = sign(me%steq,me%step)
      fi = 0.0_wp
      Icode = 1
      oradik = 0.0_wp
      oterm = 0.0_wp
      stp = r2*me%steq
      z = p(3,2) + stp
      stp = stp/0.75_wp
      p(8,1) = step2*(p(1,1)*p(4,1)+p(2,1)*p(5,1))
      p(8,2) = step2*(p(1,2)*p(4,2)+p(2,2)*p(5,2))
      !*****MAIN LOOP (FIELD LINE TRACING)
      main: DO n = 3 , max_loop_index
         !*****CORRECTOR (FIELD LINE TRACING)
         p(1,n) = p(1,n-1) + step12*(5.0_wp*p(4,n)+8.0_wp*p(4,n-1)-p(4,n-2))
         p(2,n) = p(2,n-1) + step12*(5.0_wp*p(5,n)+8.0_wp*p(5,n-1)-p(5,n-2))
         !*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION
         !*****OF SLOWLY VARYING QUANTITIES
         p(8,n) = step2*(p(1,n)*p(4,n)+p(2,n)*p(5,n))
         c0 = p(1,n-1)**2 + p(2,n-1)**2
         c1 = p(8,n-1)
         c2 = (p(8,n)-p(8,n-2))*0.25_wp
         c3 = (p(8,n)+p(8,n-2)-c1-c1)/6.0_wp
         d0 = p(6,n-1)
         d1 = (p(6,n)-p(6,n-2))*0.5_wp
         d2 = (p(6,n)+p(6,n-2)-d0-d0)*0.5_wp
         e0 = p(7,n-1)
         e1 = (p(7,n)-p(7,n-2))*0.5_wp
         e2 = (p(7,n)+p(7,n-2)-e0-e0)*0.5_wp
         inner: DO
            !*****INNER LOOP (FOR QUADRATURE)
            t = (z-p(3,n-1))/me%step
            IF ( t>1.0_wp ) THEN
               !*****PREDICTOR (FIELD LINE TRACING)
               p(1,n+1) = p(1,n) + step12*(23.0_wp*p(4,n)-16.0_wp*p(4,n-1)+5.0_wp*p(4,n-2))
               p(2,n+1) = p(2,n) + step12*(23.0_wp*p(5,n)-16.0_wp*p(5,n-1)+5.0_wp*p(5,n-2))
               p(3,n+1) = p(3,n) + me%step
               call me%stoer(p(1,n+1),bq3,r3)
               !*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
               IF ( bq3<bequ ) THEN
                  iequ = n + 1
                  bequ = bq3
               ENDIF
               EXIT inner
            ELSE
               hli = 0.5_wp*(((c3*t+c2)*t+c1)*t+c0)
               zq = z*z
               r = hli + sqrt(hli*hli+zq)
               IF ( r<=rmin ) THEN
                  !*****APPROXIMATION FOR HIGH VALUES OF L.
                  Icode = 3
                  t = -p(3,n-1)/me%step
                  Fl = 1.0_wp/(abs(((c3*t+c2)*t+c1)*t+c0)+1.0e-15_wp)
                  RETURN
               ENDIF
               rq = r*r
               ff = sqrt(1.0_wp+3.0_wp*zq/rq)
               radik = B0 - ((d2*t+d1)*t+d0)*r*rq*ff
               IF ( r>rmax ) THEN
                  Icode = 2
                  radik = radik - 12.0_wp*(r-rmax)**2
               ENDIF
               IF ( radik+radik<=oradik ) EXIT main
               term = sqrt(radik)*ff*((e2*t+e1)*t+e0)/(rq+zq)
               fi = fi + stp*(oterm+term)
               oradik = radik
               oterm = term
               stp = r*me%steq
               z = z + stp
            ENDIF
         ENDDO inner
      ENDDO main
      IF ( iequ<2 ) iequ = 2
      me%Sp(1) = p(1,iequ-1)
      me%Sp(2) = p(2,iequ-1)
      me%Sp(3) = p(3,iequ-1)
      IF ( oradik>=1.0e-15_wp ) fi = fi + stp/0.75_wp*oterm*oradik/(oradik-radik)
      !
      !-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
      !-- because 1E-38 is the minimal allowable arg. for ALOG in our envir.
      !-- D. Bilitza, Nov 87.
      !
      fi = 0.5_wp*abs(fi)/sqrt(B0) + 1.0e-12_wp
      !*****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR.
      !
      !-- Correct dipole moment is used here. D. Bilitza, Nov 87.
      !
      dimob0 = Dimo/B0
      arg1 = log(fi)
      arg2 = log(dimob0)
!       arg = FI*FI*FI/DIMOB0
!       if(abs(arg)>88.0_wp) arg=88.0_wp
      xx = 3*arg1 - arg2
      IF ( xx>23.0_wp ) THEN
         gg = xx - 3.0460681_wp
      ELSEIF ( xx>11.7_wp ) THEN
         gg = (((((2.8212095E-8_wp*xx-3.8049276E-6_wp)*xx+&
                   2.170224E-4_wp)*xx-6.7310339E-3_wp)*xx+&
                   1.2038224E-1_wp)*xx-1.8461796E-1_wp)*xx + 2.0007187_wp
      ELSEIF ( xx>+3.0_wp ) THEN
         gg = ((((((((6.3271665E-10_wp*xx-3.958306E-8_wp)*xx+&
                      9.9766148E-07_wp)*xx-1.2531932E-5_wp)*xx+&
                      7.9451313E-5_wp)*xx-3.2077032E-4_wp)*xx+&
                      2.1680398E-3_wp)*xx+1.2817956E-2_wp)*xx+&
                      4.3510529E-1_wp)*xx + 6.222355E-1_wp
      ELSEIF ( xx>-3.0_wp ) THEN
         gg = ((((((((2.6047023E-10_wp*xx+2.3028767E-9_wp)*xx-&
                      2.1997983E-8_wp)*xx-5.3977642E-7_wp)*xx-&
                      3.3408822E-6_wp)*xx+3.8379917E-5_wp)*xx+&
                      1.1784234E-3_wp)*xx+1.4492441E-2_wp)*xx+&
                      4.3352788E-1_wp)*xx + 6.228644E-1_wp
      ELSEIF ( xx>-22.0_wp ) THEN
         gg = ((((((((-8.1537735E-14_wp*xx+8.3232531E-13_wp)*xx+&
                       1.0066362E-9_wp)*xx+8.1048663E-8_wp)*xx+&
                       3.2916354E-6_wp)*xx+8.2711096E-5_wp)*xx+&
                       1.3714667E-3_wp)*xx+1.5017245E-2_wp)*xx+&
                       4.3432642E-1_wp)*xx + 6.2337691E-1_wp
      ELSE
         gg = 3.33338E-1_wp*xx + 3.0062102E-1_wp
      ENDIF
      Fl = exp(log((1.0_wp+exp(gg))*dimob0)/3.0_wp)
      RETURN
   END subroutine spag_block_1

END subroutine shellg

!*****************************************************************************************
!>
!  subroutine used for field line tracing in [[shellg]]
!  calls entry point [[feldi]] in geomagnetic field subroutine [[feldg]]

subroutine stoer(me,p,bq,r)

   class(shellig_type),intent(inout) :: me
   real(wp),dimension(7),intent(inout) :: p
   real(wp),intent(out) :: bq
   real(wp),intent(out) :: r

   real(wp) :: dr , dsq , dx , dxm , dy , dym , dz , &
               dzm , fli , q , rq , wr , xm , ym , zm

!*****XM,YM,ZM ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES
   zm = P(3)
   fli = P(1)*P(1) + P(2)*P(2) + 1.0e-15_wp
   R = 0.5_wp*(fli+sqrt(fli*fli+(zm+zm)**2))
   rq = R*R
   wr = sqrt(R)
   xm = P(1)*wr
   ym = P(2)*wr
!*****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM
   me%Xi(1) = xm*u(1,1) + ym*u(1,2) + zm*u(1,3)
   me%Xi(2) = xm*u(2,1) + ym*u(2,2) + zm*u(2,3)
   me%Xi(3) = xm*u(3,1) + zm*u(3,3)
!*****COMPUTE DERIVATIVES
! Changed from CALL FELDI(XI,H); XI, H are in COMMON block; results
! are the same; dkb Feb 1998.
! JW : feb 2024 : xi, h now class variables.
   CALL feldi(me)
   q = me%H(1)/rq
   dx = me%H(3) + me%H(3) + q*me%Xi(1)
   dy = me%H(4) + me%H(4) + q*me%Xi(2)
   dz = me%H(2) + me%H(2) + q*me%Xi(3)
!*****TRANSFORM BACK TO GEOMAGNETIC CO-ORDINATE SYSTEM
   dxm = u(1,1)*dx + u(2,1)*dy + u(3,1)*dz
   dym = u(1,2)*dx + u(2,2)*dy
   dzm = u(1,3)*dx + u(2,3)*dy + u(3,3)*dz
   dr = (xm*dxm+ym*dym+zm*dzm)/R
!*****FORM SLOWLY VARYING EXPRESSIONS
   P(4) = (wr*dxm-0.5_wp*P(1)*dr)/(R*dzm)
   P(5) = (wr*dym-0.5_wp*P(2)*dr)/(R*dzm)
   dsq = rq*(dxm*dxm+dym*dym+dzm*dzm)
   Bq = dsq*rq*rq
   P(6) = sqrt(dsq/(rq+3.0_wp*zm*zm))
   P(7) = P(6)*(rq+zm*zm)/(rq*dzm)
END subroutine stoer

!*****************************************************************************************
!>
!  Calculates earth magnetic field from spherical harmonics model
!
!### Reference
! ref: g. kluge, european space operations centre, internal note 61,
!      1970.
!
!### History
!  * changes (d. bilitza, nov 87):
!   - field coefficients in binary data files instead of block data
!   - calculates dipol moment

subroutine feldg(me,glat,glon,alt,bnorth,beast,bdown,babs)

   class(shellig_type),intent(inout) :: me
   real(wp),intent(in) :: glat  !! geodetic latitude in degrees (north)
   real(wp),intent(in) :: glon  !! geodetic longitude in degrees (east)
   real(wp),intent(in) :: alt   !! altitude in km above sea level
   real(wp),intent(out) :: bnorth, beast, bdown !! components of the field with respect
                                                !! to the local geodetic coordinate system, with axis
                                                !! pointing in the tangential plane to the north, east
                                                !! and downward.
   real(wp),intent(out) :: Babs !! magnetic field strength in gauss

   real(wp) :: b(3) , brho , bxxx , &
               byyy , bzzz , cp , ct , d , f , rho , &
               rlat , rlon , rq , s , sp , st , t , v(3) , x , xxx , &
               y , yyy , z , zzz
   integer :: i , ih , ihmax , il , imax , is , k , last , m

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

   !*****entry point  feldc  to be used with cartesian co-ordinates
   !        v(3)  cartesian coordinates in earth radii (6371.2 km)
   !            x-axis pointing to equator at 0 longitude
   !            y-axis pointing to equator at 90 long.
   !            z-axis pointing to north pole
         entry feldc(me,v,b)
         is=2
         xxx=v(1)
         yyy=v(2)
         zzz=v(3)

   10    rq=1.0_wp/(xxx*xxx+yyy*yyy+zzz*zzz)
         me%xi(1)=xxx*rq
         me%xi(2)=yyy*rq
         me%xi(3)=zzz*rq
         goto 20

   !*****entry point  feldi  used for l computation
         entry feldi(me)
         is=3
   20    ihmax=me%nmax*me%nmax+1
         last=ihmax+me%nmax+me%nmax
         imax=me%nmax+me%nmax-1
         do i=ihmax,last
             me%h(i)=me%g(i)
         end do
         do k=1,3,2
            i=imax
            ih=ihmax
            do
               il=ih-i
               f=2.0_wp/real(i-k+2, wp)
               x=me%xi(1)*f
               y=me%xi(2)*f
               z=me%xi(3)*(f+f)
               i=i-2
               if ((i-1)>=0) then
                  if ((i-1)>0) then
                     do m=3,i,2
                        me%h(il+m+1)=me%g(il+m+1)+z*me%h(ih+m+1)+x*(me%h(ih+m+3)-&
                                     me%h(ih+m-1))-y*(me%h(ih+m+2)+me%h(ih+m-2))
                        me%h(il+m)=me%g(il+m)+z*me%h(ih+m)+x*(me%h(ih+m+2)-&
                                   me%h(ih+m-2))+y*(me%h(ih+m+3)+me%h(ih+m-1))
                     end do
                  end if
                  me%h(il+2)=me%g(il+2)+z*me%h(ih+2)+x*me%h(ih+4)-y*(me%h(ih+3)+me%h(ih))
                  me%h(il+1)=me%g(il+1)+z*me%h(ih+1)+y*me%h(ih+4)+x*(me%h(ih+3)-me%h(ih))
               end if
               me%h(il)=me%g(il)+z*me%h(ih)+2.0_wp*(x*me%h(ih+1)+y*me%h(ih+2))
               ih=il
               if (i<k) exit
            end do
         end do

         if (is==3) return

         s=0.5_wp*me%h(1)+2.0_wp*(me%h(2)*me%xi(3)+me%h(3)*me%xi(1)+me%h(4)*me%xi(2))
         t=(rq+rq)*sqrt(rq)
         bxxx=t*(me%h(3)-s*xxx)
         byyy=t*(me%h(4)-s*yyy)
         bzzz=t*(me%h(2)-s*zzz)
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

!*****************************************************************************************
!>
!  Determines coefficients and dipol moment from IGRF models
!
!### Author
!  * D. BILITZA, NSSDC, GSFC, CODE 633, GREENBELT, MD 20771,
!    (301) 286-9536 NOV 1987.
!
!### History
!  * corrected for 2000 update - dkb- 5/31/2000
!  * updated to IGRF-2000 version -dkb- 5/31/2000
!  * updated to IGRF-2005 version -dkb- 3/24/2000

   subroutine feldcof(me,year,dimo)

   class(shellig_type),intent(inout) :: me
   real(wp),intent(in) :: year !! decimal year for which geomagnetic field is to
                               !! be calculated (e.g.:1995.5 for day 185 of 1995)
   real(wp),intent(out) :: dimo !! geomagnetic dipol moment in gauss (normalized
                                !! to earth's radius) at the time (year)

   real(wp) :: dte1 , dte2 , erad , gha(144) , sqrt2
   integer :: i , ier , j , l , m , n , nmax1 , nmax2, iyea
   character(len=filename_len) :: fil2
   real(wp) :: x , f0 , f !! these were double precision in original
                          !! code while everything else was single precision

   ! changed to conform with IGRF 45-95, also FILMOD, DTEMOD arrays +1
   character(len=filename_len),dimension(17),parameter :: filmod = [&
         'dgrf1945.dat ' , 'dgrf1950.dat ' , 'dgrf1955.dat ' , 'dgrf1960.dat ' , &
         'dgrf1965.dat ' , 'dgrf1970.dat ' , 'dgrf1975.dat ' , 'dgrf1980.dat ' , &
         'dgrf1985.dat ' , 'dgrf1990.dat ' , 'dgrf1995.dat ' , 'dgrf2000.dat ' , &
         'dgrf2005.dat ' , 'dgrf2010.dat ' , 'dgrf2015.dat ' , 'igrf2020.dat ' , &
         'igrf2020s.dat']
   real(wp),dimension(17),parameter :: dtemod = [1945.0_wp , 1950.0_wp , 1955.0_wp , &
                                                 1960.0_wp , 1965.0_wp , 1970.0_wp , &
                                                 1975.0_wp , 1980.0_wp , 1985.0_wp , &
                                                 1990.0_wp , 1995.0_wp , 2000.0_wp , &
                                                 2005.0_wp , 2010.0_wp , 2015.0_wp , &
                                                 2020.0_wp , 2025.0_wp]
   integer,parameter :: numye = size(dtemod)-1 ! number of 5-year priods represented by IGRF
   integer,parameter :: is = 0 !! * is=0 for schmidt normalization
                               !! * is=1 gauss normalization

   logical :: read_file

   !-- determine igrf-years for input-year
   me%time = year
   iyea = int(year/5.0_wp)*5
   read_file = iyea /= me%iyea  ! if we have to read the file
   me%iyea = iyea
   l = (me%iyea-1945)/5 + 1
   if ( l<1 ) l = 1
   if ( l>numye ) l = numye
   dte1 = dtemod(l)
   me%name = filmod(l)
   dte2 = dtemod(l+1)
   fil2 = filmod(l+1)
   if (read_file) then
      ! get igrf coefficients for the boundary years
      ! [if they have not ready been loaded]
      call me%getshc(me%name,nmax1,erad,me%g,ier)
      if ( ier/=0 ) error stop 'error reading file: '//trim(me%name)
      call me%getshc(fil2,nmax2,erad,me%gh2,ier)
      if ( ier/=0 ) error stop 'error reading file: '//trim(fil2)
   end if
   !-- determine igrf coefficients for year
   if ( l<=numye-1 ) then
      call me%intershc(year,dte1,nmax1,me%g,dte2,nmax2,me%gh2,me%nmax,gha)
   else
      call me%extrashc(year,dte1,nmax1,me%g,nmax2,me%gh2,me%nmax,gha)
   endif
   !-- determine magnetic dipol moment and coeffiecients g
   f0 = 0.0_wp
   do j = 1 , 3
      f = gha(j)*1.0e-5_wp
      f0 = f0 + f*f
   enddo
   dimo = sqrt(f0)

   me%g(1) = 0.0_wp
   i = 2
   f0 = 1.0e-5_wp
   if ( is==0 ) f0 = -f0
   sqrt2 = sqrt(2.0_wp)

   do n = 1 , me%nmax
      x = n
      f0 = f0*x*x/(4.0_wp*x-2.0_wp)
      if ( is==0 ) f0 = f0*(2.0_wp*x-1.0_wp)/x
      f = f0*0.5_wp
      if ( is==0 ) f = f*sqrt2
      me%g(i) = gha(i-1)*f0
      i = i + 1
      do m = 1 , n
         f = f*(x+m)/(x-m+1.0_wp)
         if ( is==0 ) f = f*sqrt((x-m+1.0_wp)/(x+m))
         me%g(i) = gha(i-1)*f
         me%g(i+1) = gha(i)*f
         i = i + 2
      enddo
   enddo

end subroutine feldcof

!*****************************************************************************************
!>
!  Reads spherical harmonic coefficients from the specified
!  file into an array.
!
!### Author
!  * Version 1.01, A. Zunde, USGS, MS 964,
!    Box 25046 Federal Center, Denver, CO  80225

subroutine getshc(me,Fspec,Nmax,Erad,Gh,Ier)

   class(shellig_type),intent(inout) :: me
   character(len=*),intent(in) :: Fspec !! File specification
   integer,intent(out) :: Nmax !! Maximum degree and order of model
   real(wp),intent(out) :: Erad !! Earth's radius associated with the spherical
                                !! harmonic coefficients, in the same units as
                                !! elevation
   real(wp),dimension(*),intent(out) :: Gh !! Schmidt quasi-normal internal spherical
                                           !! harmonic coefficients
   integer,intent(out) :: Ier !! Error number:
                              !!
                              !!  * 0, no error
                              !!  * -2, records out of order
                              !!  * FORTRAN run-time error number

   integer :: iu !! logical unit number
   real(wp) :: g , h
   integer :: i , m , mm , n , nn

   read_file : block
      ! ---------------------------------------------------------------
      !  Open coefficient file. Read past first header record.
      !  Read degree and order of model and Earth's radius.
      ! ---------------------------------------------------------------
      OPEN (newunit=Iu,FILE=Fspec,STATUS='OLD',IOSTAT=Ier)
      if (Ier/=0) then
         write(*,*) 'Error opening file: '//trim(fspec)
         exit read_file
      end if
      READ (Iu,*,IOSTAT=Ier)
      if (Ier/=0) exit read_file
      READ (Iu,*,IOSTAT=Ier) Nmax , Erad
      if (Ier/=0) exit read_file

      ! ---------------------------------------------------------------
      !  Read the coefficient file, arranged as follows:
      !
      !          N     M     G     H
      !          ----------------------
      !            /   1     0    GH(1)  -
      !           /  1     1    GH(2) GH(3)
      !          /  2     0    GH(4)  -
      !         /  2     1    GH(5) GH(6)
      !      NMAX*(NMAX+3)/2   /  2     2    GH(7) GH(8)
      !         records    \  3     0    GH(9)  -
      !         \      .     .     .     .
      !          \  .     .     .     .
      !      NMAX*(NMAX+2)     \  .     .     .     .
      !      elements in GH      \  NMAX  NMAX   .     .
      !
      !  N and M are, respectively, the degree and order of the
      !  coefficient.
      ! ---------------------------------------------------------------
      i = 0
      main: DO nn = 1 , Nmax
         DO mm = 0 , nn
            READ (Iu,*,IOSTAT=Ier) n , m , g , h
            if (Ier/=0) exit main
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

   end block read_file

   CLOSE (Iu)

END subroutine getshc

!*****************************************************************************************
!>
!  Interpolates linearly, in time, between two spherical
!  harmonic models.
!
!  The coefficients (GH) of the resulting model, at date
!  DATE, are computed by linearly interpolating between the
!  coefficients of the earlier model (GH1), at date DTE1,
!  and those of the later model (GH2), at date DTE2. If one
!  model is smaller than the other, the interpolation is
!  performed with the missing coefficients assumed to be 0.
!
!### Author
!  * Version 1.01, A. Zunde
!    USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225

subroutine intershc(me,date,dte1,nmax1,gh1,dte2,nmax2,gh2,nmax,gh)

   class(shellig_type),intent(inout) :: me
   real(wp),intent(in) :: date !! Date of resulting model (in decimal year)
   real(wp),intent(in) :: dte1 !! Date of earlier model
   integer,intent(in) :: nmax1 !! Maximum degree and order of earlier model
   real(wp),intent(in) :: gh1(*) !! Schmidt quasi-normal internal spherical harmonic coefficients of earlier model
   real(wp),intent(in) :: dte2 !! Date of later model
   integer,intent(in) :: nmax2 !! Maximum degree and order of later model
   real(wp),intent(in) :: gh2(*) !! Schmidt quasi-normal internal spherical harmonic coefficients of later model
   real(wp),intent(out) :: gh(*) !! Coefficients of resulting model
   integer,intent(out) :: nmax !! Maximum degree and order of resulting model

   real(wp) :: factor
   integer :: i , k , l

   factor = (date-dte1)/(dte2-dte1)

   if ( nmax1==nmax2 ) then
      k = nmax1*(nmax1+2)
      nmax = nmax1
   elseif ( nmax1>nmax2 ) then
      k = nmax2*(nmax2+2)
      l = nmax1*(nmax1+2)
      do i = k + 1 , l
         gh(i) = gh1(i) + factor*(-gh1(i))
      enddo
      nmax = nmax1
   else
      k = nmax1*(nmax1+2)
      l = nmax2*(nmax2+2)
      do i = k + 1 , l
         gh(i) = factor*gh2(i)
      enddo
      nmax = nmax2
   endif

   do i = 1 , k
      gh(i) = gh1(i) + factor*(gh2(i)-gh1(i))
   enddo

end subroutine intershc

!*****************************************************************************************
!>
!  Extrapolates linearly a spherical harmonic model with a
!  rate-of-change model.
!
!  The coefficients (GH) of the resulting model, at date
!  DATE, are computed by linearly extrapolating the coef-
!  ficients of the base model (GH1), at date DTE1, using
!  those of the rate-of-change model (GH2), at date DTE2. If
!  one model is smaller than the other, the extrapolation is
!  performed with the missing coefficients assumed to be 0.
!
!### Author
!  * Version 1.01, A. Zunde
!    USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225

subroutine extrashc(me,date,dte1,nmax1,gh1,nmax2,gh2,nmax,gh)

   class(shellig_type),intent(inout) :: me
   real(wp),intent(in) :: date   !! Date of resulting model (in decimal year)
   real(wp),intent(in) :: dte1   !! Date of base model
   integer,intent(in)  :: nmax1  !! Maximum degree and order of base model
   real(wp),intent(in) :: gh1(*) !! Schmidt quasi-normal internal spherical harmonic coefficients of base model
   integer,intent(in)  :: nmax2  !! Maximum degree and order of rate-of-change model
   real(wp),intent(in) :: gh2(*) !! Schmidt quasi-normal internal spherical harmonic coefficients of rate-of-change model
   real(wp),intent(out) :: gh(*) !! Coefficients of resulting model
   integer,intent(out) :: nmax   !! Maximum degree and order of resulting model

   real(wp) :: factor
   integer :: i , k , l

   factor = (date-dte1)

   if ( nmax1==nmax2 ) then
      k = nmax1*(nmax1+2)
      nmax = nmax1
   elseif ( nmax1>nmax2 ) then
      k = nmax2*(nmax2+2)
      l = nmax1*(nmax1+2)
      do i = k + 1 , l
         gh(i) = gh1(i)
      enddo
      nmax = nmax1
   else
      k = nmax1*(nmax1+2)
      l = nmax2*(nmax2+2)
      do i = k + 1 , l
         gh(i) = factor*gh2(i)
      enddo
      nmax = nmax2
   endif

   do i = 1 , k
      gh(i) = gh1(i) + factor*gh2(i)
   enddo

end subroutine extrashc

end module SHELLIG_module