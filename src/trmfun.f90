!*==trara1.f90 processed by SPAG 8.01MH 17:13 30 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
! TRMFUN.FOR	1987
!
!********************************************************************
!*************** SUBROUTINES, FUNCTIONS *****************************
!********************************************************************
!******************* TRARA1, TRARA2 *********************************
!********************************************************************
!
SUBROUTINE trara1(Descr,Map,Fl,Bb0,E,F,N)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL Bb0 , E , e0 , e1 , e2 , escale , F , f0 , f1 , f2 , Fistep , Fl , fscale , trara2 , xnl
   INTEGER i0 , i1 , i2 , i3 , ie , l3 , Map , N , nb , nl
!*** End of declarations inserted by SPAG
!***********************************************************************
!*** TRARA1 FINDS PARTICLE FLUXES FOR GIVEN ENERGIES, MAGNETIC FIELD ***
!*** STRENGTH AND L-VALUE. FUNCTION TRARA2 IS USED TO INTERPOLATE IN ***
!*** B-L-SPACE.                                                      ***
!***   INPUT: DESCR(8)   HEADER OF SPECIFIED TRAPPED RADITION MODEL  ***
!***          MAP(...)   MAP OF TRAPPED RADITION MODEL               ***
!***                     (DESCR AND MAP ARE EXPLAINED AT THE BEGIN   ***
!***                     OF THE MAIN PROGRAM MODEL)                  ***
!***          N          NUMBER OF ENERGIES                          ***
!***          E(N)       ARRAY OF ENERGIES IN MEV                    ***
!***          FL         L-VALUE                                     ***
!***          BB0        =B/B0  MAGNETIC FIELD STRENGTH NORMALIZED   ***
!***                     TO FIELD STRENGTH AT MAGNETIC EQUATOR       ***
!***  OUTPUT: F(N)       DECADIC LOGARITHM OF INTEGRAL FLUXES IN     ***
!***                     PARTICLES/(CM*CM*SEC)                       ***
!***********************************************************************
   LOGICAL s0 , s1 , s2
   DIMENSION E(N) , F(N) , Map(*)
   INTEGER Descr(8)
   COMMON /tra2  / Fistep
   DATA f1 , f2/1.001 , 1.002/
!
   Fistep = Descr(7)/Descr(2)
   escale = Descr(4)
   fscale = Descr(7)
   xnl = amin1(15.6,abs(Fl))
   nl = xnl*Descr(5)
   IF ( Bb0<1. ) Bb0 = 1.
   nb = (Bb0-1.)*Descr(6)
!
! I2 IS THE NUMBER OF ELEMENTS IN THE FLUX MAP FOR THE FIRST ENERGY.
! I3 IS THE INDEX OF THE LAST ELEMENT OF THE SECOND ENERGY MAP.
! L3 IS THE LENGTH OF THE MAP FOR THE THIRD ENERGY.
! E1 IS THE ENERGY OF THE FIRST ENERGY MAP (UNSCALED)
! E2 IS THE ENERGY OF THE SECOND ENERGY MAP (UNSCALED)
!
   i1 = 0
   i2 = Map(1)
   i3 = i2 + Map(i2+1)
   l3 = Map(i3+1)
   e1 = Map(i1+2)/escale
   e2 = Map(i2+2)/escale
!
! S0, S1, S2 ARE LOGICAL VARIABLES WHICH INDICATE WHETHER THE FLUX FOR
! A PARTICULAR E, B, L POINT HAS ALREADY BEEN FOUND IN A PREVIOUS CALL
! TO FUNCTION TRARA2. IF NOT, S.. =.TRUE.
!
   s1 = .TRUE.
   s2 = .TRUE.
!
!			ENERGY LOOP
!
   DO ie = 1 , N
!
! FOR EACH ENERGY E(I) FIND THE SUCCESSIVE ENERGIES E0,E1,E2 IN
! MODEL MAP, WHICH OBEY  E0 < E1 < E(I) < E2 .
!
      DO WHILE ( (E(ie)>e2) .AND. (l3/=0) )
         i0 = i1
         i1 = i2
         i2 = i3
         i3 = i3 + l3
         l3 = Map(i3+1)
         e0 = e1
         e1 = e2
         e2 = Map(i2+2)/escale
         s0 = s1
         s1 = s2
         s2 = .TRUE.
         f0 = f1
         f1 = f2
      ENDDO
!
! CALL TRARA2 TO INTERPOLATE THE FLUX-MAPS FOR E1,E2 IN L-B/B0-
! SPACE TO FIND FLUXES F1,F2 [IF THEY HAVE NOT ALREADY BEEN
! CALCULATED FOR A PREVIOUS E(I)].
!
      IF ( s1 ) f1 = trara2(Map(i1+3),nl,nb)/fscale
      IF ( s2 ) f2 = trara2(Map(i2+3),nl,nb)/fscale
      s1 = .FALSE.
      s2 = .FALSE.
!
! FINALLY, INTERPOLATE IN ENERGY.
!
      F(ie) = f1 + (f2-f1)*(E(ie)-e1)/(e2-e1)
      IF ( f2<=0.0 ) THEN
         IF ( i1/=0 ) THEN
!
! --------- SPECIAL INTERPOLATION ---------------------------------
! IF THE FLUX FOR THE SECOND ENERGY CANNOT BE FOUND (I.E. F2=0.0),
! AND THE ZEROTH ENERGY MAP HAS BEEN DEFINED (I.E. I1 NOT EQUAL 0),
! THEN INTERPOLATE USING THE FLUX MAPS FOR THE ZEROTH AND FIRST
! ENERGY AND CHOOSE THE MINIMUM OF THIS INTERPOLATIONS AND THE
! INTERPOLATION THAT WAS DONE WITH F2=0.
!
            IF ( s0 ) f0 = trara2(Map(i0+3),nl,nb)/fscale
            s0 = .FALSE.
            F(ie) = amin1(F(ie),f0+(f1-f0)*(E(ie)-e0)/(e1-e0))
         ENDIF
      ENDIF
!
! THE LOGARITHMIC FLUX IS ALWAYS KEPT GREATER OR EQUAL ZERO.
!
      F(ie) = amax1(F(ie),0.)
   ENDDO
END SUBROUTINE trara1
!*==trara2.f90 processed by SPAG 8.01MH 17:13 30 Jan 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
FUNCTION trara2(Map,Il,Ib)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   REAL dfl , fincr1 , fincr2 , Fistep , fkb , fkb1 , fkb2 , fkbj1 , fkbj2 , fkbm , fll1 , fll2 , flog , flog1 , flog2 , flogm ,   &
      & fnb , fnl , sl1 , sl2
   REAL trara2
   INTEGER i1 , i2 , Ib , Il , itime , j1 , j2 , kt , l1 , l2 , Map
!*** End of declarations inserted by SPAG
!*****************************************************************
!***  TRARA2 INTERPOLATES LINEARLY IN L-B/B0-MAP TO OBTAIN     ***
!***  THE LOGARITHM OF INTEGRAL FLUX AT GIVEN L AND B/B0.      ***
!***    INPUT: MAP(..) IS SUB-MAP (FOR SPECIFIC ENERGY) OF     ***
!***                   TRAPPED RADIATION MODEL MAP             ***
!***           IL      SCALED L-VALUE                          ***
!***           IB      SCALED B/B0-1                           ***
!***   OUTPUT: TRARA2  SCALED LOGARITHM OF PARTICLE FLUX       ***
!*****************************************************************
!***  SEE MAIN PROGRAM 'MODEL' FOR EXPLANATION OF MAP FORMAT   ***
!***  SCALING FACTORS.                                         ***
!***  THE STEPSIZE FOR THE PARAMETERIZATION OF THE LOGARITHM   ***
!***  OF FLUX IS OBTAINED FROM 'COMMON/TRA2/'.                 ***
!*****************************************************************
   DIMENSION Map(*)
   COMMON /tra2  / Fistep
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
         fnl = Il
         fnb = Ib
         itime = 0
         i2 = 0
         SPAG_Loop_1_1: DO
!
! FIND CONSECUTIVE SUB-SUB-MAPS FOR SCALED L-VALUES LS1,LS2,
! WITH IL LESS OR EQUAL LS2.  L1,L2 ARE LENGTHS OF SUB-SUB-MAPS.
! I1,I2 ARE INDECES OF FIRST ELEMENTS MINUS 1.
!
            l2 = Map(i2+1)
            IF ( Map(i2+2)<=Il ) THEN
               i1 = i2
               l1 = l2
               i2 = i2 + l2
!
! IF SUB-SUB-MAPS ARE EMPTY, I. E. LENGTH LESS 4, THAN TRARA2=0
!
            ELSEIF ( (l1<4) .AND. (l2<4) ) THEN
               trara2 = 0.
               RETURN
            ELSE
!
! IF FLOG2 LESS FLOG1, THAN LS2 FIRST MAP AND LS1 SECOND MAP
!
               IF ( Map(i2+3)<=Map(i1+3) ) EXIT SPAG_Loop_1_1
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO SPAG_Loop_1_1
         spag_nextblock_1 = 2
      CASE (2)
         kt = i1
         i1 = i2
         i2 = kt
         kt = l1
         l1 = l2
         l2 = kt
         spag_nextblock_1 = 3
      CASE (3)
!
! DETERMINE INTERPOLATE IN SCALED L-VALUE
!
         fll1 = Map(i1+2)
         fll2 = Map(i2+2)
         dfl = (fnl-fll1)/(fll2-fll1)
         flog1 = Map(i1+3)
         flog2 = Map(i2+3)
         fkb1 = 0.
         fkb2 = 0.
         IF ( l1>=4 ) THEN
!
! B/B0 LOOP
!
            DO j2 = 4 , l2
               fincr2 = Map(i2+j2)
               IF ( fkb2+fincr2>fnb ) GOTO 10
               fkb2 = fkb2 + fincr2
               flog2 = flog2 - Fistep
            ENDDO
            itime = itime + 1
            IF ( itime==1 ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            trara2 = 0.
            RETURN
 10         IF ( itime/=1 ) THEN
               IF ( j2==4 ) THEN
                  spag_nextblock_1 = 4
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               sl2 = flog2/fkb2
               DO j1 = 4 , l1
                  fincr1 = Map(i1+j1)
                  fkb1 = fkb1 + fincr1
                  flog1 = flog1 - Fistep
                  fkbj1 = ((flog1/Fistep)*fincr1+fkb1)/((fincr1/Fistep)*sl2+1.)
                  IF ( fkbj1<=fkb1 ) GOTO 15
               ENDDO
               IF ( fkbj1<=fkb2 ) THEN
                  trara2 = 0.
                  RETURN
               ENDIF
 15            IF ( fkbj1<=fkb2 ) THEN
                  fkbm = fkbj1 + (fkb2-fkbj1)*dfl
                  flogm = fkbm*sl2
                  flog2 = flog2 - Fistep
                  fkb2 = fkb2 + fincr2
                  sl1 = flog1/fkb1
                  sl2 = flog2/fkb2
                  spag_nextblock_1 = 5
                  CYCLE SPAG_DispatchLoop_1
               ELSE
                  fkb1 = 0.
               ENDIF
            ENDIF
            fkb2 = 0.
         ENDIF
         j2 = 4
         fincr2 = Map(i2+j2)
         flog2 = Map(i2+3)
         flog1 = Map(i1+3)
         spag_nextblock_1 = 4
      CASE (4)
         flogm = flog1 + (flog2-flog1)*dfl
         fkbm = 0.
         fkb2 = fkb2 + fincr2
         flog2 = flog2 - Fistep
         sl2 = flog2/fkb2
         IF ( l1<4 ) THEN
            fincr1 = 0.
            sl1 = -900000.
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ELSE
            j1 = 4
            fincr1 = Map(i1+j1)
            fkb1 = fkb1 + fincr1
            flog1 = flog1 - Fistep
            sl1 = flog1/fkb1
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
         DO WHILE ( sl1>=sl2 )
            fkbj2 = ((flog2/Fistep)*fincr2+fkb2)/((fincr2/Fistep)*sl1+1.)
            fkb = fkb1 + (fkbj2-fkb1)*dfl
            flog = fkb*sl1
            IF ( fkb>=fnb ) THEN
               spag_nextblock_1 = 7
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            fkbm = fkb
            flogm = flog
            IF ( j1>=l1 ) THEN
               trara2 = 0.
               RETURN
            ELSE
               j1 = j1 + 1
               fincr1 = Map(i1+j1)
               flog1 = flog1 - Fistep
               fkb1 = fkb1 + fincr1
               sl1 = flog1/fkb1
            ENDIF
         ENDDO
         spag_nextblock_1 = 6
      CASE (6)
         fkbj1 = ((flog1/Fistep)*fincr1+fkb1)/((fincr1/Fistep)*sl2+1.)
         fkb = fkbj1 + (fkb2-fkbj1)*dfl
         flog = fkb*sl2
         IF ( fkb<fnb ) THEN
            fkbm = fkb
            flogm = flog
            IF ( j2>=l2 ) THEN
               trara2 = 0.
               RETURN
            ELSE
               j2 = j2 + 1
               fincr2 = Map(i2+j2)
               flog2 = flog2 - Fistep
               fkb2 = fkb2 + fincr2
               sl2 = flog2/fkb2
               spag_nextblock_1 = 5
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
         spag_nextblock_1 = 7
      CASE (7)
         IF ( fkb<fkbm+1.E-10 ) THEN
            trara2 = 0.
         ELSE
            trara2 = flogm + (flog-flogm)*((fnb-fkbm)/(fkb-fkbm))
            trara2 = amax1(trara2,0.)
            RETURN
         ENDIF
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END FUNCTION trara2
