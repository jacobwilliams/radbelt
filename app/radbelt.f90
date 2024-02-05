PROGRAM spag_program_1

   use radbelt_module
   use radbelt_kinds_module

   IMPLICIT NONE
   REAL(wp) af , bb0 , bbeg , bend , blv , bstep , df , e , ebeg , eda , &
            ediff , eend , ei , estep , fl , flux , vbeg , vend , vstep , &
            xl
   REAL(wp) xmax , xmin , xx
   INTEGER i , ib , ibl , ibltab , ie , iei , ihead , il , inde , itest , &
           itt , iuaeap , iuout , j , jagnr , jpara , jtab , k , l ,&
           map
   INTEGER maxi , mini , monito , mtype , n , nb , ne , nl , nmap , nn
! RADBELT.FOR   SEPTEMBER 88
!
!***********************************************************************
!*                                                                     *
!*                    TRAPPED RADIATION MODELS                         *
!*                         Program RADBELT                             *
!*                                                                     *
!***  Dieter Bilitza  ******************************  25 March 1988  ***
!***********************************************************************
!***************  Goddard Space Flight Center, code 633  ***************
!**  National Space Science Data Center, Greenbelt, MD 20771, U.S.A.  **
!***********************************************************************
!**********************  NSSDC-ID  PT-14A  *****************************
!***********************************************************************
!***   This program  gives an example of how to use NSSDC's Trapped  ***
!***   Radiation Models. It determines omnidirectional integral and  ***
!***   differential electron or proton fluxes for given energies,    ***
!***   L-values, magnetic field strengths and map-type.              ***
!***   The program will ask you for:                                 ***
!***      NE     number of energies you want displayed               ***
!***      E(1),... E(NE)   energies   or   EBEGIN,EEND,ESTEP  begin, ***
!***             end and stepsize of energy range                    ***
!***      NL     number of L-values                                  ***
!***      L(1),... L(NL)   L-values   or   LBEGIN,LEND,LSTEP         ***
!***      NB     number of B/B0-values                               ***
!***      B/B0(1),... B/B0(NL)   B/B0-values   or  BBEGIN,BEND,BSTEP ***
!***      MTYPE  map type:  1 AP8MAX   2 AP8MIN  3 AE4MAX  4 AE4MIN  ***
!***                        5 AEI7HI   6 AEI7LO  7 AE8MAX  8 AE8MIN  ***
!***      JTAB     output options: integral or differential fluxes     ***
!***             versus L or B/B0                                    ***
!***   The program interpolates the NSSDC model map in B/B0, L       ***
!***   (function TRARA2), and energy (subroutine TRARA1).            ***
!***   The program opens the map as binary data file (e.g.           ***
!***   AE8MAX.BIN) on unit 15. Make sure that the map you want to    ***
!***   use is available under the correct name in your account       ***
!***   ( see MTYPE for available models and names ).                 ***
!***********************************************************************
!***********************************************************************
!**************************  NOMENCLATURE  *****************************
!********  omnidirectional flux is flux per unit time and      *********
!********      unit sphere surface.                            *********
!********  integral flux for energy E is flux per unit time    *********
!********      and unit surface of particles with energies     *********
!********      greater than E.                                 *********
!********  differential flux for energy E and energy range DE  *********
!********      is average flux per unit time, unit surface and *********
!********      unit energy of particles with energies between  *********
!********      E-DE/2 and E+DE/2.                              *********
!***********************************************************************
!******************************  UNITS  ********************************
!********                 energies: MeV                        *********
!********          integral fluxes: particles/(cm*cm*sec)      *********
!********      differential fluxes: particles/(cm*cm*sec*MeV)  *********
!********  magnetic field strength: Gauss                      *********
!***********************************************************************
!***********************************************************************
!*************  DESCRIPTION OF MODEL DATA FILE FORMAT  *****************
!***********************************************************************
!***  THE FILE CONSISTS OF A HEADER ARRAY (IHEAD(8)) AND A MODEL MAP ***
!***  ARRAY (MAP(...)). ALL ELEMENTS ARE INTEGER.                    ***
!***                                                                 ***
!***  IHEAD(1)   MODEL MAP TYPE (SEE ABOVE)                          ***
!***       (2)   INCREMENTS PER DECADE OF LOGARITHMIC FLUX           ***
!***       (3)   EPOCH OF MODEL                                      ***
!***       (4)   SCALE FACTOR FOR ENERGY; E/MEV=E(MAP)/IHEAD(4)      ***
!***                    =6400 (AE-8),  =100 (AP-8) ***
!***       (5)   SCALE FACTOR FOR L-VALUE =2100 (AE-8), =2048 (AP-8) ***
!***       (6)   SCALE FACTOR FOR B/B0    =1024 (AE-8), =2048 (AP-8) ***
!***       (7)   SCALE FACTOR FOR LOGARITHM OF FLUXES =1024 (AE,AP-8)***
!***       (8)   NUMBER OF ELEMENTS IN MAP =13548 (AE8MAX),        ***
!***        =13168 (AE8MIN), =6509 (AP8MAX), =6688 (AP8MIN)    ***
!***                                                                 ***
!***  LAYOUT OF MAP:                                                 ***
!***      MAP CONSISTS OF SEVERAL VARIABLE-LENGTH SUB-MAPS, EACH     ***
!***      FOR A DIFFERENT ENERGY. EACH SUB-MAP CONSISTS OF SEVERAL   ***
!***      VARIABLE-LENGTH SUB-SUB-MAPS EACH FOR A DIFFERENT L-VALUE. ***
!***      EACH SUB-SUB-MAP CONTAINS THE CURVE LOG(F) [DECADIC        ***
!***      LOGARITHM OF OMNIDIRECTIONAL INTEGRAL PARTICLE FLUX]       ***
!***      VERSUS B/B0 [MAGNETIC FIELD STRENGTH NORMALIZED TO THE     ***
!***      EQUATORIAL VALUE]. THE CURVE IS PARAMETERIZED BY USING     ***
!***      EQUAL INCREMENTS IN LOG(F); THE NUMBER OF INCREMENTS       ***
!***      PER DECADE IS LISTED IN THE HEADER ARRAY [IHEAD(2)]:       ***
!***                                                                 ***
!***         I     B(I)/B(0)   (B(I)-B(I-1))/B(0)   LOG(F(I))        ***
!***       ----------------------------------------------------      ***
!***         0        1                 -              Y             ***
!***         1     B(1)/B(0)   (B(1)-B(0))/B(0)      Y-1/IHEAD(2)    ***
!***         2     B(2)/B(0)   (B(2)-B(1))/B(0)      Y-2/IHEAD(2)    ***
!***         .       ....            .....            ....           ***
!***                                                                 ***
!***      THE SUB-SUB-MAP CONTAINS THE EQUATORIAL FLUX LOGARITHM Y   ***
!***      AND THE B/B0-INCREMENTS (THIRD COLUMN) MULTIPLIED BY       ***
!***      THEIR CORRESPONDING SCALE VALUES ( IHEAD(7) AND (8) ).     ***
!***                                                                 ***
!***  MAP(1)  NUMBER OF ELEMENTS IN SUB-MAP                          ***
!***  MAP(2)  ENERGY FOR THIS SUB-MAP; MAP(2)=E/MEV*IHEAD(4)         ***
!***    MAP(3)  NUMBER OF ELEMENTS IN SUB-SUB-MAP                    ***
!***    MAP(4)  L-VALUE FOR THIS SUB-SUB-MAP; MAP(4)=L*IHEAD(5)      ***
!***    MAP(5)  LOGARITHM OF FLUX AT EQUATOR; MAP(5)=LOG(F0)*IHEAD(7)***
!***      MAP(6)  =(B1-B0)/B0; B1 IS THE MAGNETIC FIELD STRENGTH     ***
!***              THAT CORRESPONDS TO LOG(F1)=LOG(F0)-1/IHEAD(2)     ***
!***      MAP(7)  =(B2-B1)/B0; LOG(F2)=LOG(F1)-1/IHEAD(2)            ***
!***       ...              ....                                     ***
!***      MAP(L)  LAST ELEMENT IN SUB-SUB-MAP; L=MAP(3)+2            ***
!***    MAP(I)  NUMBER OF ELEMENTS IN NEXT SUB-SUB-MAP; I=L+1        ***
!***       ...              ....                                     ***
!***     ...                ....                                     ***
!***    MAP( )  NUMBER OF ELEMENTS IN LAST SUB-SUB-MAP               ***
!***       ...              ....                                     ***
!***      MAP(K)  LAST ELEMENT IN SUB-MAP; K=MAP(1)                  ***
!***  MAP(J)  NUMBER OF ELEMENTS IN NEXT SUB-MAP; J=MAP(1)+1         ***
!***     ...                ....                                     ***
!***       ...              ....                                     ***
!***   ...                  ....                                     ***
!***  MAP( )  NUMBER OF ELEMENTS IN LAST SUB-MAP                     ***
!***     ...                ....                                     ***
!***       ...              ....                                     ***
!***      MAP(M)  LAST ELEMENT OF MAP; M=IHEAD(8)                    ***
!***********************************************************************
!**                        ENERGY/MEV GRID                           ***
!**  AE-8:    0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
!**           3.5   4.0   4.5   5.0   5.5   6.0  6.5  7.0 (18 GRID P)***
!**  AE-5,6:  0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
!**           4.0   4.5   (12 GRID POINTS)                           ***
!**  AE-4:    0.04  0.1   0.3   0.5   1.0   2.0  2.5  3.0  3.5  4.0  ***
!**           4.1   4.25  4.35  4.5   4.65  4.85 (16 GRID POINTS)    ***
!**                                                                  ***
!**                         L-VALUE GRID                             ***
!**           BEGIN  STEP  END   STEP  END   STEP  END   STEP  END   ***
!**  AE-8:     1.2   0.05  1.5    0.1  2.0    0.2  2.4    0.1  3.0   ***
!**            3.0   0.2   3.4    0.1  3.6    0.2  4.4    0.1  4.6   ***
!**            4.6   0.2   5.0    0.5* 8.0    1.0 12.0 (43 GRID P.)  ***
!**                    * 6.6 INSTEAD OF 6.5                          ***
!**  AE-5,6:   1.2   0.05  1.5    0.1  2.0    0.2  2.8 (15 GRID P.)  ***
!**  AE-4:     2.8   0.2   4.0    0.5  6.0    0.6  6.6    0.4  7.0   ***
!**            7.0   1.0  11.0   (16 GRID POINTS)                    ***
!***********************************************************************
   INTEGER agnr , egnr
   DIMENSION xl(2,10) , e(7) , flux(7) , af(7,10,10) , ihead(8) , map(20000) , df(5,10,10) , eda(3)
   LOGICAL notbeg
   CHARACTER*1 ite(5,10,10)
   CHARACTER*4 bltex , lbtex
   CHARACTER*5 untex
   CHARACTER*6 name , mname(8)
   CHARACTER*9 particle
   CHARACTER*10 fname
   CHARACTER*12 idtex
   INTEGER :: spag_nextblock_1
   DATA mname/'AP8MAX' , 'AP8MIN' , 'AE4MAX' , 'AE4MIN' , 'AEI7HI' , 'AEI7LO' , 'AE8MAX' , 'AE8MIN'/
   DATA mtype , ne , ebeg , eend , estep , nl , vbeg , vend , vstep , nb , bbeg , bend , bstep , jagnr , inde , jtab/7 , 0 , 1. ,  &
      & 6. , 1. , 0 , 1. , 10. , 1. , 0 , 1 , 100 , 10 , 0 , 1 , 1/
   DATA e , xl , af , df/1227*0.0/ , eda/.05 , .1 , .2/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         DO i = 1 , 5
            DO k = 1 , 10
               DO l = 1 , 10
                  ite(i,k,l) = ' '
               ENDDO
            ENDDO
         ENDDO
!
!      I/O UNIT NUMBERS
!
         egnr = 5               ! INPUT
         monito = 6             ! MONITOR
         iuaeap = 15            ! MODEL COEFFICIENTS INPUT
         iuout = 16             ! OUTPUT (OUTPUT.AEP)
!
!      INPUT OF PARAMETERS
!
         notbeg = .FALSE.
!-----------1. window: introduction--------------------------------
         WRITE (monito,99001)
99001    FORMAT (1X/////8X,'**** AE-8, AP-8 TRAPPED PARTICLE MODELS',' ****'//4X,                                                  &
                &'This program allows you to produce tables of integral'/4X,'or differential particle flux for up to 6 energies'/4X&
               & ,'varying in L-value or normalized magnetic field strength.'/4X,                                                  &
                &'In each of the following windows you will be asked to'/4X,'enter one or more parameters, defining the conditions'&
               & /4X,'for your model tables.'/4X,'In each window the current value(s) is (are) displayed'/4X,                      &
                &'in the right upper corner. You can choose the current'/4X,'value(s) by entering / at the prompt.'/4X,            &
                &'If you enter a wrong character or a value outside the'/4X,'allowed parameter range, the program will ask you for'&
               & /4X,'a new entry.'/4X,'After your tables are displayed you can exit or'/4X,                                       &
                &'calculate tables for new input parameters.'/1X,26('*'),' GOOD LUCK ',26('*'))
         spag_nextblock_1 = 4
         CYCLE SPAG_DispatchLoop_1
      CASE (2)
!--------------------2. window: which parameter to change-----------
         WRITE (monito,99002)
99002    FORMAT (1X/////' WHICH PARAMETER DO YOU WANT TO CHANGE ?'/1X,60('-')/'     0  NO MORE CHANGES, CALCULATE PROFILES'/       &
                &'     1  STORE OR DISPLAY',13X,'4  L-VALUES'/'     2  MODEL',24X,'5  B/B0'/'     3  ENERGIES'/1X,60('-')          &
                &/' ENTER NUMBER')
         maxi = 5
         mini = 0
         spag_nextblock_1 = 3
      CASE (3)
         READ (egnr,*,ERR=20,END=99999) jpara
         IF ( (jpara<=maxi) .AND. (jpara>=mini) ) THEN
            IF ( jpara+1==1 ) THEN
               spag_nextblock_1 = 23
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( jpara+1==3 ) THEN
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( jpara+1==4 ) THEN
               spag_nextblock_1 = 8
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( jpara+1==5 ) THEN
               spag_nextblock_1 = 13
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( jpara+1==6 ) THEN
               spag_nextblock_1 = 18
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 20      WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 3
         CYCLE SPAG_DispatchLoop_1
      CASE (4)
!---------------3. window: display options--------------------------
         WRITE (monito,99003) jagnr
99003    FORMAT (1X/' DO YOU WANT YOUR TABLES',34X,'#',I1,'#'/5X,'DISPLAYED ON YOUR MONITOR:  ENTER  0  AT PROMPT'/5X,             &
                &'STORED IN FILE OUTPUT.AEP:  ENTER  1  AT PROMPT'/5X,'DISPLAYED AND STORED:       ENTER  2  AT PROMPT')
         maxi = 2
         mini = 0
         WRITE (monito,99021)
         spag_nextblock_1 = 5
      CASE (5)
         READ (egnr,*,ERR=40,END=99999) jagnr
         IF ( (jagnr<=maxi) .AND. (jagnr>=mini) ) THEN
            agnr = monito
            IF ( jagnr>0 ) THEN
               OPEN (UNIT=iuout,FILE='OUTPUT.AEP',STATUS='NEW',FORM='FORMATTED')
               IF ( jagnr==1 ) agnr = iuout
            ENDIF
            IF ( .NOT.(notbeg) ) THEN
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 40      WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 5
         CYCLE SPAG_DispatchLoop_1
      CASE (6)
!---------------4. window: model type-------------------------------
         WRITE (monito,99004) mtype
99004    FORMAT (1X/' MTYPE',18X,'* choose model *',18X,'#',I1,'#'/1X,60('-')/' 1   AP8MAX, protons, solar maximum'/' 2   AP8MIN,',&
                &' protons, solar minimum'/' 3   AE4MAX, electrons,',' solar maximum'/' 4   AE4MIN, electrons, solar minimum'/     &
                &' 5   AEI7HI, electrons, high estimate'/' 6   AEI7LO,',' electrons, low estimate'/' 7   AE8MAX, electrons, solar',&
                &' maximum'/' 8   AE8MIN, electrons, solar minimum'//' You must have the corrosponding binary data file in'/       &
                &' your account e.g., AE8MAX.BIN, if you enter MTYPE=7.')
         WRITE (monito,99021)
         maxi = 8
         mini = 1
         spag_nextblock_1 = 7
      CASE (7)
         READ (egnr,*,ERR=60,END=99999) mtype
         IF ( (mtype>=mini) .AND. (mtype<=maxi) ) THEN
            name = mname(mtype)
            WRITE (fname,99005) name
! 99005       FORMAT (A6,'.BIN')
!             OPEN (iuaeap,FILE=fname,STATUS='OLD',ERR=80,FORM='UNFORMATTED')
!             READ (iuaeap) ihead
!             nmap = ihead(8)
!             READ (iuaeap) (map(i),i=1,nmap)
 
! ASCIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! When using the ASCII coefficient files instead of the binary
! coefficient files, one should replace the preceding 5 statements
! with the following 6 statements
!
99005   FORMAT(A6,'.ASC')
      OPEN(IUAEAP,FILE=FNAME,STATUS='OLD',ERR=80,FORM='FORMATTED')
      READ(IUAEAP,1301) IHEAD
      NMAP=IHEAD(8)
      READ(IUAEAP,1301) (MAP(I),I=1,NMAP)
1301   FORMAT(1X,12I6)
 
            CLOSE (iuaeap)
            IF ( mtype<3 ) THEN
               particle = 'PROTONS'
            ELSE
               particle = 'ELECTRONS'
            ENDIF
            IF ( .NOT.(notbeg) ) THEN
               spag_nextblock_1 = 8
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 60      WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 7
         CYCLE SPAG_DispatchLoop_1
 80      WRITE (monito,99006) fname
99006    FORMAT (' Data file ',A10,' is not in your directory,',' try other data file')
         spag_nextblock_1 = 7
         CYCLE SPAG_DispatchLoop_1
      CASE (8)
!-------------------5. window: number of energies----------------------
         WRITE (monito,99007) ne
99007    FORMAT (1X/////' NE',10X,'* Number of Energies *',8X,'<max=6>',7X,'#',I2,'#'/1X,60('-')                                   &
                &/' If you enter 0, you will be asked for the BEGIN, END'/' and STEPSIZE of your energy grid.')
         WRITE (monito,99021)
         maxi = 6
         mini = 0
         spag_nextblock_1 = 9
      CASE (9)
         READ (egnr,*,ERR=100,END=99999) ne
         IF ( (ne>=mini) .AND. (ne<=maxi) ) THEN
            xmax = 7.
            xmin = 0.1
            IF ( mtype<3 ) xmax = 400
            IF ( mtype>2 ) xmin = 0.04
            IF ( ne>0 ) THEN
!--------------------6. window (b): input energies-----------------
               WRITE (monito,99008) ne , (e(j),j=1,6)
99008          FORMAT (1X//' E(1), ... E(',I1,')     * MeV *'//' #',F7.2,5(',',F7.2),'#'/1X,60('-')                                &
                      &/' Enter energies in ascending order.')
               WRITE (monito,99021)
               spag_nextblock_1 = 11
               CYCLE SPAG_DispatchLoop_1
            ELSE
!------------------6. window (a): input energy grid-------------------
               WRITE (monito,99009) ebeg , eend , estep
99009          FORMAT (1X//' BEGIN, END, STEPSIZE    * MeV *     #',F6.2,', ',F6.2,', ',F6.2,'#'/1X,60('-')/3X,                    &
                      &'Begin, end, and stepsize of energy grid in MeV.'/3X,                                                       &
                      &'The number of grid points is restricted to 10 or less !')
               WRITE (monito,99021)
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
 100     WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 9
         CYCLE SPAG_DispatchLoop_1
      CASE (10)
         READ (egnr,*,ERR=120,END=99999) ebeg , eend , estep
         ne = int((eend-ebeg)/estep) + 1
         IF ( ne>6 ) ne = 6
         IF ( (ebeg>=xmin) .AND. (eend<=xmax) ) THEN
            e(1) = ebeg
            DO ie = 2 , ne
               e(ie) = e(ie-1) + estep
            ENDDO
            spag_nextblock_1 = 12
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 120     WRITE (monito,99011) xmin , xmax
         spag_nextblock_1 = 10
         CYCLE SPAG_DispatchLoop_1
      CASE (11)
         READ (egnr,*,ERR=140,END=99999) (e(ie),ie=1,ne)
         IF ( (e(1)>=xmin) .AND. (e(ne)<=xmax) ) THEN
            spag_nextblock_1 = 12
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 140     WRITE (monito,99011) xmin , xmax
         spag_nextblock_1 = 11
         CYCLE SPAG_DispatchLoop_1
      CASE (12)
         DO ei = ne + 1 , 6
            e(ei) = 0.0
         ENDDO
         IF ( notbeg ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 13
      CASE (13)
!-------------------7. window: number of L-values--------------------
         WRITE (monito,99012) nl
99012    FORMAT (1X///' NL',10X,'* Number of L-values *',9X,'<max=10>',4X,'#',I3,'#'/1X,60('-')                                    &
                &/' If you enter 0, you will be asked for the BEGIN, END'/' and STEPSIZE of your L-value grid.')
         WRITE (monito,99021)
         maxi = 10
         mini = 0
         spag_nextblock_1 = 14
      CASE (14)
         READ (egnr,*,ERR=160,END=99999) nl
         IF ( (nl>=mini) .AND. (nl<=maxi) ) THEN
            xmin = 1.
            xmax = 15.
            IF ( nl>0 ) THEN
!----------------8. window (b): L-values-------------------------
               WRITE (monito,99013) nl , (xl(1,j),j=1,10)
99013          FORMAT (1X//' L(1), ... L(',I2,')     * L-values *'//' #',F7.2,5(',',F7.2)/2X,F7.2,3(',',F7.2),'#')
               WRITE (monito,99021)
               spag_nextblock_1 = 16
               CYCLE SPAG_DispatchLoop_1
            ELSE
               WRITE (monito,99014) vbeg , vend , vstep
!----------------8. window (a): L grid---------------------------------
99014          FORMAT (1X//' BEGIN, END, STEPSIZE     * L *      #',F6.2,', ',F6.2,', ',F6.2,'#'/1X,60('-')/3X,                    &
                      &'Begin, end, and stepsize of L-value grid.'/3X,'The number of grid points is restricted to 10 or less !')
               WRITE (monito,99021)
               spag_nextblock_1 = 15
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
 160     WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 14
         CYCLE SPAG_DispatchLoop_1
      CASE (15)
         READ (egnr,*,ERR=180,END=99999) vbeg , vend , vstep
         IF ( (vbeg>=xmin) .AND. (vend<=xmax) ) THEN
            nl = int((vend-vbeg)/vstep) + 1
            IF ( nl>10 ) nl = 10
            xl(1,1) = vbeg
            DO ie = 2 , nl
               xl(1,ie) = xl(1,ie-1) + vstep
            ENDDO
            spag_nextblock_1 = 17
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 180     WRITE (monito,99011) xmin , xmax
         spag_nextblock_1 = 15
         CYCLE SPAG_DispatchLoop_1
      CASE (16)
         READ (egnr,*,ERR=200,END=99999) (xl(1,ie),ie=1,nl)
         IF ( (xl(1,1)>=xmin) .AND. (xl(1,nl)<=xmax) ) THEN
            spag_nextblock_1 = 17
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 200     WRITE (monito,99011) xmin , xmax
         spag_nextblock_1 = 16
         CYCLE SPAG_DispatchLoop_1
      CASE (17)
         DO ie = nl + 1 , 10
            xl(1,ie) = 0.0
         ENDDO
         IF ( notbeg ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 18
      CASE (18)
!-----------------9. window: number of B/B0 values------------------
         WRITE (monito,99015) nb
99015    FORMAT (1X/////' NB',9X,'* Number of B/B0 values *',8X,'<max=10>',3X,'#',I3,'#'/1X,60('-')                                &
                &/' If you enter 0, you will be asked for the BEGIN, END'/' and STEPSIZE of your B/B0 value grid.')
         WRITE (monito,99021)
         maxi = 10
         mini = 0
         spag_nextblock_1 = 19
      CASE (19)
         READ (egnr,*,ERR=220,END=99999) nb
         IF ( (nb>=mini) .AND. (nb<=maxi) ) THEN
            xmax = 9999.9
            xmin = 1.
            IF ( nb>0 ) THEN
!-----------10. window (b): B/B0 values------------------------
               WRITE (monito,99016) nb , (xl(2,j),j=1,10)
99016          FORMAT (1X//' B/B0(1), ... B/B0(',I2,')'//' #',F7.2,5(',',F7.2)/2X,F7.2,3(',',F7.2),'#')
               WRITE (monito,99021)
               spag_nextblock_1 = 21
               CYCLE SPAG_DispatchLoop_1
            ELSE
!------------10. window (a): B/B0 grid--------------------------------
               WRITE (monito,99017) bbeg , bend , bstep
99017          FORMAT (1X//' BEGIN, END, STEPSIZE    * B/B0 *    #',F7.2,',',F7.2,',',F7.2,'#'/1X,60('-')/3X,                      &
                      &'Begin, end, and stepsize of B/B0 grid.'/3X,'The number of grid points is restricted to 10 or less !')
               WRITE (monito,99021)
               spag_nextblock_1 = 20
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
 220     WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 19
         CYCLE SPAG_DispatchLoop_1
      CASE (20)
         READ (egnr,*,ERR=240,END=99999) bbeg , bend , bstep
         nb = int((bend-bbeg)/bstep) + 1
         IF ( nb>10 ) nb = 10
         IF ( (bbeg>=xmin) .AND. (bend<=xmax) ) THEN
            xl(2,1) = bbeg
            DO ie = 2 , nb
               xl(2,ie) = xl(2,ie-1) + bstep
            ENDDO
            spag_nextblock_1 = 22
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 240     WRITE (monito,99011) xmin , xmax
         spag_nextblock_1 = 20
         CYCLE SPAG_DispatchLoop_1
      CASE (21)
         READ (egnr,*,ERR=260,END=99999) (xl(2,ie),ie=1,nb)
         IF ( (xl(2,1)>=xmin) .AND. (xl(2,nb)<=xmax) ) THEN
            spag_nextblock_1 = 22
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 260     WRITE (monito,99011) xmin , xmax
         spag_nextblock_1 = 21
         CYCLE SPAG_DispatchLoop_1
      CASE (22)
         DO ie = nb + 1 , 10
            xl(2,ie) = 0.0
         ENDDO
         IF ( notbeg ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         notbeg = .TRUE.
         spag_nextblock_1 = 23
      CASE (23)
!----------------THE L-VALUE LOOP-------------------------------
         DO l = 1 , nl
            fl = xl(1,l)
!----------------THE B LOOP-------------------------------------
            DO i = 1 , nb
               bb0 = xl(2,i)
               CALL trara1(ihead,map,fl,bb0,e,flux,ne)
!----------------THE ENERGY LOOP--------------------------------
               DO k = 1 , ne
                  af(k,l,i) = 0.0
                  IF ( flux(k)>0.0 ) af(k,l,i) = 10.**flux(k)
                  IF ( k/=1 ) THEN
                     df(k-1,l,i) = abs(af(k,l,i)-af(k-1,l,i))/(e(k)-e(k-1))
                     IF ( af(k,l,i)<=0.0 ) df(k-1,l,i) = 0.0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!-----------------TESTS VALIDITY OF DIFFERENTIAL FLUX------------
         DO k = 1 , ne - 1
            ei = e(k+1)
            itest = 0
            ediff = ei - e(k)
!-----------------IS ENERGY INTERVALL LARGE ENOUGH ?-------------
!            IF ((EI.GT.0.25).AND.(EDIFF.LT.0.1999)) ITEST=1
!            IF ((EI.LE.0.10).AND.(EDIFF.LT.0.0499)) ITEST=1
!            IF(((EI.LE.0.25).AND.(EDIFF.LT.0.0999)
!     &            .AND.(EI.GT.0.1))) ITEST=1
            iei = 1
            IF ( ei>0.10 ) iei = 2
            IF ( ei>0.25 ) iei = 3
            IF ( ediff<eda(iei) ) itest = 1
            DO l = 1 , nl
               itt = 0
               IF ( xl(1,l)<1.2 ) itt = 1
               DO i = 1 , nb
                  ite(k,l,i) = ' '
                  IF ( itest+itt/=0 ) ite(k,l,i) = '?'
               ENDDO
            ENDDO
         ENDDO
!------------------TABLE OUTPUT------------------------------------
         WRITE (monito,99018)
99018    FORMAT (1X/)
         spag_nextblock_1 = 24
      CASE (24)
         WRITE (monito,99019) jtab
99019    FORMAT (' TABLE OPTIONS:',6X,'integral',12X,'differential flux',5X,'#',I1,'#'/21X,40('-')/2X,                             &
                &'0  EXIT            1  E * L TABLE      3  E * L TABLE'/2X,                                                       &
                &'5  NEW PARAMETER   2  E * B/B0 TABLE   4  E * B/B0 TABLE')
         WRITE (monito,99021)
         mini = 0
         maxi = 5
         spag_nextblock_1 = 25
      CASE (25)
         READ (egnr,*,ERR=280,END=99999) jtab
         IF ( (jtab<=maxi) .AND. (jtab>=mini) ) THEN
            IF ( jtab==0 ) STOP
            IF ( jtab==5 ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( jtab>2 ) THEN
               untex = '*MeV]'
               idtex = 'differential'
            ELSE
               untex = ']    '
               idtex = 'integral'
            ENDIF
            IF ( (jtab==1) .OR. (jtab==3) ) THEN
               ibl = 2
               bltex = 'B/B0'
               lbtex = ' L  '
               n = nl
               nn = nb
            ELSE
               ibl = 1
               bltex = ' L  '
               lbtex = 'B/B0'
               n = nb
               nn = nl
            ENDIF
            ibltab = 1
            IF ( ibl==1 ) ibltab = 2
            IF ( nn==1 ) THEN
               spag_nextblock_1 = 27
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            WRITE (monito,99022) bltex , inde , bltex , (xl(ibl,j),j=1,10)
99022       FORMAT (1X/' FOR WHICH ',A4,' ?         * enter index number,',' not value*',9X,'#',I2,'#'/1X,70('-')/' Nr:',5X,       &
                   &'1     2     3     4     5     6     7   ','  8     9    10'/1X,A4,':',10F6.1)
            WRITE (monito,99020)
99020       FORMAT (1X,70('-')/,' enter /, to continue with current value(s);',' Ctrl Z, to exit')
            mini = 1
            maxi = 10
            spag_nextblock_1 = 26
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 280     WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 25
         CYCLE SPAG_DispatchLoop_1
      CASE (26)
         READ (egnr,*,ERR=300,END=99999) inde
         IF ( (inde<=maxi) .AND. (inde>=mini) ) THEN
            spag_nextblock_1 = 27
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 300     WRITE (monito,99010) mini , maxi
         spag_nextblock_1 = 26
         CYCLE SPAG_DispatchLoop_1
      CASE (27)
         WRITE (agnr,99023) idtex , particle , untex , lbtex , (e(j),j=1,6)
         IF ( jagnr==2 ) WRITE (iuout,99023) idtex , particle , untex , lbtex , (e(j),j=1,6)
         blv = xl(ibl,inde)
         DO i = 1 , n
            il = inde
            ib = inde
            IF ( ibl==1 ) THEN
               ib = i
               xx = xl(2,i)
            ELSE
               il = i
               xx = xl(1,i)
            ENDIF
            IF ( jtab<3 ) THEN
               WRITE (agnr,99024) xx , (af(j,il,ib),j=1,6)
               IF ( jagnr==2 ) WRITE (iuout,99024) xx , (af(j,il,ib),j=1,6)
            ELSE
               WRITE (agnr,99025) xx , df(1,il,ib) , ite(1,il,ib) , df(2,il,ib) , ite(2,il,ib) , df(3,il,ib) , ite(3,il,ib) ,      &
                                & df(4,il,ib) , ite(4,il,ib) , df(5,il,ib) , ite(5,il,ib)
               IF ( jagnr==2 ) WRITE (iuout,99025) xx , df(1,il,ib) , ite(1,il,ib) , df(2,il,ib) , ite(2,il,ib) , df(3,il,ib) ,    &
                                    & ite(3,il,ib) , df(4,il,ib) , ite(4,il,ib) , df(5,il,ib) , ite(5,il,ib)
            ENDIF
         ENDDO
         WRITE (agnr,99026) name , bltex , blv
         IF ( jagnr==2 ) WRITE (iuout,99026) name , bltex , blv
         spag_nextblock_1 = 24
         CYCLE SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99010 FORMAT (' Your input is ouside the value range:',I4,' to',I4,', try again')
99011 FORMAT (' Your input is ouside the value range:',F8.2,' to ',F8.2,', try again')
99021 FORMAT (1X,60('-')/,' enter /, to continue with current value(s);',' Ctrl Z, to exit')
99023 FORMAT (1X/25X,A12,' flux [',A9,'/cm*cm*sec',A5/1X,A4,'\\ E/MeV',F7.2,5F10.2/1X,'------\\',62('-'))
99024 FORMAT (1X,F7.2,' I  ',1PE9.3,5(1X,1PE9.3))
99025 FORMAT (1X,F7.2,' I  ',5X,5(1PE9.2,A1))
99026 FORMAT (1X,20('-'),' flux values below 10 are not reliable ',11('-')/1X,'Model: ',A6,10X,A4,'=',F8.2/1X,70('*'))
99999 END PROGRAM spag_program_1
