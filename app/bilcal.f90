!*==spag_program_1.f90 processed by SPAG 8.01MH 22:19  8 Feb 2024
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
PROGRAM spag_program_1

   use radbelt_kinds_module
   use shellig_module

   IMPLICIT NONE

   type(shellig_type) :: igrf
!*** Start of declarations inserted by SPAG
   REAL(wp) alog2 , Aquad , bab1 , babs , bbx , bdel , bdown , &
            beast , beq , bequ , bnorth , Bquad , bvar , dec , dimo , dip , Era ,   &
            evar , height , rr0
   REAL(wp) svar , Umr , vamax , vamin , varb , vare , xcor , xl , xmax , xmin , xvar , year
   INTEGER iall , ibbb , icode , imax , imin , istart , iswit , ivar , ivarnr , jagnr , lanz , lfd , monito
!*** End of declarations inserted by SPAG
! BILCAL, VERSION 3.0, AUGUST 1995
!
!mm/dd/yy
! 1/25/92-DKB-Modified for use with the IGRF-91 coefficients, which
!	      were provided by R. Langel, GSFC.
! 2/ 5/92-DKB-Reduce variable-name: INITI(ALI)ZE
! 3/25/96-DKB-Modified for use with the IGRF-95 coefficients, which
!	      were provided by R. Langel, GSFC.
! 6/ 6/00-DKB-Modified for use with IGRF-2000 coefficients.
!11/14/01-DKB-Add IMIN=0 above 4927 READ(...)  [Rui Pereira]
!04/25/05-DKB-IBBB instead of IBB in data statem.  [Alexey Petrov]
!
!*****************************************************************
!**************** IGRF MAGNETIC FIELD MODEL  *********************
!**************** SHELLG L-VALUE CALCULATION *********************
!*****************************************************************
!*****************************************************************
!*** THIS PROGRAM PRODUCES PROFILES OF:                        ***
!***      GEOMAGNETIC FIELD STRENGTH (GAUSS) 		       ***
!***      L-VALUE 					       ***
!*****************************************************************
!*** FOR SPECIFIED:                                            ***
!***      YEAR (DECIMAL YEAR, E.G., 1995.5 FOR MID 1995)       ***
!***      GEODATIC LATITUDE AND LONGITUDE (DEGREE)             ***
!***      ALTITUDE (KM)                                        ***
!*****************************************************************
!*****************************************************************
!*     --------------------ADDRESS--------------------------     *
!*     I  DR. DIETER BILITZA  (301)513-1664      	   I     *
!*     I  GSFC, NSSDC, CODE 933, GREENBELT, MD 20771, USA  I     *
!*     I  SPAN:     NSSDCA::BILITZA, NSSDC::BILITZA        I     *
!*     I  BITNET:   BILITZA%NSSDCA.SPAN@DFTNIC.BITNET      I     *
!*     -----------------------------------------------------     *
!*****************************************************************
!*****************************************************************
!*****************************************************************
   INTEGER egnr , agnr , ognr
   REAL(wp) lati , longi
   CHARACTER*4 itext(4)
   CHARACTER*7 itb
   LOGICAL notbeg , val
   DIMENSION xvar(4) , vare(4) , varb(4)
   COMMON /gener / Umr , Era , Aquad , Bquad
   INTEGER :: spag_nextblock_1
!
   DATA itext/'LATI' , 'LONG' , 'H/km' , 'YEAR'/
   DATA lati , longi , height , year , ivar , bvar , evar , svar , &
        ibbb , jagnr/45.1 , 293.1 , 100 , 1985.5 , 3 , 100 , 1000 , &
        100 , 0 , 2/
!### year limit modified
   DATA varb/ - 90.0_wp , -360.0_wp , 0.00000_wp , 1940.0_wp/
   DATA vare/ + 90.0_wp , +360.0_wp , 30000.0_wp , 2010.0_wp/
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
         !CALL initize
         alog2 = log(2.0_wp)
         istart = 1
!
! FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
! EGNR=INPUT, MONITO=MONITOR, KONSOL=MESSAGES......................
! AGNR=DISPLAY, OGNR=FILE OUTPUT...................................
!
         egnr = 5
         monito = 6
         agnr = 6
         ognr = 16
         WRITE (monito,99001)
99001    FORMAT (1X/////4X,54('*')/4X,'****** IGRF GEOMAGNETIC FIELD MODEL 1945 - 2005 ******'/4X,                                 &
                &'***********  SHELLG L-VALUE CALCULATION  *************'/1X,60('*')                                               &
                &/'  This program allows you to produce B and L ','profiles in '/'  latitude, longitude, year or altitude.'/       &
                &'  In each of the following windows you will be ','asked to enter'/'  one or more values, defining the conditions'&
               & ,' for your tables.'/'  In each window the current value(s) is',                                                  &
                &' (are) shown in the right'/'  upper corner (#...#). You can ',                                                   &
                &'choose the current values by'/'  entering / at the prompt.'/                                                     &
                &'  If you enter a wrong character or a value outside the ',                                                       &
                &'allowed'/'  parameter range, the program will ask you for a',                                                    &
                &' new entry.'/'  After your tables are displayed, you can ',                                                      &
                &'change any parameter'/'  you wish to change and create a ',                                                      &
                &'new profile.'/'  You can leave the program at any point ','by entering Ctrl Z.'/1X,25('*'),' GOOD LUCK ',25('*'))
         notbeg = .FALSE.
         spag_nextblock_1 = 4
         CYCLE SPAG_DispatchLoop_1
      CASE (2)
!---------------START ENTERING PARAMETERS----------------------------
         istart = istart + 1
!---------------WINDOW 1: WHICH PARAMETER CHANGE ?-------------------
         WRITE (monito,99002) lati , longi , itext(ivar) , height , year , bvar , evar , svar
99002    FORMAT (1X//' **** WHICH PARAMETER DO YOU WANT TO CHANGE?'/1X,60('-')/' 0  NO FURTHER CHANGES, CALCULATE PROFILE'/        &
                &' 1  LATITUDE  #',F6.1,'#',7X,'5  DISPLAY OR STORE'/' 2  LONGITUDE #',F6.1,'#',7X,'6  SELECTION OF VARIABLE  #',  &
               & A4,'#'/' 3  ALTITUDE  #',F8.1,'#',5X,'7  VARIABLE RANGE'/' 4  YEAR      #',F6.1,'#',11X,'#',F8.1,',',F8.1,',',    &
               & F8.1,'#'/29X,'8  B OR B/B0'/1X,60('-')/' ENTER NUMBER')
         imin = 0
         imax = 8
         spag_nextblock_1 = 3
      CASE (3)
         READ (egnr,*,ERR=20,END=99999) iswit
         IF ( (iswit>=imin) .AND. (iswit<=imax) ) THEN
            IF ( iswit+1==1 ) THEN
               spag_nextblock_1 = 23
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( iswit+1==2 ) THEN
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( iswit+1==3 ) THEN
               spag_nextblock_1 = 13
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( iswit+1==4 ) THEN
               spag_nextblock_1 = 16
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( iswit+1==5 ) THEN
               spag_nextblock_1 = 19
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( iswit+1==7 ) THEN
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( iswit+1==8 ) THEN
               spag_nextblock_1 = 8
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            IF ( iswit+1==9 ) THEN
               spag_nextblock_1 = 21
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 20      WRITE (monito,99003) imin , imax
         spag_nextblock_1 = 3
         CYCLE SPAG_DispatchLoop_1
      CASE (4)
!--------------WINDOW 2: DISPLAY OPTIONS--------------------------
         WRITE (monito,99004) jagnr
99004    FORMAT (/' DO YOU WANT YOUR PROFILES',32X,'#',I1,'#'/5X,'DISPLAYED ON YOUR MONITOR: ENTER  0  AT PROMPT'/5X,              &
                &'STORED IN FILE OUTPUT.IGR: ENTER  1  AT PROMPT'/5X,'DISPLAYED AND STORED:      ENTER  2  AT PROMPT')
         WRITE (monito,99006)
         imax = 2
         imin = 0
         spag_nextblock_1 = 5
      CASE (5)
         READ (egnr,*,ERR=40,END=99999) jagnr
         IF ( (jagnr>=imin) .AND. (jagnr<=imax) ) THEN
            ivarnr = 0
            IF ( jagnr>0 ) OPEN (UNIT=ognr,FILE='OUTPUT.IGR',STATUS='NEW',FORM='FORMATTED')
            IF ( jagnr==1 ) agnr = ognr
            IF ( .NOT.(notbeg) ) THEN
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 40      WRITE (monito,99003) imin , imax
         spag_nextblock_1 = 5
         CYCLE SPAG_DispatchLoop_1
      CASE (6)
!---------------WINDOW 3: SELECT VARIABLE------------------------
         WRITE (monito,99005) ivar
99005    FORMAT (1X//' SELECT YOUR VARIABLE:',31X,'#LAST:',I1,'#'//' 1  LATITUDE       3  ALTITUDE'/' 2  LONGITUDE      4  YEAR')
         WRITE (monito,99006)
         imin = 1
         imax = 4
         spag_nextblock_1 = 7
      CASE (7)
         READ (egnr,*,ERR=60,END=99999) ivar
         IF ( (ivar>=imin) .AND. (ivar<=imax) ) THEN
            spag_nextblock_1 = 8
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 60      WRITE (monito,99003) imin , imax
         spag_nextblock_1 = 7
         CYCLE SPAG_DispatchLoop_1
      CASE (8)
!--------------WINDOW 4: SELECT VARIABLE RANGE---------------------
         WRITE (monito,99007) bvar , evar , svar
99007    FORMAT (1X//' CHOOSE YOUR VARIABLE RANGE:',5X,' BEGIN, END, ','STEPWIDTH ?'/32X,'#',F8.1,',',F8.1,',',F8.1,'#')
         WRITE (monito,99006)
         vamin = varb(ivar)
         vamax = vare(ivar)
         spag_nextblock_1 = 9
      CASE (9)
         READ (egnr,*,ERR=80,END=99999) bvar , evar , svar
         IF ( (bvar>=vamin) .AND. (evar<=vamax) ) THEN
            lanz = int((evar-bvar)/svar) + 1
            IF ( notbeg ) THEN
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            ivarnr = ivarnr + 1
            IF ( ivarnr/=ivar ) THEN
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 12
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 80      WRITE (monito,99008) vamin , vamax
         spag_nextblock_1 = 9
         CYCLE SPAG_DispatchLoop_1
      CASE (10)
!--------------WINDOW 5: LATITUDE-----------------------------------
         WRITE (monito,99009) lati
99009    FORMAT (1X//1X,'GEOD LATITUDE ?   !NORTH!    [DEGREE,DECIMAL]',8X,'#',F5.1,'#')
         WRITE (monito,99006)
         xmax = vare(1)
         xmin = varb(1)
         spag_nextblock_1 = 11
      CASE (11)
         READ (egnr,*,ERR=100,END=99999) lati
         IF ( (lati>=xmin) .AND. (lati<=xmax) ) THEN
            IF ( .NOT.(notbeg) ) THEN
               spag_nextblock_1 = 12
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 100     WRITE (monito,99008) xmin , xmax
         spag_nextblock_1 = 11
         CYCLE SPAG_DispatchLoop_1
      CASE (12)
         ivarnr = ivarnr + 1
         IF ( ivarnr==ivar ) THEN
            spag_nextblock_1 = 15
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 13
      CASE (13)
!---------------WINDOW 6: LONGITUDE---------------------------------
         WRITE (monito,99010) longi
99010    FORMAT (1X//1X,'GEOD LONGITUDE ?   !EAST!    [DEGREE,DECIMAL]',7X,'#',F6.1,'#')
         WRITE (monito,99006)
         xmax = vare(2)
         xmin = varb(2)
         spag_nextblock_1 = 14
      CASE (14)
         READ (egnr,*,ERR=120,END=99999) longi
         IF ( (longi>=xmin) .AND. (longi<=xmax) ) THEN
            IF ( .NOT.(notbeg) ) THEN
               spag_nextblock_1 = 15
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 120     WRITE (monito,99008) xmin , xmax
         spag_nextblock_1 = 14
         CYCLE SPAG_DispatchLoop_1
      CASE (15)
         ivarnr = ivarnr + 1
         IF ( ivarnr==ivar ) THEN
            spag_nextblock_1 = 18
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 16
      CASE (16)
!---------------WINDOW 7: ALTITUDE---------------------------------
         WRITE (monito,99011) height
99011    FORMAT (1X//1X,'ALTITUDE ?    [KM]',33X,'#',F7.1,'#')
         WRITE (monito,99006)
         xmax = vare(3)
         xmin = varb(3)
         spag_nextblock_1 = 17
      CASE (17)
         READ (egnr,*,ERR=140,END=99999) height
         IF ( (height>=xmin) .AND. (height<=xmax) ) THEN
            IF ( .NOT.(notbeg) ) THEN
               spag_nextblock_1 = 18
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 140     WRITE (monito,99008) xmin , xmax
         spag_nextblock_1 = 17
         CYCLE SPAG_DispatchLoop_1
      CASE (18)
         ivarnr = ivarnr + 1
         IF ( ivarnr==ivar ) THEN
            spag_nextblock_1 = 21
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 19
      CASE (19)
!----------------WINDOW 8: YEAR------------------------------------
         WRITE (monito,99012) year
99012    FORMAT (1X//' YEAR(EPOCH) ?',9X,'*decimal*',9X,'#',F6.1,'#')
         WRITE (monito,99006)
         xmax = vare(4)
         xmin = varb(4)
         spag_nextblock_1 = 20
      CASE (20)
         READ (egnr,*,ERR=160,END=99999) year
         IF ( (year>=xmin) .AND. (year<=xmax) ) THEN
            IF ( .NOT.(notbeg) ) THEN
               spag_nextblock_1 = 21
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 160     WRITE (monito,99008) xmin , xmax
         spag_nextblock_1 = 20
         CYCLE SPAG_DispatchLoop_1
      CASE (21)
!----------------WINDOW 9: ABSOLUTE OR NORMALIZED B--------------
         WRITE (monito,99013) ibbb
99013    FORMAT (1X//' OUTPUT OPTION: B OR B/B0 ?',19X,'#',I1,'#'//4X,'if you enter 0, the absolute magnetic field strength'/4X,   &
                &'will be listed, otherwise the field strength normalized'/4X,                                                     &
                &'to the field strength at the magnetic equator is listed')
         WRITE (monito,99006)
         spag_nextblock_1 = 22
      CASE (22)
         READ (egnr,*,ERR=180,END=99999) ibbb
         IF ( ibbb/=0 ) THEN
            itb = '  B/B0 '
         ELSE
            itb = 'B/Gauss'
         ENDIF
         IF ( .NOT.(notbeg) ) THEN
            spag_nextblock_1 = 23
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 2
         CYCLE SPAG_DispatchLoop_1
 180     WRITE (monito,99014)
99014    FORMAT (' Your input should be a integer value'/' try again')
         spag_nextblock_1 = 22
         CYCLE SPAG_DispatchLoop_1
      CASE (23)
!----------------CALCULATE PROFILES-----------------------------------
         WRITE (agnr,99015) itext(ivar) , itb
         IF ( jagnr==2 ) WRITE (ognr,99015) itext(ivar) , itb
         xvar(1) = lati
         xvar(2) = longi
         xvar(3) = height
         xvar(4) = year
         lfd = 0
         xvar(ivar) = bvar - svar
         SPAG_Loop_1_1: DO
            xvar(ivar) = xvar(ivar) + svar
            lfd = lfd + 1
            lati = xvar(1)
            longi = xvar(2)
            height = xvar(3)
            year = xvar(4)
            IF ( (ivar>=4) .OR. (lfd<=1) ) CALL igrf%feldcof(year,dimo)
            CALL igrf%feldg(lati,longi,height,bnorth,beast,bdown,babs)
            CALL igrf%shellg(lati,longi,height,dimo,xl,icode,bab1)
            IF ( iabs(icode)>9 ) THEN
               WRITE (monito,99016) icode
99016          FORMAT (' ICODE=',I10,' is set to 2')
               icode = 2
            ENDIF
            IF ( ibbb/=0 ) THEN
               bequ = dimo/(xl*xl*xl)
               IF ( icode==1 ) THEN
                  bdel = 1.0E-3_wp
                  CALL igrf%findb0(0.05_wp,bdel,val,beq,rr0)
                  IF ( val ) bequ = beq
               ENDIF
            ENDIF
            dip = asin(bdown/babs)/Umr
            dec = asin(beast/sqrt(beast*beast+bnorth*bnorth))/Umr
            xcor = xvar(ivar)
            IF ( ibbb==0 ) THEN
               WRITE (agnr,99017) xcor , dimo , babs , bnorth , beast , bdown , dip , dec , xl , icode
               IF ( jagnr==2 ) WRITE (ognr,99017) xcor , dimo , babs , bnorth , beast , bdown , dip , dec , xl , icode
            ELSE
               bbx = babs/bequ
               IF ( bbx>9999.999 ) bbx = 9999.999
               WRITE (agnr,99018) xcor , dimo , bbx , bnorth , beast , bdown , dip , dec , xl , icode
               IF ( jagnr==2 ) WRITE (ognr,99018) xcor , dimo , bbx , bnorth , beast , bdown , dip , dec , xl , icode
            ENDIF
            IF ( xcor>=evar ) THEN
               WRITE (agnr,99019) lati , longi , height , year
               IF ( jagnr==2 ) WRITE (ognr,99019) lati , longi , height , year
               IF ( height>5000.0 ) THEN
                  WRITE (agnr,99020)
                  IF ( jagnr==2 ) WRITE (ognr,99020)
               ENDIF
! ### year limits corrected
               IF ( (year<1945.0_wp) .OR. (year>2005.0_wp) ) THEN
                  WRITE (agnr,99021)
                  IF ( jagnr==2 ) WRITE (ognr,99021)
               ENDIF
!-----------------LAST WINDOW: CONTINUE ?-----------------------
               WRITE (monito,99022)
99022          FORMAT (1X/' **** DO YOU WANT TO CONTINUE?'/1X,60('-')/' "0"  QUIT AND EXIT        "1"  NEW PARAMETERS'/1X,60('-'))
               imin = 0
               imax = 1
               EXIT SPAG_Loop_1_1
            ENDIF
         ENDDO SPAG_Loop_1_1
         spag_nextblock_1 = 24
      CASE (24)
         READ (egnr,*,ERR=200,END=99999) iall
         IF ( (iall>=imin) .AND. (iall<=imax) ) THEN
            notbeg = .TRUE.
            IF ( iall/=1 ) STOP
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
 200     WRITE (monito,99003) imin , imax
         spag_nextblock_1 = 24
         CYCLE SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99003 FORMAT (' Your input is outside the value range:',I2,' to',I2/' try again')
99006 FORMAT (1X,60('-')/' Enter / to use previous value(s) ','(see # .. #); Ctrl Z to exit')
99008 FORMAT (' Your input is outside the value range:',F8.1,' to',F8.1/' try again')
99015 FORMAT (1X////////////5X,A4,'   DIMO  ',A7,' B-NORTH  B-EAST  B-DOWN ','   DIP    DEC  L-VALUE C')
99017 FORMAT (1X,F8.2,F8.4,4(1X,F7.5),2F7.1,F8.3,I3)
99018 FORMAT (1X,F8.2,F8.4,F8.3,3(1X,F7.5),2F7.1,F8.3,I3)
! ### edition date corrected
99019 FORMAT (1X,'------- International Geomagnetic Reference Field',' --- Edition 2000 -------'/' LATI=',F7.1,'  LONGI=',F6.1,    &
             &'  I   DIMO is Dipol   I   C=1  L and B0 correct'/'  ALT=',F7.1,'   YEAR=',F6.1,'  I  Moment in Gauss',              &
             &'  I    =2  wrong,  =3  approx.'/1X,74('-'))
99020 FORMAT (' !! REMINDER: this field model does not',' include external sources !!')
! ### timeperiod corrected
99021 FORMAT (' !! REMINDER: Recommended time period is 1945',' to 2005 !!')
99999 END PROGRAM spag_program_1
