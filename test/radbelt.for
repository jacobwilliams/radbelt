C RADBELT.FOR	SEPTEMBER 88
C
C***********************************************************************
C*                                                                     *
C*                    TRAPPED RADIATION MODELS                         *
C*                         Program RADBELT                             *
C*                                                                     *
C***  Dieter Bilitza  ******************************  25 March 1988  ***
C***********************************************************************
C***************  Goddard Space Flight Center, code 633  ***************
C**  National Space Science Data Center, Greenbelt, MD 20771, U.S.A.  **
C***********************************************************************
C**********************  NSSDC-ID  PT-14A  *****************************
C***********************************************************************
C***   This program  gives an example of how to use NSSDC's Trapped  ***
C***   Radiation Models. It determines omnidirectional integral and  ***
C***   differential electron or proton fluxes for given energies,    ***
C***   L-values, magnetic field strengths and map-type.              ***
C***   The program will ask you for:                                 ***
C***      NE 	 number of energies you want displayed               ***
C***      E(1),... E(NE)   energies   or   EBEGIN,EEND,ESTEP  begin, *** 
C***             end and stepsize of energy range                    ***
C***      NL     number of L-values                                  ***
C***      L(1),... L(NL)   L-values   or   LBEGIN,LEND,LSTEP         ***
C***      NB     number of B/B0-values                               ***
C***      B/B0(1),... B/B0(NL)   B/B0-values   or  BBEGIN,BEND,BSTEP ***
C***      MTYPE  map type:  1 AP8MAX   2 AP8MIN  3 AE4MAX  4 AE4MIN  ***
C***                        5 AEI7HI   6 AEI7LO  7 AE8MAX  8 AE8MIN  ***
C***      JTAB 	 output options: integral or differential fluxes     ***
C***             versus L or B/B0                                    ***
C***   The program interpolates the NSSDC model map in B/B0, L       ***
C***   (function TRARA2), and energy (subroutine TRARA1).            ***
C***   The program opens the map as binary data file (e.g.           ***
C***   AE8MAX.BIN) on unit 15. Make sure that the map you want to    ***
C***   use is available under the correct name in your account       ***
C***   ( see MTYPE for available models and names ).                 ***
C***********************************************************************
C***********************************************************************
C**************************  NOMENCLATURE  *****************************
C********  omnidirectional flux is flux per unit time and      *********
C********      unit sphere surface.                            *********
C********  integral flux for energy E is flux per unit time    *********
C********      and unit surface of particles with energies     *********
C********      greater than E.                                 *********
C********  differential flux for energy E and energy range DE  *********
C********      is average flux per unit time, unit surface and *********
C********      unit energy of particles with energies between  *********
C********      E-DE/2 and E+DE/2.                              *********
C***********************************************************************
C******************************  UNITS  ********************************
C********                 energies: MeV                        *********
C********          integral fluxes: particles/(cm*cm*sec)      ********* 
C********      differential fluxes: particles/(cm*cm*sec*MeV)  *********
C********  magnetic field strength: Gauss                      *********
C***********************************************************************
C***********************************************************************
C*************  DESCRIPTION OF MODEL DATA FILE FORMAT  *****************
C***********************************************************************
C***  THE FILE CONSISTS OF A HEADER ARRAY (IHEAD(8)) AND A MODEL MAP ***
C***  ARRAY (MAP(...)). ALL ELEMENTS ARE INTEGER.                    ***
C***                                                                 ***
C***  IHEAD(1)   MODEL MAP TYPE (SEE ABOVE)                          ***
C***       (2)   INCREMENTS PER DECADE OF LOGARITHMIC FLUX           ***
C***       (3)   EPOCH OF MODEL                                      ***
C***       (4)   SCALE FACTOR FOR ENERGY; E/MEV=E(MAP)/IHEAD(4)      ***
C***		   			  =6400 (AE-8),  =100 (AP-8) ***
C***       (5)   SCALE FACTOR FOR L-VALUE =2100 (AE-8), =2048 (AP-8) ***
C***       (6)   SCALE FACTOR FOR B/B0    =1024 (AE-8), =2048 (AP-8) ***
C***       (7)   SCALE FACTOR FOR LOGARITHM OF FLUXES =1024 (AE,AP-8)***
C***       (8)   NUMBER OF ELEMENTS IN MAP =13548 (AE8MAX),	     ***
C***		  =13168 (AE8MIN), =6509 (AP8MAX), =6688 (AP8MIN)    ***
C***                                                                 ***
C***  LAYOUT OF MAP:                                                 ***
C***      MAP CONSISTS OF SEVERAL VARIABLE-LENGTH SUB-MAPS, EACH     ***
C***      FOR A DIFFERENT ENERGY. EACH SUB-MAP CONSISTS OF SEVERAL   ***
C***      VARIABLE-LENGTH SUB-SUB-MAPS EACH FOR A DIFFERENT L-VALUE. ***
C***      EACH SUB-SUB-MAP CONTAINS THE CURVE LOG(F) [DECADIC        ***
C***      LOGARITHM OF OMNIDIRECTIONAL INTEGRAL PARTICLE FLUX]       ***
C***      VERSUS B/B0 [MAGNETIC FIELD STRENGTH NORMALIZED TO THE     ***
C***      EQUATORIAL VALUE]. THE CURVE IS PARAMETERIZED BY USING     ***
C***      EQUAL INCREMENTS IN LOG(F); THE NUMBER OF INCREMENTS       ***
C***      PER DECADE IS LISTED IN THE HEADER ARRAY [IHEAD(2)]:       ***
C***                                                                 ***
C***         I     B(I)/B(0)   (B(I)-B(I-1))/B(0)   LOG(F(I))        ***
C***       ----------------------------------------------------      ***
C***         0        1                 -              Y             ***
C***         1     B(1)/B(0)   (B(1)-B(0))/B(0)      Y-1/IHEAD(2)    ***
C***         2     B(2)/B(0)   (B(2)-B(1))/B(0)      Y-2/IHEAD(2)    ***
C***         .       ....            .....            ....           ***
C***                                                                 ***
C***      THE SUB-SUB-MAP CONTAINS THE EQUATORIAL FLUX LOGARITHM Y   ***
C***      AND THE B/B0-INCREMENTS (THIRD COLUMN) MULTIPLIED BY       ***
C***      THEIR CORRESPONDING SCALE VALUES ( IHEAD(7) AND (8) ).     ***
C***                                                                 ***
C***  MAP(1)  NUMBER OF ELEMENTS IN SUB-MAP                          ***
C***  MAP(2)  ENERGY FOR THIS SUB-MAP; MAP(2)=E/MEV*IHEAD(4)         ***
C***    MAP(3)  NUMBER OF ELEMENTS IN SUB-SUB-MAP                    ***
C***    MAP(4)  L-VALUE FOR THIS SUB-SUB-MAP; MAP(4)=L*IHEAD(5)      ***
C***    MAP(5)  LOGARITHM OF FLUX AT EQUATOR; MAP(5)=LOG(F0)*IHEAD(7)***
C***      MAP(6)  =(B1-B0)/B0; B1 IS THE MAGNETIC FIELD STRENGTH     ***
C***              THAT CORRESPONDS TO LOG(F1)=LOG(F0)-1/IHEAD(2)     ***
C***      MAP(7)  =(B2-B1)/B0; LOG(F2)=LOG(F1)-1/IHEAD(2)            ***
C***       ...              ....                                     ***
C***      MAP(L)  LAST ELEMENT IN SUB-SUB-MAP; L=MAP(3)+2            ***
C***    MAP(I)  NUMBER OF ELEMENTS IN NEXT SUB-SUB-MAP; I=L+1        ***
C***       ...              ....                                     ***
C***     ...                ....                                     ***
C***    MAP( )  NUMBER OF ELEMENTS IN LAST SUB-SUB-MAP               ***
C***       ...              ....                                     ***
C***      MAP(K)  LAST ELEMENT IN SUB-MAP; K=MAP(1)                  ***
C***  MAP(J)  NUMBER OF ELEMENTS IN NEXT SUB-MAP; J=MAP(1)+1         ***
C***     ...                ....                                     ***
C***       ...              ....                                     ***
C***   ...                  ....                                     ***   
C***  MAP( )  NUMBER OF ELEMENTS IN LAST SUB-MAP                     ***
C***     ...                ....                                     ***
C***       ...              ....                                     ***
C***      MAP(M)  LAST ELEMENT OF MAP; M=IHEAD(8)                    ***
C***********************************************************************
C**                        ENERGY/MEV GRID                           ***
C**  AE-8:    0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
C**           3.5   4.0   4.5   5.0   5.5   6.0  6.5  7.0 (18 GRID P)***
C**  AE-5,6:  0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
C**           4.0   4.5   (12 GRID POINTS)                           ***
C**  AE-4:    0.04  0.1   0.3   0.5   1.0   2.0  2.5  3.0  3.5  4.0  ***
C**           4.1   4.25  4.35  4.5   4.65  4.85 (16 GRID POINTS)    ***
C**                                                                  ***
C**                         L-VALUE GRID                             ***
C**           BEGIN  STEP  END   STEP  END   STEP  END   STEP  END   ***
C**  AE-8:     1.2   0.05  1.5    0.1  2.0    0.2  2.4    0.1  3.0   ***
C**            3.0   0.2   3.4    0.1  3.6    0.2  4.4    0.1  4.6   ***
C**            4.6   0.2   5.0    0.5* 8.0    1.0 12.0 (43 GRID P.)  ***
C**                    * 6.6 INSTEAD OF 6.5                          ***
C**  AE-5,6:   1.2   0.05  1.5    0.1  2.0    0.2  2.8 (15 GRID P.)  ***
C**  AE-4:     2.8   0.2   4.0    0.5  6.0    0.6  6.6    0.4  7.0   ***
C**            7.0   1.0  11.0   (16 GRID POINTS)                    ***
C***********************************************************************
      INTEGER		AGNR,EGNR
      DIMENSION    	XL(2,10),E(7),FLUX(7),AF(7,10,10),
     &			IHEAD(8),MAP(20000),DF(5,10,10),EDA(3)
      LOGICAL		NOTBEG
      CHARACTER*1	ITE(5,10,10)
      CHARACTER*4	BLTEX,LBTEX
      CHARACTER*5	UNTEX
      CHARACTER*6  	NAME,MNAME(8)
      CHARACTER*9	PARTICLE
      CHARACTER*10 	FNAME
      CHARACTER*12 	IDTEX
      DATA MNAME	/'AP8MAX','AP8MIN','AE4MAX','AE4MIN','AEI7HI',
     &  		'AEI7LO','AE8MAX','AE8MIN'/
      DATA MTYPE,NE,EBEG,EEND,ESTEP,NL,VBEG,VEND,VSTEP,NB,BBEG,BEND,
     &	   BSTEP,JAGNR,INDE,JTAB
     &			/7,0,1.,6.,1.,0,1.,10.,1.,0,1,100,10,0,1,1/
      DATA E,XL,AF,DF	/1227*0.0/,EDA/.05,.1,.2/
C
	DO 1567 I=1,5
	DO 1567 K=1,10
	DO 1567 L=1,10
1567	ITE(I,K,L)=' '
C
C		I/O UNIT NUMBERS
C
	EGNR=5			! INPUT
	MONITO=6		! MONITOR
	IUAEAP=15		! MODEL COEFFICIENTS INPUT
	IUOUT=16		! OUTPUT (OUTPUT.AEP)
C
C		INPUT OF PARAMETERS
C
	NOTBEG=.FALSE.
C-----------1. window: introduction--------------------------------
	WRITE(MONITO,5100)
5100	FORMAT(1X/////8X,'**** AE-8, AP-8 TRAPPED PARTICLE MODELS',
     &	 ' ****'//4X,
     &	 'This program allows you to produce tables of integral'/4X, 
     &   'or differential particle flux for up to 6 energies'/4X,
     &   'varying in L-value or normalized magnetic field strength.'/4X,
     &   'In each of the following windows you will be asked to'/4X,
     &   'enter one or more parameters, defining the conditions'/4X,
     &   'for your model tables.'/4X,
     &   'In each window the current value(s) is (are) displayed'/4X,
     &   'in the right upper corner. You can choose the current'/4X,
     &   'value(s) by entering / at the prompt.'/4X,
     &   'If you enter a wrong character or a value outside the'/4X,
     &   'allowed parameter range, the program will ask you for'/4X,
     &   'a new entry.'/4X,
     &   'After your tables are displayed you can exit or'/4X,
     &   'calculate tables for new input parameters.'/1X,26('*'),
     &   ' GOOD LUCK ',26('*'))
	GOTO 1092
C--------------------2. window: which parameter to change-----------
7001	WRITE(MONITO,7002) 
7002	FORMAT(1X/////' WHICH PARAMETER DO YOU WANT TO CHANGE ?'/
     &	  1X,60('-')/'     0  NO MORE CHANGES, CALCULATE PROFILES'/
     &	  '     1  STORE OR DISPLAY',13X,'4  L-VALUES'/
     &	  '     2  MODEL',24X,'5  B/B0'/
     &    '     3  ENERGIES'/1X,60('-')/' ENTER NUMBER') 
	MAXI=5
	MINI=0
7004	READ(EGNR,*,ERR=7005,END=290) JPARA
	IF((JPARA.LE.MAXI).AND.(JPARA.GE.MINI)) GOTO 7006
7005	WRITE(MONITO,1804) MINI,MAXI
	GOTO 7004
7006	GOTO (7014,1092,7010,7011,7012,7013) JPARA+1
c---------------3. window: display options--------------------------
1092	WRITE(MONITO,6123) JAGNR
6123	FORMAT(1X/' DO YOU WANT YOUR TABLES',34X,'#',I1,'#'/5X,
     &	  'DISPLAYED ON YOUR MONITOR:  ENTER  0  AT PROMPT'/5X,
     &    'STORED IN FILE OUTPUT.AEP:  ENTER  1  AT PROMPT'/5X,
     &    'DISPLAYED AND STORED:       ENTER  2  AT PROMPT')
	MAXI=2
	MINI=0
	WRITE(MONITO,8634)
8105	READ(EGNR,*,ERR=8106,END=290) JAGNR
	IF((JAGNR.LE.MAXI).AND.(JAGNR.GE.MINI)) GOTO 8107
8106	WRITE(MONITO,1804) MINI,MAXI
	GOTO 8105
8107	AGNR=MONITO
	IF(JAGNR.GT.0) THEN
	  OPEN(UNIT=IUOUT,FILE='OUTPUT.AEP',STATUS='NEW',
     &		FORM='FORMATTED')
	  IF(JAGNR.EQ.1) AGNR=IUOUT
	  ENDIF
	IF(NOTBEG) GOTO 7001
c---------------4. window: model type-------------------------------
7010	WRITE(MONITO,5555) MTYPE
5555	FORMAT(1X/' MTYPE',18X,'* choose model *',18X,'#',I1,'#'/
     &	  1X,60('-')/
     &    ' 1   AP8MAX, protons, solar maximum'/' 2   AP8MIN,',
     &    ' protons, solar minimum'/' 3   AE4MAX, electrons,',
     &    ' solar maximum'/' 4   AE4MIN, electrons, solar minimum'/
     &    ' 5   AEI7HI, electrons, high estimate'/' 6   AEI7LO,',
     &    ' electrons, low estimate'/' 7   AE8MAX, electrons, solar',
     &    ' maximum'/' 8   AE8MIN, electrons, solar minimum'//
     &    ' You must have the corrosponding binary data file in'/
     &    ' your account e.g., AE8MAX.BIN, if you enter MTYPE=7.')
      	WRITE(MONITO,8634)
	MAXI=8
	MINI=1
1820	READ(EGNR,*,ERR=1821,END=290) MTYPE
	IF((MTYPE.GE.MINI).AND.(MTYPE.LE.MAXI)) GOTO 1823
1821    	WRITE(MONITO,1804) MINI,MAXI
		GOTO 1820
1822	WRITE(MONITO,1824) FNAME
1824	FORMAT(' Data file ',A10,' is not in your directory,',
     &		' try other data file')
	GOTO 1820
1823	NAME=MNAME(MTYPE)
	WRITE(FNAME,1129) NAME
1129	FORMAT(A6,'.BIN')
	OPEN(IUAEAP,FILE=FNAME,STATUS='OLD',ERR=1822,FORM='UNFORMATTED')
	READ(IUAEAP) IHEAD
	NMAP=IHEAD(8)
	READ(IUAEAP) (MAP(I),I=1,NMAP)

C ASCIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C When using the ASCII coefficient files instead of the binary 
C coefficient files, one should replace the preceding 5 statements
C with the following 6 statements
C
C1129	FORMAT(A6,'.ASC')
C	OPEN(IUAEAP,FILE=FNAME,STATUS='OLD',ERR=1822,FORM='FORMATTED')
C	READ(IUAEAP,1301) IHEAD
C	NMAP=IHEAD(8)
C	READ(IUAEAP,1301) (MAP(I),I=1,NMAP)
C1301	FORMAT(1X,12I6)

 	CLOSE(IUAEAP)
	IF(MTYPE.LT.3) THEN
		PARTICLE='PROTONS'
	ELSE
		PARTICLE='ELECTRONS'
	ENDIF
	IF(NOTBEG) GOTO 7001
C-------------------5. window: number of energies----------------------
7011  WRITE(MONITO,5111) NE    
5111  FORMAT(1X/////' NE',10X,'* Number of Energies *',8X,
     &	'<max=6>',7X,'#',I2,'#'/1X,60('-')/
     &  ' If you enter 0, you will be asked for the BEGIN, END'/
     &  ' and STEPSIZE of your energy grid.')
      WRITE(MONITO,8634)
      MAXI=6
      MINI=0
1801  READ(EGNR,*,ERR=1802,END=290) NE
      IF((NE.GE.MINI).AND.(NE.LE.MAXI)) GOTO 1803
1802  WRITE(MONITO,1804) MINI,MAXI
1804  FORMAT(' Your input is ouside the value range:',I4,' to',I4,
     & 	', try again')
      GOTO 1801
1803	XMAX=7.
	XMIN=0.1
	IF(MTYPE.LT.3) XMAX=400
	IF(MTYPE.GT.2) XMIN=0.04
      IF(NE.GT.0) GOTO 1899
c------------------6. window (a): input energy grid-------------------
	WRITE(MONITO,511) EBEG,EEND,ESTEP
511	FORMAT(1X//' BEGIN, END, STEPSIZE    * MeV *     #',
     &	 F6.2,', ',F6.2,', ',F6.2,'#'/1X,60('-')/3X,
     &	 'Begin, end, and stepsize of energy grid in MeV.'/3X,
     &   'The number of grid points is restricted to 10 or less !')
	WRITE(MONITO,8634)
1810	READ(EGNR,*,ERR=1811,END=290) EBEG,EEND,ESTEP
	NE=INT((EEND-EBEG)/ESTEP)+1
	IF(NE.GT.6) NE=6
	IF((EBEG.GE.XMIN).AND.(EEND.LE.XMAX)) GOTO 1812
1811	  WRITE(MONITO,1813) XMIN,XMAX
1813  	  FORMAT(' Your input is ouside the value range:',F8.2,
     &	   ' to ',F8.2,', try again')
	  GOTO 1810
1812	E(1)=EBEG
	DO 512 IE=2,NE
512		E(IE)=E(IE-1)+ESTEP
    	GOTO 1898
c--------------------6. window (b): input energies-----------------
1899	WRITE(MONITO,513) NE,(E(J),J=1,6)
513	FORMAT(1X//' E(1), ... E(',I1,')     * MeV *'//' #',
     &	  F7.2,5(',',F7.2),'#'/1X,60('-')/
     &	  ' Enter energies in ascending order.')
	WRITE(MONITO,8634)
1816	READ(EGNR,*,ERR=1817,END=290) (E(IE),IE=1,NE)
	IF((E(1).GE.XMIN).AND.(E(NE).LE.XMAX)) GOTO 1898
1817	WRITE(MONITO,1813) XMIN,XMAX
	GOTO 1816
1898	DO 1779 EI=NE+1,6
1779		E(EI)=0.0
	IF(NOTBEG) GOTO 7001
C-------------------7. window: number of L-values--------------------
7012  WRITE(MONITO,5311) NL
5311  FORMAT(1X///' NL',10X,'* Number of L-values *',9X,'<max=10>',
     &	4X,'#',I3,'#'/1X,60('-')/
     & 	' If you enter 0, you will be asked for the BEGIN, END'/
     &  ' and STEPSIZE of your L-value grid.')
      WRITE(MONITO,8634)
      MAXI=10
      MINI=0
2801  READ(EGNR,*,ERR=2802,END=290) NL
      IF((NL.GE.MINI).AND.(NL.LE.MAXI)) GOTO 2803
2802  WRITE(MONITO,1804) MINI,MAXI
      GOTO 2801
2803  	XMIN=1.
	XMAX=15.
	IF(NL.GT.0) GOTO 2899
	WRITE(MONITO,3111) VBEG,VEND,VSTEP
c----------------8. window (a): L grid---------------------------------
3111	FORMAT(1X//' BEGIN, END, STEPSIZE     * L *      #',
     &	 F6.2,', ',F6.2,', ',F6.2,'#'/1X,60('-')/3X,
     &	 'Begin, end, and stepsize of L-value grid.'/3X,
     &   'The number of grid points is restricted to 10 or less !')
	WRITE(MONITO,8634)
2810	READ(EGNR,*,ERR=2811,END=290) VBEG,VEND,VSTEP
	IF((VBEG.GE.XMIN).AND.(VEND.LE.XMAX)) GOTO 2812
2811	  WRITE(MONITO,1813) XMIN,XMAX
	  GOTO 2810
2812	NL=INT((VEND-VBEG)/VSTEP)+1
	IF(NL.GT.10) NL=10
	XL(1,1)=VBEG
	DO 5122 IE=2,NL
5122		XL(1,IE)=XL(1,IE-1)+VSTEP
	GOTO 2898
C----------------8. window (b): L-values-------------------------
2899	WRITE(MONITO,5133) NL,(XL(1,J),J=1,10)
5133	FORMAT(1X//' L(1), ... L(',I2,')     * L-values *'//
     &    ' #',F7.2,5(',',F7.2)/2X,F7.2,3(',',F7.2),'#')
	WRITE(MONITO,8634)
2816	READ(EGNR,*,ERR=2817,END=290) (XL(1,IE),IE=1,NL)
	IF((XL(1,1).GE.XMIN).AND.(XL(1,NL).LE.XMAX)) GOTO 2898
2817	WRITE(MONITO,1813) XMIN,XMAX
	GOTO 2816
2898	DO 2798 IE=NL+1,10
2798		XL(1,IE)=0.0
	IF(NOTBEG) GOTO 7001
C-----------------9. window: number of B/B0 values------------------
7013  WRITE(MONITO,5551) NB    
5551  FORMAT(1X/////' NB',9X,'* Number of B/B0 values *',8X,
     &	'<max=10>',3X,'#',I3,'#'/1X,60('-')/
     & 	' If you enter 0, you will be asked for the BEGIN, END'/
     &  ' and STEPSIZE of your B/B0 value grid.')
      WRITE(MONITO,8634)
      MAXI=10
      MINI=0
3801  READ(EGNR,*,ERR=3802,END=290) NB
      IF((NB.GE.MINI).AND.(NB.LE.MAXI)) GOTO 3803
3802  WRITE(MONITO,1804) MINI,MAXI
      GOTO 3801
3803	XMAX=9999.9
	XMIN=1.
      IF(NB.GT.0) GOTO 3899
C------------10. window (a): B/B0 grid--------------------------------
	WRITE(MONITO,5113) BBEG,BEND,BSTEP
5113	FORMAT(1X//' BEGIN, END, STEPSIZE    * B/B0 *    #',
     &	 F7.2,',',F7.2,',',F7.2,'#'/1X,60('-')/3X,
     &	 'Begin, end, and stepsize of B/B0 grid.'/3X,
     &   'The number of grid points is restricted to 10 or less !')
	WRITE(MONITO,8634)
3810	READ(EGNR,*,ERR=3811,END=290) BBEG,BEND,BSTEP
	NB=INT((BEND-BBEG)/BSTEP)+1
	IF(NB.GT.10) NB=10
	IF((BBEG.GE.XMIN).AND.(BEND.LE.XMAX)) GOTO 3812
3811	  WRITE(MONITO,1813) XMIN,XMAX
	  GOTO 3810
3812	XL(2,1)=BBEG
	DO 5123 IE=2,NB
5123		XL(2,IE)=XL(2,IE-1)+BSTEP
    	GOTO 3898
C-----------10. window (b): B/B0 values------------------------
3899	WRITE(MONITO,6133) NB,(XL(2,J),J=1,10)
6133	FORMAT(1X//' B/B0(1), ... B/B0(',I2,')'//
     &		' #',F7.2,5(',',F7.2)/2X,F7.2,3(',',F7.2),'#')
	WRITE(MONITO,8634)
3816	READ(EGNR,*,ERR=3817,END=290) (XL(2,IE),IE=1,NB)
	IF((XL(2,1).GE.XMIN).AND.(XL(2,NB).LE.XMAX)) GOTO 3898
3817	WRITE(MONITO,1813) XMIN,XMAX
	GOTO 3816
3898	DO 3798 IE=NB+1,10
3798		XL(2,IE)=0.0	
	IF(NOTBEG) GOTO 7001
	NOTBEG=.TRUE.
C----------------THE L-VALUE LOOP-------------------------------
7014  DO 280 L=1,NL  
      	FL=XL(1,L) 
C----------------THE B LOOP------------------------------------- 
      	DO 160 I=1,NB
		BB0=XL(2,I)
      		CALL TRARA1(IHEAD,MAP,FL,BB0,E,FLUX,NE)
C----------------THE ENERGY LOOP--------------------------------
		DO 1110 K=1,NE
			AF(K,L,I)=0.0
	 		IF(FLUX(K).GT.0.0) AF(K,L,I)=10.**FLUX(K) 
			IF(K.EQ.1) GOTO 1110
			DF(K-1,L,I)=ABS(AF(K,L,I)-AF(K-1,L,I))/
     &				(E(K)-E(K-1))
			IF(AF(K,L,I).LE.0.0) DF(K-1,L,I)=0.0
1110			CONTINUE
160  		CONTINUE            
280	CONTINUE
C-----------------TESTS VALIDITY OF DIFFERENTIAL FLUX------------
	DO 7788 K=1,NE-1
		EI=E(K+1)
		ITEST=0
    		EDIFF=EI-E(K)     
C-----------------IS ENERGY INTERVALL LARGE ENOUGH ?-------------
C      		IF ((EI.GT.0.25).AND.(EDIFF.LT.0.1999)) ITEST=1
C      		IF ((EI.LE.0.10).AND.(EDIFF.LT.0.0499)) ITEST=1
C      		IF(((EI.LE.0.25).AND.(EDIFF.LT.0.0999)
C     &				.AND.(EI.GT.0.1))) ITEST=1
		IEI=1      
      		IF(EI.GT.0.10) IEI=2
      		IF(EI.GT.0.25) IEI=3
      		IF(EDIFF.LT.EDA(IEI)) ITEST=1
		DO 7788 L=1,NL
			ITT=0
	      		IF(XL(1,L).LT.1.2) ITT=1
			DO 7788 I=1,NB
      				ITE(K,L,I)=' '
				IF(ITEST+ITT.NE.0) ITE(K,L,I)='?'
7788				CONTINUE
C------------------TABLE OUTPUT------------------------------------
	WRITE(MONITO,9101)
9101	FORMAT(1X/)
1097	WRITE(MONITO,9000) JTAB
9000	FORMAT(' TABLE OPTIONS:',6X,'integral',12X,'differential flux',
     &	  5X,'#',I1,'#'/21X,40('-')/
     & 	  2X,'0  EXIT            1  E * L TABLE      3  E * L TABLE'/
     &    2X,'5  NEW PARAMETER   2  E * B/B0 TABLE   4  E * B/B0 TABLE')
	WRITE(MONITO,8634)
8735	FORMAT(1X,70('-')/,' enter /, to continue with current value(s);',
     &	  ' Ctrl Z, to exit')
8634	FORMAT(1X,60('-')/,' enter /, to continue with current value(s);',
     &	  ' Ctrl Z, to exit')
	MINI=0
	MAXI=5
9003	READ(EGNR,*,ERR=9001,END=290) JTAB
	IF((JTAB.LE.MAXI).AND.(JTAB.GE.MINI)) GOTO 9002
9001	WRITE(MONITO,1804) MINI,MAXI
	GOTO 9003
9002	IF(JTAB.EQ.0) GOTO 290
	IF(JTAB.EQ.5) GOTO 7001
	IF(JTAB.GT.2) THEN
		UNTEX='*MeV]'
		IDTEX='differential'
	ELSE
		UNTEX=']    '
		IDTEX='integral'
	ENDIF
	IF((JTAB.EQ.1).OR.(JTAB.EQ.3)) THEN
		IBL=2
		BLTEX='B/B0'
		LBTEX=' L  '
		N=NL
		NN=NB
	ELSE
		IBL=1
		BLTEX=' L  '
		LBTEX='B/B0'
		N=NB
		NN=NL
	ENDIF
	IBLTAB=1
	IF(IBL.EQ.1) IBLTAB=2
	IF(NN.EQ.1) GOTO 9112
	WRITE(MONITO,9104) BLTEX,INDE,BLTEX,(XL(IBL,J),J=1,10)
9104	FORMAT(1X/' FOR WHICH ',A4,' ?         * enter index number,',
     &	  ' not value*',9X,'#',I2,'#'/1X,70('-')/' Nr:',5X,
     &	  '1     2     3     4     5     6     7   ',
     &	  '  8     9    10'/1X,A4,':',10F6.1)
	WRITE(MONITO,8735)
	MINI=1
	MAXI=10
9113	READ(EGNR,*,ERR=9111,END=290) INDE
	IF((INDE.LE.MAXI).AND.(INDE.GE.MINI)) GOTO 9112
9111	WRITE(MONITO,1804) MINI,MAXI
	GOTO 9113
9112    WRITE(AGNR,1091) IDTEX,PARTICLE,UNTEX,LBTEX,(E(J),J=1,6)
	IF(JAGNR.EQ.2) 
     &	 WRITE(IUOUT,1091) IDTEX,PARTICLE,UNTEX,LBTEX,(E(J),J=1,6)
1091	FORMAT(1X/25X,A12,' flux [',A9,'/cm*cm*sec',A5/
     &	  1X,A4,'\\ E/MeV',F7.2,5F10.2/1X,'------\\',62('-'))
	BLV=XL(IBL,INDE)
	DO 1094 I=1,N
	IL=INDE
	IB=INDE          
	IF(IBL.EQ.1) THEN
		IB=I
		XX=XL(2,I)
	ELSE
		IL=I
		XX=XL(1,I)
	ENDIF
	IF(JTAB.LT.3) THEN
		WRITE(AGNR,3092) XX,(AF(J,IL,IB),J=1,6)
		IF(JAGNR.EQ.2) WRITE(IUOUT,3092) XX,(AF(J,IL,IB),J=1,6)
	ELSE
		WRITE(AGNR,1093) XX,DF(1,IL,IB),ITE(1,IL,IB),
     &		 DF(2,IL,IB),ITE(2,IL,IB),DF(3,IL,IB),ITE(3,IL,IB),
     &		 DF(4,IL,IB),ITE(4,IL,IB),DF(5,IL,IB),ITE(5,IL,IB)
		IF(JAGNR.EQ.2) 
     &		 WRITE(IUOUT,1093) XX,DF(1,IL,IB),ITE(1,IL,IB),
     &		 DF(2,IL,IB),ITE(2,IL,IB),DF(3,IL,IB),ITE(3,IL,IB),
     &		 DF(4,IL,IB),ITE(4,IL,IB),DF(5,IL,IB),ITE(5,IL,IB)
	ENDIF
3092	FORMAT(1X,F7.2,' I  ',1PE9.3,5(1X,1PE9.3))
1093	FORMAT(1X,F7.2,' I  ',5X,5(1PE9.2,A1))
1094	CONTINUE
	WRITE(AGNR,1095) NAME,BLTEX,BLV
	IF(JAGNR.EQ.2) WRITE(IUOUT,1095) NAME,BLTEX,BLV
1095	FORMAT(1X,20('-'),' flux values below 10 are not reliable ',
     &	  11('-')/1X,'Model: ',A6,10X,A4,'=',F8.2/1X,70('*'))
	GOTO 1097
 290  CONTINUE 
      STOP     
      END
