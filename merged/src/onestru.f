	subroutine onestru
 	implicit real*8(a-h,o-z)
	real*8 masch2,masbo2,masto2,lambda2
	CHARACTER*20 PARM(20)
	character*10 projectile(6)
	CHARACTER*10 CIBLE 
	COMMON/OUTPUT_UNITS/io,ier,itest
        COMMON/ALFA/MASCH2,MASBO2,MASTO2,LAMBDA2
        DIMENSION CALC(8,20,25),CALCPI(8,20,25)
     1,PAR(40),PARPI(40)
	DIMENSION VALUE(20)
       include "conv.f"
       include "projec.f"
       include "qcd.f"
       include "phoscale.f"
c
c *** common blocks for use with pdflib tables
      COMMON/W50511/NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
      COMMON/W50512/QCDL4,QCDL5
      COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX
      include "pdfmod.f"
c   structure functions conventions
c   iord=0                     leading order
c   iord=1             next to leading order
c   iconv=0             universal convention
c   iconv=1          quark distribution = f2
c
c        IORD=1,iconv=0 standard choice
c
         data projectile/'proton','antiproton','pi+','pi-',
     1    ' ','nucleus'/
c  ------- choice of structure functions
c

c
 	
	save
      write(io,102)
102   FORMAT(' ')
c      IF(ISTRUC.EQ.1) write(io,*) ' STRUCTURE FUNCTIONS=FROM MRS98 '
c      IF(ISTRUC.EQ.2) write(io,*) ' STRUCTURE FUNCTIONS=FROM MRS99'      
c      IF(ISTRUC.EQ.3) write(io,*) ' STRUCTURE FUNCTIONS FROM CTEQ5M'
      IF(ISTRUC.EQ.4) write(io,*) ' STRUCTURE FUNCTIONS FROM PDFLIB'
      IF(ISTRUC.EQ.5) write(io,*) ' STRUCTURE FUNCTIONS=FROM TABLES '


C
C  ----->>>>>>>>>>>>>>>>> INIT PDFLIB <<<<<<<<<<<<<<<<<<<<<< ------
C
      IF(ISTRUC.EQ.4) THEN
          PARM(1)='NPTYPE'
          VALUE(1)=NPTYPRO
          PARM(2)='NGROUP'
          VALUE(2)=NGROPRO
          PARM(3)='NSET'
          VALUE(3)=NSETPRO
cc          PARM(4)='NFL'
cc          VALUE(4)=NF
C         PARM(6)='HO'
C         VALUE(6)=1
cc         PARM(5)='TMAS'
c          VALUE(5)=1.D10
cc	  VALUE(5)=DSQRT(MASTO2)
c         PARM(9)='QCDL4'
c         PARM(10)='XMIN'
c         PARM(11)='XMAX'
c         PARM(12)='Q2MIN'
c         PARM(13)='Q2MAX'
c setting of lambda4 for cteq4   as it does not work correctly
	  if(value(2).eq.4) then
	  value(4)=0.D0
CCTEQ4M=34		
	   if(value(3).eq.34) then
	   PARM(4)='QCDL4'
	   value(4)=0.296D0
	   
CCTEQ4A1=35
	   elseif(value(3).eq.35) then
	   PARM(4)='QCDL4'
	   value(4)=0.213D0
CCTEQ4A2=36	
	   elseif(value(3).eq.36) then
	   PARM(4)='QCDL4'
	   value(4)=0.253D0
CCTEQ4A4?=38	
	   elseif(value(3).eq.38) then
	   PARM(4)='QCDL4'
	   value(4)=0.344D0
CCTEQ4A5?=39
	   elseif(value(3).eq.39) then
	   PARM(4)='QCDL4'
	   value(4)=0.399D0
CCTEQ4HJ=40	
	   elseif(value(3).eq.40) then
	   PARM(4)='QCDL4'
	   value(4)=0.302D0
	  endif
	  endif
c no nuclear effect
          if(value(1).ne.4) then 
	  CALL PDFSET(PARM,VALUE)
	  OWLAM=QCDL4
c nuclear effect	  
	  elseif(value(1).eq.4) then
	  owlam=0.0d0
	  if(value(2).eq.4 .and.value(4).gt.0.1D0) owlam=value(4)         
	  PARM(1)='NPTYPE'
          VALUE(1)=1
          PARM(2)='NGROUP'
          VALUE(2)=NGROPRO
          PARM(3)='NSET'
          VALUE(3)=NSETPRO
	 PARM(4)='NATYPE'
	 VALUE(4)=4
	 PARM(5)='NAGROUP'
	 VALUE(5)=1
	 PARM(6)='NASET'
	 VALUE(6)=1
	 CALL PDFSET(PARM,VALUE)
	   if(owlam.lt.0.1d0) owlam=qcdl4
	  endif
          OWLAM2=OWLAM**2
	  
          ALBD=OWLAM
          ALBD2=OWLAM2
	  lambda2 = albd2
	  
          Q0S=Q2MIN
          Q02PR=Q0S
          Q02PI=Q0S
      write(io,*) ' THE PROTON STRUCT FUNCT. ARE FROM THE PDFLIB TABLES:
     1 NGROUP = ',NGROTAR,' NSET =',NSETTAR

      write(io,105) QCDL4
105   FORMAT('  THE PROTON : LAMBDA=',F7.3 )
c
c     the pion structure functions from pdflib
c
          IF(IPRO.GE.7.OR.IPRO.EQ.3.OR.IPRO.EQ.4)  THEN
          VALUE(1)=NPTYPRO
          VALUE(2)=NGROPRO
          VALUE(3)=NSETPRO
          CALL PDFSET(PARM,VALUE)
c         NPTYPRO=VALUE(1)
c          NGROPRO=VALUE(2)
c          NSETPRO=VALUE(3)
          Q02PI=Q0S
          write(io,102)
      write(io,*) ' THE PION   STRUCT FUNCT. ARE FROM THE PDFLIB TABLES:
     1 NGROUP = ',NGROPRO,' NSET ',NSETPRO
          write(io,106) QCDL4
106   FORMAT('  THE PION   : LAMBDA=',F7.3 )
           write(io,102)
           LAMDEL=ABS(OWLAM-QCDL4)
           Q0SDEL=ABS(Q0S-Q2MIN)
           IF(LAMDEL.GT.1.D-3) write(io,107)
           IF(Q0SDEL.GT.1.D-3)  write(io,108)
107   FORMAT(' LAMBDA VALUES NOT COMPATIBLE IN PION AND NUCL PDFLIB')
108   FORMAT(' Q2MIN  VALUES NOT COMPATIBLE IN PION AND NUCL PDFLIB')
           ENDIF
      ENDIF
C
C  ----->>>>>>>>>>>>>>>>> INIT D.I.T. <<<<<<<<<<<<<<<<<<<<<< ------
C
      IF(ISTRUC.EQ.5) THEN
c
c ----- input structure function tables for the nucleon
c      OPEN(UNIT=12,NAME='MSL23G40.DREYAN',TYPE='OLD',READONLY)
       OPEN(UNIT=12,FILE='MSL23G40.DREYAN',STATUS='OLD')
c ----- input structure function tables for the pion
       OPEN(UNIT=22,FILE='LPIOG194.DREYAN',STATUS='OLD')
        IF(IORD.EQ.0) THEN
          READ(11,2) PAR
          READ(11,2) CALC
          IF(IPRO.GE.7.OR.IPRO.EQ.3.OR.IPRO.EQ.4)
     1    READ(21,2) PARPI
          IF(IPRO.GE.7.OR.IPRO.EQ.3.OR.IPRO.EQ.4)
     1    READ(21,2) CALCPI
        ELSE
          READ(12,2) PAR
          READ(12,2) CALC
          IF(IPRO.GE.7.OR.IPRO.EQ.3.OR.IPRO.EQ.4)
     1    READ(22,2) PARPI
          IF(IPRO.GE.7.OR.IPRO.EQ.3.OR.IPRO.EQ.4)
     1    READ(22,2) CALCPI
        ENDIF
2       FORMAT(8D15.8)
C
        OWLAM=PAR(1)
        OWLAM2=OWLAM**2
        ALBD2=OWLAM2
        ALBD=PAR(1)
C
      write(io,*) ' THE PARAMETERS OF THE TABLE PARTON DISTRIBUTIONS'
      write(io,103) PAR(1),PAR(13),PAR(10)
      IF(IPRO.GE.7.OR.IPRO.EQ.3.OR.IPRO.EQ.4)
     1write(io,104) PARPI(1),PARPI(10),PARPI(8),PARPI(2),PARPI(3)
103   FORMAT('  THE PROTON: LAMBDA=',F7.3,', ETA-GL.=',F7.3
     1,', ETA-SEA=',F7.3)
104   FORMAT('  THE PION  : LAMBDA=',F7.3,', ETA-GL.=',F7.3
     1,', ETA-S=',F7.3,', ETA-1=',F7.3,', ETA-2=',F7.3)
C   VALUE OF QO**2 IN STRUCTURE FUNCTIONS OF PIONS
        Q02PR=2.D0
        Q02PI=Q02PR
      ENDIF
C
C  ----->>>>>>>>>>>>>>>>> INIT MRS98/MRS99<<<<<<<<<<<<<<<<<<<<<< ------
C
      IF(ISTRUC.EQ.1.or.ISTRUC.EQ.2) THEN
       mmrs98=nsetpro
      if(mmrs98.eq.1.or.mmrs98.eq.2.or.mmrs98.eq.3) then
       ALBD=0.300D0
      elseif(mmrs98.eq.4) then
       ALBD=0.229D0
      elseif(mmrs98.eq.5) then
       ALBD=0.383D0
      endif
C   VALUE OF LAMBDA AND QO**2 IN STRUCTURE FUNCTIONS
C
        ALBD2=ALBD**2
        Q0S=1.25d0
	lambda2=albd2

        OWLAM=ALBD
        OWLAM2=OWLAM**2

        Q02PR=Q0S
        Q02PI=Q02PR	
      write(io,*) ' THE PROTON STRUCT FUNCT. ARE FROM MRS',
     1ngropro,' PACKAGE:',
     1 'MODE = ',mmrs98,' LAMBDA4 =',albd,' GeV'
c      write(io,105) QCDL4
c105   FORMAT('  THE PROTON : LAMBDA=',F7.3 )
      ENDIF
C
C  ----->>>>>>>>>>>>>>>>> INIT CTEQ5M <<<<<<<<<<<<<<<<<<<<<< ------
C
      IF(ISTRUC.EQ.3) THEN
        iset=nsetpro
c	call SetCtq5(Iset)
	albd=0.326d0
	if(Iset.eq.3) then
	  albd=0.192D0
	elseif(Iset.eq.7) then
	  albd=0.309D0
	endif  
   
C   VALUE OF LAMBDA AND QO**2 IN STRUCTURE FUNCTIONS
C
        ALBD2=ALBD**2
        Q0S=1.25d0
	
	lambda2 = albd2


        OWLAM=ALBD
        OWLAM2=OWLAM**2

        Q02PR=Q0S
        Q02PI=Q02PR	
      write(io,*) ' THE PROTON STRUCT FUNCT. ARE FROM CTEQ5M PACKAGE:',
     1 'MODE = ',Iset,' LAMBDA4 =',albd
c      write(io,105) QCDL4
c105   FORMAT('  THE PROTON : LAMBDA=',F7.3 )
 	

      ENDIF
C
C
C  ----->>>>>>>>>>>>>>>>> END OF INIT S.F. <<<<<<<<<<<<<<<<<<<<<< ------
C
C     VALUE OF LAMBDA IN ALFAS
      RLAM=OWLAM
      RLAM2=RLAM**2
C
      IF(ACIB.eq.1D0 .and. ZCIB.eq.1D0) then
        CIBLE='PROTON'
      ELSE
        CIBLE='NUCLEUS'
      ENDIF		

      write(io,101)
101   FORMAT(//)
       write(io,*) '   PROJECTILE=',projectile(ipro)
       if(ipro.eq.6) then
       write(io,*)' with  Z= ',ZPRO,' A= ',APRO
       write(io,*) '  cross sections are given per nucleon without ',
     1 'corrections for nuclear effects' 
       endif
       if (CIBLE.EQ.'PROTON') then
       write(io,*) '   TARGET=',CIBLE
       else		 
       write(io,*) '   TARGET=',CIBLE,' with  Z= ',ZCIB,' A= ',ACIB
       write(io,*) '  cross sections are given per nucleon without ',
     1 'corrections for nuclear effects' 
       endif
      write(io,102)
      write(io,*) ' ORDER OF STR. FUN.=',IORD,'  CONVENTION OF STR. FUN.
     1=',ICONV
      write(io,*) ' LAMBDA IN STR. FUN.=',OWLAM,' GeV LAMBDA IN ALFAS='
     2,RLAM,' GeV'
c      write(io,102)
      write(io,*) ' Q02 VALUES : FOR PROTON Q02=',Q02PR,'GeV^2  FOR PION Q02
     1=',Q02PI,'GeV^2'
      write(io,102)
      RETURN
      END
