	subroutine pioincl(csbor,csect)
C*===================================================================
C
C                            PI0INCL
C                            =======
C
C PHENOMENOLOGY PROGRAM TO CALCULATE HADRON - HADRON CROSS-SECTIONS
C NUMERICALLY, STARTING FROM MATRIX ELEMENTS AT O(ALFAS**3).
C FOR INCLUSIVE PI0 PRODUCTION.
C                      ALL ENERGIES IN GEV
C  returns cross section csbor, born  and csect, NLL 
C  according to data passed via commons
C==================================================================  */
      IMPLICIT REAL*8(A-H,L-Z)
     
C******* COMMUNICATION WITH SRUCTURE FUNCTIONS *********************
      COMMON/PION/IPION
C*******************************************************************
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/HASCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/FACT/ZQ
      COMMON/YRANGE/YMIN,YMAX
      COMMON/PTVALUE/PT
      common/output_units/io,ier,itest
C******* COMMUNICATION WITH  param ********************************
      COMMON/EVOPT/PTEVO(19)
      COMMON/VALU/JMAR,IPT
      COMMON/XSECT/ISIGM,ILOOP,IHIOR
      COMMON/CHOI/HC2,ZRAP,ZRS,ZAL,CM,CMU,CMP
c      COMMON/ABFOW1/IMODE
      COMMON/ALFA/MASCH2,MASBO2,MASTO2,LAMBDA2
      COMMON/EXPERIENCE/IEXP
      include "pdfchoix.f"
      common/INTCHOIX/integy
      common/baprecision/accu,jtmx1,jtmx2,jcall1,jcall2,jcall3
C******** COMMUNICATION WITH bases**********************************
	COMMON/BPARM1/XL(50),XU(50),IDIM,IWILD,IG(50),ICALL
	COMMON/BPARM2/ACC1,ACC2,ITMX1,ITMX2
C*******************************************************************
	logical integy
	data ifirst/0/
	save
c	save ifirst
      EXTERNAL CDEL,DPLUS
      PI=4.*DATAN(1.D0)

      CALL FICT
 
      RAP=ZRAP 
      RS=ZRS
      CQ=ZQ
      AL=ZAL
C
      if(ifirst.eq.0) then
      WRITE (itest,55) CMU,CM,CMP
 55   FORMAT(3X,'MU=PT*',E15.6,2X,'M=PT*',E15.6,2X,'MP=PT*',E15.6,/)
      WRITE (itest,35) RS
 35   FORMAT(3X,'SQUARE ROOT OF S=',E15.6,/)
      endif 
C 
      GS=RS**2
      N=3.D0
      VC=(N)**2-1.D0
      CF=4.D0/3.D0
      V1=VC**2/N
      V2=VC/N
      V3=(N**4-1.D0)/2.D0/N**2
      V4=VC**2/2.D0/N**2

        ETA=2.D0*PT/RS
        PT2=PT**2*CM**2
        pt2 = dmax1(pt2,2.1d0)

cc
        M=DSQRT(PT2)
        MU=CMU*PT
        MU = DMAX1(MU,DSQRT(2.1D0))
        MP=PT*CMP
        MP = DMAX1(MP,DSQRT(2.1D0))
c        GTR=NF/2.D0
        MU2=MU**2

c        WRITE (8,*) '  BORN     BORN+HIGHER ORDER '
C  
C  COMMENT *******************************************************
C        J0=16==>  PROCESSUS : G  G   ---> QJ
C        J0=15==>  PROCESSUS : G  G   ---> G
C        J0=14==>  PROCESSUS : QI G   ---> G
C        J0=13==>  PROCESSUS : QI G   ---> QI
C        J0=12==>  PROCESSUS : QI QBI ---> G
C        J0=11==>  PROCESSUS : QI QBI ---> QI
C        J0=10==>  PROCESSUS : QI G   ---> QBI
C        J0=9 ==>  PROCESSUS : QI G   ---> QBK
C        J0=8 ==>  PROCESSUS : QI G   ---> QK
C        J0=7 ==>  PROCESSUS : QI QI  ---> G
C        J0=6 ==>  PROCESSUS : QI QI  ---> QI
C        J0=5 ==>  PROCESSUS : QI QBI ---> QK
C        J0=4 ==>  PROCESSUS : QI QBK ---> G
C        J0=3 ==>  PROCESSUS : QI QBK ---> QI
C        J0=2 ==>  PROCESSUS : QI QK  ---> G
C        J0=1 ==>  PROCESSUS : QI QK  ---> QI
C*******************************************************************
C      IFLAG=1 (RESP 2) CORRESPONDs to ALPHAS 1 loop (2 loopsS)
C      ICHOI=1 (RESP 2) CORRESPONDs to BORN ( BORN+PART IN DELTA(1-W))
       IFLAG=ILOOP	
       if(ifirst.eq.0) then
	  write(itest,*) ' MU**2 ',MU2,' ALFAS ',ALFAS(IFLAG,MU2)
	  ifirst=1
       endif	  
        ICHOI=1
	CALL BSINIT
	ACC1 = ACCU
	ACC2 = ACCU
	ITMX1 = JTMX1
	ITMX2 = JTMX2
	ICALL = JCALL1
	if (integy) then
	  IDIM = 3
	  IWILD = 3
	else
	  IDIM = 2
	  IWILD = 2
	endif
	DO I=1,IDIM
	  XL(I) = 0.D0
	  XU(I) = 1.D0
	ENDDO
c	write(itest,*)'pioincl before bases cdel jcall1 mp',mp
	CALL BASES (CDEL,RES,SD,CTIME,IT1,IT2)
c	write(itest,*)'pioincl after bases cdel jcal1 mp',mp
        RBOR=RES*HC2/PI/GS
        SDBOR=SD*HC2/PI/GS
C  RBOR EST E*DSIGMA/D3P POUR LE TERME DE BORN
        IF (IHIOR.EQ.2) THEN
          ICHOI=2
	  CALL BSINIT
	  ACC1 = ACCU
	  ACC2 = ACCU
	  ITMX1 = JTMX1
	  ITMX2 = JTMX2
	  ICALL = JCALL2
	  if (integy) then
	    IDIM = 3
	    IWILD = 3
	  else
	    IDIM = 2
	    IWILD = 2
	  endif
	  DO I=1,IDIM
	    XL(I) = 0.D0
	    XU(I) = 1.D0
	  ENDDO
c	  write(itest,*)'pioincl before bases cdel jcall2 mp',mp
	  CALL BASES (CDEL,RES,SD,CTIME,IT1,IT2)
c	  write(itest,*)'pioincl after bases cdel jcall2 mp',mp
	  SECH1=RES
	  SDECH1=SD
c          WRITE (8,*) 'FROM BASIS CDEL SECH1=',SECH1,SDECH1
c
	  CALL BSINIT
	  ACC1 = ACCU
	  ACC2 = ACCU
	  ITMX1 = JTMX1
	  ITMX2 = JTMX2
	  ICALL = JCALL3
	  if (integy) then
	    IDIM = 4
	    IWILD = 4
	  else
	    IDIM = 3
	    IWILD = 3
	  endif
	  DO I=1,IDIM
	    XL(I) = 0.D0
	    XU(I) = 1.D0
	  ENDDO
c	  write(itest,*)'pioincl before bases dplus mp',mp
	  CALL BASES (DPLUS,RES,SD,CTIME,IT1,IT2)
c	  write(itest,*)'pioincl after bases dplus mp',mp
	  SDXES1=SD
          XES1=RES
c          WRITE (8,*) 'FROM BASIS DPLUS XES1=',XES1,SDXES1
c	  WRITE (8,*) 'TIME ',CTIME
	  ERRORNLO=SDECH1+SDXES1
          SECH1=SECH1+XES1
          SIGTO=SECH1*HC2/PI/GS
          ERRORNLO=ERRORNLO*HC2/PI/GS
C SIGTO EST E*DSIGMA/D3P POUR LE TERME  BORN+CORRECTIONS RADIATIVES
        ELSE IF (IHIOR.EQ.1) THEN
          SIGTO=0.D0
        ENDIF
        IF (ISIGM.EQ.1) THEN
C  POUR OBTENIR DSIGMA/DY/DPT2
          FACIN=PI
        ELSE IF (ISIGM.EQ.2) THEN
          FACIN=1.D0
        ELSE IF (ISIGM.EQ.3) THEN
C  POUR OBTENIR DSIGMA/DPT/DY
          FACIN=PI*2.D0*PT
        ENDIF
C
 	if (integy) then
c	  write(8,*)' RBOR*FACIN ',RBOR*FACIN
c	  write(8,*)' SIGTO*FACIN ',SIGTO*FACIN
c	  write(8,*)' ymax-ymin ',ymax-ymin
          CSBOR=FACIN*RBOR/(ymax-ymin)
          CSECT=FACIN*SIGTO/(ymax-ymin)
          ERRORNLO=FACIN*ERRORNLO/(ymax-ymin)
          ERRORLL=FACIN*SDBOR/(ymax-ymin)
	else
          CSBOR=FACIN*RBOR
          CSECT=FACIN*SIGTO
          ERRORNLO=FACIN*ERRORNLO
          ERRORLL=FACIN*SDBOR
	endif
c 
	 write(itest,*) 'PT ',PT,' CSBOR ',CSBOR,' CSECT ',CSECT,
     1   ' CSECT-CSBOR ',CSECT-CSBOR,' ERROR On LL',ERRORLL,
     2	' ERROR On NLO',ERRORNLO 
c         write(8,112) CSBOR,CSECT 
c 112     format(1x,2(D12.5,1x))
c        write(9,113)pt,csbor,csect
c 113    format(1x,E12.5,1x,D12.5,1x,D12.5)
c         WRITE (8,211) CSBOR,CSECT
c         write (9,*) pt,csect-csbor,error
c	 write(8,*) 'PT MP**2 AA XLM2(I2) CSBOR CSECT '
c	 write (8,*) pt,mp**2,aa,xlm2(i2),CSBOR,CSECT

      return
 211  FORMAT(2(2X,D12.6))
 1    FORMAT(3X,'DS/DY/DPT2',//)
 2    FORMAT(3X,'E DS/D3P=',//)
 3    FORMAT(3X,'DS/DY/DPT=',//)
      END
  
