      DOUBLE PRECISION FUNCTION CDEL(XX)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/HASCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      COMMON/HADR/IH1,IH2
      COMMON/FCTD/GRRT,GRRC
      COMMON/YRANGE/YMIN,YMAX
      COMMON/PTVALUE/PT
      COMMON/YVALUE/RAPIDITE
      common/INTCHOIX/integy
      include "protar.f"
      logical integy
      DIMENSION XX(10)
c      data ifirst/0/
      save 
      if (integy) then
        RAP = YMIN + (YMAX-YMIN)*XX(1)
      else
        RAP = RAPIDITE
      endif
      RS = DSQRT(GS)
      GV=1.D0-PT/RS*DEXP(-RAP)
      GW=PT*DEXP(RAP)/(RS-PT*DEXP(-RAP))
c      if(ifirst.eq.0) then
c      write(65,*)' IN CDEL RS ',RS,' PT ',PT,' RAP ',RAP,
c     1 ' YMIN ',YMIN, 'YMAX',YMAX
c      ifirst=1
c      endif
      X3MIN=1.D0-GV+GV*GW
      X3MAX=1.D0
      if (integy) then
        X3=X3MIN+(X3MAX-X3MIN)*XX(2)
      else
        X3=X3MIN+(X3MAX-X3MIN)*XX(1)
      endif
      VMAX=1.D0-(1.D0-GV)/X3
      VMIN=GV*GW/X3
      if (integy) then
        V = VMIN + (VMAX-VMIN)*XX(3)
      else
        V = VMIN + (VMAX-VMIN)*XX(2)
      endif
      A=GV*GW/V/X3
      UN=1.D0
      GHD=0.D0
      IF(A.GE.UN) THEN
        GOTO 11
      ENDIF
c      write(65,*)' IN CDEL GV',GV,' GW ',GW,' V ',V,' X3 ',X3
      BX1=GV*GW/V/X3
c      write(65,*)' IN CDEL BX1 from above is ',BX1
      BX2=(1.D0-GV)/(1.D0-V)/X3
      MU2=MU**2
      M2=M**2
      MP2=MP**2
c      write(65,*)' IN CDEL CALL STRU BX1 ',BX1, ' BX2 ',BX2
	iproject=1
	itarget=2
      CALL STRU(BX1,BX2,X3,IH1,IH2,M2,MP2,GRRT1,GRRT2,GRRT3,GRRT4
     #   ,GRRT5,GRRT6,GRRT7,GRRT8,GRRT9,GRRT10,GRRT11,GRRT12,GRRT13
     #   ,GRRT14,GRRT15,GRRT16)
c      write(65,*)' IN CDEL CALL STRU BX2 ',BX2, ' BX1 ',BX1
	iproject=2
	itarget=1
      CALL STRU(BX2,BX1,X3,IH2,IH1,M2,MP2,GRRC1,GRRC2,GRRC3,GRRC4
     #   ,GRRC5,GRRC6,GRRC7,GRRC8,GRRC9,GRRC10,GRRC11,GRRC12,GRRC13
     #   ,GRRC14,GRRC15,GRRC16)
      DO J0=1,16
        IF (J0.EQ.16.OR.J0.EQ.15) THEN
          CC=VC**2
        ELSE IF (J0.EQ.14.OR.J0.EQ.13.OR.J0.EQ.10.OR.J0.EQ.9.OR.J0
     #     .EQ.8) THEN
          CC=VC*N
        ELSE
          CC=N**2
        ENDIF
        IF (J0.EQ.16) THEN
          GRRT=GRRT16
          GRRC=GRRC16
        ELSE IF (J0.EQ.15) THEN
          GRRT=GRRT15
          GRRC=GRRC15
        ELSE IF (J0.EQ.14) THEN
          GRRT=GRRT14
          GRRC=GRRC14
        ELSE IF (J0.EQ.13) THEN
          GRRT=GRRT13
          GRRC=GRRC13
        ELSE IF (J0.EQ.12) THEN
          GRRT=GRRT12
          GRRC=GRRC12
        ELSE IF (J0.EQ.11) THEN
          GRRT=GRRT11
          GRRC=GRRC11
        ELSE IF (J0.EQ.10) THEN
          GRRT=GRRT10
          GRRC=GRRC10
        ELSE IF (J0.EQ.9) THEN
          GRRT=GRRT9
          GRRC=GRRC9
        ELSE IF (J0.EQ.8) THEN
          GRRT=GRRT8
          GRRC=GRRC8
        ELSE IF (J0.EQ.7) THEN
          GRRT=GRRT7
          GRRC=GRRC7
        ELSE IF (J0.EQ.6) THEN
          GRRT=GRRT6
          GRRC=GRRC6
        ELSE IF (J0.EQ.5) THEN
          GRRT=GRRT5
          GRRC=GRRC5
        ELSE IF (J0.EQ.4) THEN
          GRRT=GRRT4
          GRRC=GRRC4
        ELSE IF (J0.EQ.3) THEN
          GRRT=GRRT3
          GRRC=GRRC3
        ELSE IF (J0.EQ.2) THEN
          GRRT=GRRT2
          GRRC=GRRC2
        ELSE IF (J0.EQ.1) THEN
          GRRT=GRRT1
          GRRC=GRRC1
        ENDIF
        IF (ICHOI.EQ.1) THEN
          GHD=F0(V,X3)+GHD
        ELSE IF (ICHOI.EQ.2) THEN
          GHD=(FDEL(V,X3)+FVWPL(UN,V,X3)*DLOG(1.-A)+FVLO(UN,V,X3)*
     #      (DLOG(1.-A))**2/2.D0)/(8.*CC*(1.-V))*ALFAS(IFLAG,MU2)+
     #      F0(V,X3) + GHD
        ENDIF
      ENDDO
 11   CONTINUE
      CDEL=GHD*ALFAS(IFLAG,MU2)**2
      if (integy) then
        CDEL=CDEL*(X3MAX-X3MIN)*(VMAX-VMIN)*(YMAX-YMIN)
      else
        CDEL=CDEL*(X3MAX-X3MIN)*(VMAX-VMIN)
      endif
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION DPLUS(XX)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/HASCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/EDF/J0
      COMMON/HADR/IH1,IH2
      COMMON/FCTC/GPPT,GPPC
      COMMON/FCTD/GRRT,GRRC
      COMMON/YRANGE/YMIN,YMAX
      COMMON/PTVALUE/PT
      COMMON/YVALUE/RAPIDITE
      common/INTCHOIX/integy      
      include "protar.f"
      logical integy
      DIMENSION XX(10)
c      data ifirst/0/
      save
      if (integy) then
        RAP = YMIN + (YMAX-YMIN)*XX(1)
      else
        RAP = RAPIDITE
      endif
c      if(ifirst.eq.0) then
c      write(65,*)' IN DPLUS RAP ',RAP,' YMAX ',YMAX,' YMIN ',YMIN
c      ifirst=1
c      endif
      RS = DSQRT(GS)
      GV=1.D0-PT/RS*DEXP(-RAP)
      GW=PT*DEXP(RAP)/(RS-PT*DEXP(-RAP))
      X3MIN=1.D0-GV+GV*GW
      X3MAX=1.D0
      if (integy) then
        X3=X3MIN+(X3MAX-X3MIN)*XX(2)
      else
        X3=X3MIN+(X3MAX-X3MIN)*XX(1)
      endif
      VMIN=GV*GW/X3
      VMAX=1.D0-(1.D0-GV)/X3
      if (integy) then
        V = VMIN + (VMAX-VMIN)*XX(3)
      else
        V = VMIN + (VMAX-VMIN)*XX(2)
      endif
      WMIN=GV*GW/X3/V
      WMAX=1.D0
      if (integy) then
        W=WMIN + (WMAX-WMIN)*XX(4)
      else
       W=WMIN + (WMAX-WMIN)*XX(3)
      endif
      UN=1.D0
      X1=GV*GW/V/W/X3
      X2=(1.D0-GV)/(1.D0-V)/X3
      MU2=MU**2
      M2=M**2
      MP2=MP**2
c      write(65,*)' IN DPLUS CALL STRU X1 ',X1, ' X2 ',X2	
	iproject=1
	itarget=2
      CALL STRU(X1,X2,X3,IH1,IH2,M2,MP2,GPPT1,GPPT2,GPPT3,GPPT4
     #   ,GPPT5,GPPT6,GPPT7,GPPT8,GPPT9,GPPT10,GPPT11,GPPT12,GPPT13
     #   ,GPPT14,GPPT15,GPPT16)
c      write(65,*)' IN DPLUS CALL STRU X2 ',X2, ' X1 ',X1
	iproject=2
	itarget=1
      CALL STRU(X2,X1,X3,IH2,IH1,M2,MP2,GPPC1,GPPC2,GPPC3,GPPC4
     #   ,GPPC5,GPPC6,GPPC7,GPPC8,GPPC9,GPPC10,GPPC11,GPPC12,GPPC13
     #   ,GPPC14,GPPC15,GPPC16)
      BX1=GV*GW/V/X3
      BX2=(1.D0-GV)/(1.D0-V)/X3
c      write(65,*)' IN DPLUS CALL STRU BX1 ',BX1, ' BX2 ',BX2
	iproject=1
	itarget=2
      CALL STRU(BX1,BX2,X3,IH1,IH2,M2,MP2,GRRT1,GRRT2,GRRT3,GRRT4
     #   ,GRRT5,GRRT6,GRRT7,GRRT8,GRRT9,GRRT10,GRRT11,GRRT12,GRRT13
     #   ,GRRT14,GRRT15,GRRT16)
c      write(65,*)' IN CDPLUS CALL STRU BX2 ',BX2, ' BX1 ',BX1
	iproject=2
	itarget=1
      CALL STRU(BX2,BX1,X3,IH2,IH1,M2,MP2,GRRC1,GRRC2,GRRC3,GRRC4
     #   ,GRRC5,GRRC6,GRRC7,GRRC8,GRRC9,GRRC10,GRRC11,GRRC12,GRRC13
     #   ,GRRC14,GRRC15,GRRC16)
      GHE=0.D0
      DO J0=1,16
        IF (J0.EQ.16.OR.J0.EQ.15) THEN
          CC=VC**2
        ELSE IF (J0.EQ.14.OR.J0.EQ.13.OR.J0.EQ.10.OR.J0.EQ.9.OR.J0
     #     .EQ.8) THEN
          CC=VC*N
        ELSE
          CC=N**2
        ENDIF
        IF (J0.EQ.16) THEN
          GPPT=GPPT16
          GPPC=GPPC16
          GRRT=GRRT16
          GRRC=GRRC16
        ELSE IF (J0.EQ.15) THEN
          GPPT=GPPT15
          GPPC=GPPC15
          GRRT=GRRT15
          GRRC=GRRC15
        ELSE IF (J0.EQ.14) THEN
          GPPT=GPPT14
          GPPC=GPPC14
          GRRT=GRRT14
          GRRC=GRRC14
        ELSE IF (J0.EQ.13) THEN
          GPPT=GPPT13
          GPPC=GPPC13
          GRRT=GRRT13
          GRRC=GRRC13
        ELSE IF (J0.EQ.12) THEN
          GPPT=GPPT12
          GPPC=GPPC12
          GRRT=GRRT12
          GRRC=GRRC12
        ELSE IF (J0.EQ.11) THEN
          GPPT=GPPT11
          GPPC=GPPC11
          GRRT=GRRT11
          GRRC=GRRC11
        ELSE IF (J0.EQ.10) THEN
          GPPT=GPPT10
          GPPC=GPPC10
          GRRT=GRRT10
          GRRC=GRRC10
        ELSE IF (J0.EQ.9) THEN
          GPPT=GPPT9
          GPPC=GPPC9
          GRRT=GRRT9
          GRRC=GRRC9
        ELSE IF (J0.EQ.8) THEN
          GPPT=GPPT8
          GPPC=GPPC8
          GRRT=GRRT8
          GRRC=GRRC8
        ELSE IF (J0.EQ.7) THEN
          GPPT=GPPT7
          GPPC=GPPC7
          GRRT=GRRT7
          GRRC=GRRC7
        ELSE IF (J0.EQ.6) THEN
          GPPT=GPPT6
          GPPC=GPPC6
          GRRT=GRRT6
          GRRC=GRRC6
        ELSE IF (J0.EQ.5) THEN
          GPPT=GPPT5
          GPPC=GPPC5
          GRRT=GRRT5
          GRRC=GRRC5
        ELSE IF (J0.EQ.4) THEN
          GPPT=GPPT4
          GPPC=GPPC4
          GRRT=GRRT4
          GRRC=GRRC4
        ELSE IF (J0.EQ.3) THEN
          GPPT=GPPT3
          GPPC=GPPC3
          GRRT=GRRT3
          GRRC=GRRC3
        ELSE IF (J0.EQ.2) THEN
          GPPT=GPPT2
          GPPC=GPPC2
          GRRT=GRRT2
          GRRC=GRRC2
        ELSE IF (J0.EQ.1) THEN
          GPPT=GPPT1
          GPPC=GPPC1
          GRRT=GRRT1
          GRRC=GRRC1
        ENDIF
        CPL=((FVWPL(W,V,X3)/W-FVWPL(UN,V,X3))/(UN-W)+(
     #    FVLO(W,V,X3)/W- FVLO(UN,V,X3))*DLOG(Un-W)/(Un-W))/
     #    (8.D0*CC*(Un-V) )
        CCO=FRESC(W,V,X3)/W/(8.D0*CC*(Un-V))
        GHE=CPL+CCO + GHE
      ENDDO
      if (integy) then
        DPLUS=GHE*ALFAS(IFLAG,MU2)**3*(X3MAX-X3MIN)*(VMAX-VMIN)*
     #  (WMAX-WMIN)*(YMAX-YMIN)
      else
        DPLUS=GHE*ALFAS(IFLAG,MU2)**3*(X3MAX-X3MIN)*(VMAX-VMIN)*
     #  (WMAX-WMIN)
      endif
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FDEL(V,X3)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/HASCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/FCTD/GRRT,GRRC
      BX1=GV*GW/V/X3
      BX2=(1.D0-GV)/(1.D0-V)/X3
      SHD=BX1*BX2*GS
      MU2=MU**2
      UN=1.D0
      FKEL=GRRT*AVDEL(V,SHD)
      FKELC=GRRC*(AVDEL(Un-V,SHD)+DLOG(V/(UN-V))*AVWPL(UN,UN-V,SHD)
     # +.5D0*AVLO(UN,UN-V,SHD)*(DLOG((UN-V)/V))**2)*(UN-V)/V
      FDEL=(FKEL+FKELC)/SHD
      RETURN
      END
      DOUBLE PRECISION FUNCTION FVWPL(W,V,X3)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/HASCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/FCTC/GPPT,GPPC
      COMMON/FCTD/GRRT,GRRC
      X1=GV*GW/V/W/X3
      X2=(1.D0-GV)/(1.D0-V)/X3
      VX=1.D0-V*W
      WX=(1.D0-V)/(1.D0-V*W)
      SH=X1*X2*GS
      IF (W.EQ.1.D0) THEN
        FPPT=GRRT
        FPPC=GRRC
      ELSE
        FPPT=GPPT
        FPPC=GPPC
      ENDIF
      RVWPL=AVWPL(W,V,SH)*FPPT
      RVWPLC=(AVWPL(WX,VX,SH)+AVLO(WX,VX,SH)*DLOG(V/VX))*VX/V*FPPC
      FVWPL=(RVWPL+RVWPLC)/SH
      RETURN
      END
      DOUBLE PRECISION FUNCTION FVLO(W,V,X3)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/HASCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/FCTC/GPPT,GPPC
      COMMON/FCTD/GRRT,GRRC
      X1=GV*GW/V/W/X3
      X2=(1.D0-GV)/(1.D0-V)/X3
      VX=1.D0-V*W
      WX=(1.D0-V)/(1.D0-V*W)
      SH=X1*X2*GS
      IF (W.EQ.1.D0) THEN
        FPPT=GRRT
        FPPC=GRRC
      ELSE
        FPPT=GPPT
        FPPC=GPPC
      ENDIF
      RVWLO=AVLO(W,V,SH)*FPPT
      RVWLOC=AVLO(WX,VX,SH)*FPPC*VX/V
      FVLO=(RVWLO+RVWLOC)/SH
      RETURN
      END
      DOUBLE PRECISION FUNCTION FRESC(W,V,X3)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/HASCALE/M,MP,MU
      COMMON/ORDE/IFLAG,ICHOI,AL,CQ,V1,V2,V3,V4
      COMMON/FCTC/GPPT,GPPC
      X1=GV*GW/V/W/X3
      X2=(1.D0-GV)/(1.D0-V)/X3
      VX=1.D0-V*W
      WX=(1.D0-V)/(1.D0-V*W)
      SH=X1*X2*GS
      RRESC=(STRUV(W,V,X3,SH)+AVGO(W,V))*GPPT
      RRESCC=(STRUV(WX,VX,X3,SH)+AVGO(WX,VX))*GPPC
      FRESC=(RRESC+RRESCC)/SH
      RETURN
      END
C    THIS PROGRAM ALLOWS TO MOVE FROM A FACTORIZATION SCHEME TO ANOTHER
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C    FRAGMENTATION   US(CQ=0) - MSBAR
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION FQQD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      PI=4.*DATAN(1.D0)
      PI2=PI**2
      FQQD=4.D0/3.D0*(-9.D0/2.D0+2.D0*PI2/3.D0)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FQQW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FQQW=4.D0/3.D0*(-3.D0/2.D0+2.D0*(1.D0+X**2)*
     1      LOG(X)+(1.D0-X)**2*3.D0/2.D0)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FQQL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FQQL=4.D0/3.D0*(1.D0+X**2)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FQGL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FQGL=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FQGD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FQGD=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FQGW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FQGW=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FGGL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FGGL=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FGGD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FGGD=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FGGW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FGGW=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FGQL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FGQL=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FGQD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FGQD=0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FGQW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      FGQW=0.D0
      RETURN
      END
C    THIS PROGRAM ALLOWS TO MOVE FROM A FACTORIZATION SCHEME TO ANOTHER
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  JMAR=0 OU 2    DIFFERENCE BETWEEN US(CQ=0) -  MSBAR SCHEME
C  JMAR=3         DIFFERENCE BETWEEN US(CQ=0) -  MARTINELLI
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION CGQD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/VALU/JMAR,IPT
      PI2=PI**2
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CGQD=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CGQD=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CGQD=-4.D0/3.D0*(9.D0/2.D0+PI2/3.D0)
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CGQW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CGQW=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CGQW=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CGQW=4.D0/3.D0*(-3.D0/2.D0-(1.D0+X**2)*
     1       DLOG(X)+(1.D0-X)*(3.D0+2.D0*X))
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CGQL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CGQL=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CGQL=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CGQL=4.D0/3.D0*(1.D0+X**2)
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CQQD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/VALU/JMAR,IPT
      PI2=PI**2
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CQQD=-(9.D0/2.D0+PI2/3.D0)*CF
      ELSE IF (JMAR.EQ.1) THEN
        CQQD=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CQQD=0.D0
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CQQW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CQQW=CF*(-3.D0/2.D0+(3.D0+2.D0*X)*
     1       (1.D0-X)-(1.D0+X**2)*DLOG(X))
      ELSE IF (JMAR.EQ.1) THEN
        CQQW=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CQQW=0.D0
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CQQL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CQQL=4.D0/3.D0*(1.D0+X**2)
      ELSE IF (JMAR.EQ.1) THEN
        CQQL=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CQQL=0.D0
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CQGD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CQGD=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CQGD=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CQGD=0.D0
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CQGW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CQGW=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CQGW=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CQGW=-(1.D0-X)/2.D0*(-(X**2+(1.D0-X)**2)
     1       *DLOG(X)+8.D0*X*(1.D0-X)-1.D0)
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CQGL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CQGL=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CQGL=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CQGL=-(1.D0-X)/2.D0*(X**2+(1.D0-X)**2)
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CGGD(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/VALU/JMAR,IPT
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CGGD=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CGGD=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CGGD=0.D0
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CGGW(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/VALU/JMAR,IPT
      NF=GTR*2.D0
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CGGW=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CGGW=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CGGW=2.D0*NF*(1.D0-X)/2.D0*(-(X**2+(1.D0-X)**2)
     1       *DLOG(X)+8.D0*X*(1.D0-X)-1.D0)
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION CGGL(X)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/CONS/PI,GS,GV,GW,N,GTR,CF,PT2,VC
      COMMON/VALU/JMAR,IPT
      NF=GTR*2.D0
      IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
        CGGL=0.D0
      ELSE IF (JMAR.EQ.1) THEN
        CGGL=0.D0
      ELSE IF (JMAR.EQ.3) THEN
        CGGL=2.D0*NF*(1.D0-X)/2.D0*(X**2+(1.D0-X)**2)
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION SNF(SCALE)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/ALFA/MASCH2,MASBO2,MASTO2,LAMBDA2
      IF(SCALE.GT.MASTO2) THEN
        SNF=6.D0
      ELSE IF(SCALE.GT.MASBO2) THEN
        SNF=5.D0
      ELSE IF(SCALE.GT.MASCH2) THEN
        SNF=4.D0
      ENDIF
      RETURN
      END
      SUBROUTINE FICT
      IMPLICIT REAL*8 (A-Z)
      COMMON/FACT/ZQ
      ZQ=0.D0
      END
