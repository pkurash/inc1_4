C-THIS CODE IS COMPOSED OF:
C----------- COMBINATION OF STRUCTURE FUNCTIONS 
C----------- Whereas one finds in hadlib.f:
C----------- MATRIX ELEMENTS INTEGRATED ON  PHASE SPACE
C----------- CHANGE OF FACTORISATION SCHEME
C*===================================================================
C
C                 STRUCTURE FUNCTION AT SCALE Q2
C
C==================================================================  */
      SUBROUTINE STRU(XA,XB,XC,IH1,IH2,Q2,QP2,GPPV1,GPPV2,GPPV3,GPPV4
     # ,GPPV5,GPPV6,GPPV7,GPPV8,GPPV9,GPPV10,GPPV11,GPPV12,GPPV13
     # ,GPPV14,GPPV15,GPPV16)
      IMPLICIT REAL*8(A-H,M-Z)
      COMMON/PION/IPION
      include "protar.f"
 
C  /* HADRON IH1 and IH2: =0 for PROTON
C                         =1 for ANTI-PROTON */
      IP=IH1
      IPB=IH2
      IPI=IPION      
      Q2SAVE=Q2
      QP2SAVE=QP2
      if(iproject.eq.1) then
      ITI=1
      elseif(iproject.eq.2) then 
      itI=2
      endif
      CALL FONSTRU(ITI,XA,IP,Q2,XUHA,XUBHA,XDHA,XDBHA,XSHA,
     # XCHA,XBHA,XTHA,XGPROA)      
      if(Q2SAVE.NE.Q2) then
        write (15,*) 'IN HADLIB/stru.f after fonfru Q2SAVE',Q2SAVE,
     #' differ from Q2',Q2
      endif
      if(itarget.eq.1) then
      ITI=1
      elseif(itarget.eq.2) then
      ITI=2
      endif
      CALL FONSTRU(ITI,XB,IPB,Q2,XUHB,XUBHB,XDHB,XDBHB,XSHB,
     # XCHB,XBHB,XTHB,XGPROB)
      CALL FONFRA(XC,IPI,QP2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP,XDBP,XDBBP
     # ,XDTP,XDTBP,XDGP)      
      if(QP2SAVE.NE.QP2) then
        write (*,*) 'IN HADLIB/stru.f after fonfra QP2SAVE',QP2SAVE,
     #' differ from QP2',QP2
      endif
      IF (XA.LE.1.D-5.OR.XB.LE.1.D-5) THEN
        WRITE (15,*) 'IN subroutine STRU: XA OR XB OUT OF THEIR RANGE'
	stop
      ELSE
      GPPU1=XUHA*(XDHB+XSHB+XCHB+XBHB)*XDUP
     &+XDHA*(XUHB+XSHB+XCHB+XBHB)*XDDP+
     &XSHA*(XUHB+XDHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUHB+XDHB+XSHB+XBHB)*XDCP+
     &XBHA*(XUHB+XDHB+XSHB+XCHB)*XDBP+
     &(XUHA*XDUP+XDHA*XDDP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBP)*XTHB+XTHA*
     &(XUHB+XDHB+XSHB+XCHB+XBHB)*XDTP
     &+XUBHA*(XDBHB+XSHB+XCHB+XBHB)*XDUBP
     &+XDBHA*(XUBHB+XSHB+XCHB+XBHB)*XDDBP+
     &XSHA*(XUBHB+XDBHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUBHB+XDBHB+XSHB+XBHB)*XDCP+
     &XBHA*(XUBHB+XDBHB+XSHB+XCHB)*XDBBP+
     &(XUBHA*XDUBP+XDBHA*XDDBP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBBP)*XTHB+XTHA*
     &(XUBHB+XDBHB+XSHB+XCHB+XBHB)*XDTBP
      GPPV1=GPPU1/XA/XB/XC**3
C*********************************************************************
      GPPU2=(XUHA*(XDHB+XSHB+XCHB+XBHB+XTHB)
     &+XDHA*(XSHB+XCHB+XBHB+XTHB)+
     &XSHA*(XCHB+XBHB+XTHB)+
     &XCHA*(XBHB+XTHB)+
     &XBHA*(XTHB)
     &+XUBHA*(XDBHB+XSHB+XCHB+XBHB+XTHB)
     &+XDBHA*(XSHB+XCHB+XBHB+XTHB)+
     &XSHA*(XCHB+XBHB+XTHB)+
     &XCHA*(XBHB+XTHB)+
     &XBHA*(XTHB))*XDGP
      GPPV2=GPPU2/XA/XB/XC**3
C*********************************************************************
      GPPU3=XUHA*(XDBHB+XSHB+XCHB+XBHB)*XDUP
     &+XDHA*(XUBHB+XSHB+XCHB+XBHB)*XDDP+
     &XSHA*(XUBHB+XDBHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUBHB+XDBHB+XSHB+XBHB)*XDCP+
     &XBHA*(XUBHB+XDBHB+XSHB+XCHB)*XDBP+
     &(XUHA*XDUP+XDHA*XDDP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBP)*XTHB+XTHA*
     &(XUBHB+XDBHB+XSHB+XCHB+XBHB)*XDTP
     &+XUBHA*(XDHB+XSHB+XCHB+XBHB)*XDUBP
     &+XDBHA*(XUHB+XSHB+XCHB+XBHB)*XDDBP+
     &XSHA*(XUHB+XDHB+XCHB+XBHB)*XDSP+
     &XCHA*(XUHB+XDHB+XSHB+XBHB)*XDCP+
     &XBHA*(XUHB+XDHB+XSHB+XCHB)*XDBBP+
     &(XUBHA*XDUBP+XDBHA*XDDBP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBBP)*XTHB+XTHA*
     &(XUHB+XDHB+XSHB+XCHB+XBHB)*XDTBP
      GPPV3=GPPU3/XA/XB/XC**3
C*********************************************************************
      GPPU4=(XUHA*(XDBHB+XSHB+XCHB+XBHB+XTHB)
     &+XDHA*(XUBHB+XSHB+XCHB+XBHB+XTHB)+
     &XSHA*(XUBHB+XDBHB+XCHB+XBHB+XTHB)+
     &XCHA*(XUBHB+XDBHB+XSHB+XBHB+XTHB)+
     &XBHA*(XUBHB+XDBHB+XSHB+XCHB+XTHB)+
     &XTHA*(XUBHB+XDBHB+XSHB+XCHB+XBHB))
     &*XDGP
      GPPV4=GPPU4/XA/XB/XC**3
C*********************************************************************
      GPPU5=(XDHA*XDBHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB+XTHA*XTHB)*XDUP+(XUHA
     &*XUBHB+XSHA*XSHB+XBHA*XBHB+XCHA*XCHB
     &+XTHA*XTHB)*XDDP+(XUHA*XUBHB+XDHA
     &*XDBHB+XCHA*XCHB+XBHA*XBHB+XTHA*XTHB)
     &*XDSP+(XUHA*XUBHB+XDHA*XDBHB+
     &XSHA*XSHB+XBHA*XBHB+XTHA*XTHB)*XDCP+
     &(XUHA*XUBHB+XDHA*XDBHB+XSHA*XSHB
     &+XCHA*XCHB+XTHA*XTHB)*XDBP+(XUHA*
     &XUBHB+XDHA*XDBHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB)*XDTP+
     &(XDBHA*XDHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB+XTHA*XTHB)*XDUBP+(XUBHA
     &*XUHB+XSHA*XSHB+XBHA*XBHB+XCHA*XCHB
     &+XTHA*XTHB)*XDDBP+(XUBHA*XUHB+XDBHA
     &*XDHB+XCHA*XCHB+XBHA*XBHB+XTHA*XTHB)
     &*XDSP+(XUBHA*XUHB+XDBHA*XDHB+
     &XSHA*XSHB+XBHA*XBHB+XTHA*XTHB)*XDCP+
     &(XUBHA*XUHB+XDBHA*XDHB+XSHA*XSHB
     &+XCHA*XCHB+XTHA*XTHB)*XDBBP+(XUBHA*
     &XUHB+XDBHA*XDHB+XSHA*XSHB+XCHA*
     &XCHB+XBHA*XBHB)*XDTBP
      GPPV5=GPPU5/XA/XB/XC**3
C*********************************************************************
      GPPU6=XUHA*XUHB*XDUP+XDHA*XDHB
     &*XDDP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
     &+XUBHA*XUBHB*XDUBP+XDBHA*XDBHB
     &*XDDBP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
     &+XBHA*XBHB*XDBP+XBHA*XBHB*XDBBP+
     &XTHA*XTHB*XDTP+XTHA*XTHB*XDTBP
      GPPV6=GPPU6/XA/XB/2.D0/XC**3
C*********************************************************************
      GPPU7=(XUHA*XUHB+XDHA*XDHB
     &+XSHA*XSHB+XCHA*XCHB+XBHA*XBHB*2.
     &+XUBHA*XUBHB+XDBHA*XDBHB
     &+XSHA*XSHB+XCHA*XCHB+2.*XTHA*XTHB)*XDGP
      GPPV7=GPPU7/XA/XB/2.D0/XC**3
C*********************************************************************
      GPPU8=((XDHA+XSHA+XCHA+XBHA+XTHA)*XDUP+
     &(XUHA+XSHA+XCHA+XBHA+XTHA)*XDDP+
     &(XUHA+XDHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUHA+XDHA+XSHA+XBHA+XTHA)*XDCP+
     &(XUHA+XDHA+XSHA+XCHA+XTHA)*XDBP+
     &(XUHA+XDHA+XSHA+XCHA+XBHA)*XDTP+
     &(XDBHA+XSHA+XCHA+XBHA+XTHA)*XDUP+
     &(XUBHA+XSHA+XCHA+XBHA+XTHA)*XDDP+
     &(XUBHA+XDBHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUBHA+XDBHA+XSHA+XBHA+XTHA)*XDCP+
     &(XUBHA+XDBHA+XSHA+XCHA+XTHA)*XDBP+
     &(XUBHA+XDBHA+XSHA+XCHA+XBHA)*XDTP
     &)*XGPROB
      GPPV8=GPPU8/XA/XB/XC**3
C*********************************************************************
      GPPU9=((XDHA+XSHA+XCHA+XBHA+XTHA)*XDUBP+
     &(XUHA+XSHA+XCHA+XBHA+XTHA)*XDDBP+
     &(XUHA+XDHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUHA+XDHA+XSHA+XBHA+XTHA)*XDCP+
     &(XUHA+XDHA+XSHA+XCHA+XTHA)*XDBBP+
     &(XUHA+XDHA+XSHA+XCHA+XBHA)*XDTBP+
     &(XDBHA+XSHA+XCHA+XBHA+XTHA)*XDUBP+
     &(XUBHA+XSHA+XCHA+XBHA+XTHA)*XDDBP+
     &(XUBHA+XDBHA+XCHA+XBHA+XTHA)*XDSP+
     &(XUBHA+XDBHA+XSHA+XBHA+XTHA)*XDCP+
     &(XUBHA+XDBHA+XSHA+XCHA+XTHA)*XDBBP+
     &(XUBHA+XDBHA+XSHA+XCHA+XBHA)*XDTBP
     &)*XGPROB
      GPPV9=GPPU9/XA/XB/XC**3
C*********************************************************************
      GPPU10=(XUHA*XDUBP+XDHA*XDDBP+XSHA*
     &XDSP+XCHA*XDCP+XBHA*XDBBP+XTHA*XDTBP+
     &XUBHA*XDUP+XDBHA*XDDP+XSHA*XDSP+
     &XCHA*XDCP+XBHA*XDBP+XTHA*XDTP
     &)*XGPROB
      GPPV10=GPPU10/XA/XB/XC**3
C*********************************************************************
      GPPU11=XUHA*XUBHB*XDUP+XDHA*XDBHB
     &*XDDP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
     &+XBHA*XBHB*XDBP+XTHA*XTHB*XDTP+
     &XUBHA*XUHB*XDUBP+XDBHA*XDHB
     &*XDDBP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
     &+XBHA*XBHB*XDBBP+XTHA*XTHB*XDTBP
      GPPV11=GPPU11/XA/XB/XC**3
C*********************************************************************
      GPPU12=(XUHA*XUBHB+XDHA*XDBHB+XSHA*
     &XSHB+XCHA*XCHB+XBHA*XBHB+XTHA*XTHB)*XDGP
      GPPV12=GPPU12/XA/XB/XC**3
C*********************************************************************
      GPPU13=(XUHA*XDUP+XUBHA*XDUBP+XDHA*
     &XDDP+XDBHA*XDDBP+2.*XSHA*XDSP+
     &2.*XCHA*XDCP+XBHA*(XDBP+XDBBP)+XTHA*(
     &XDTP+XDTBP))*XGPROB
      GPPV13=GPPU13/XA/XB/XC**3
C*********************************************************************
      GPPU14=(XUHA+XUBHA+XDHA+XDBHA+2.*XSHA
     &+2.*XCHA+2.*XBHA+2.*XTHA)*XGPROB*XDGP
      GPPV14=GPPU14/XA/XB/XC**3
C*********************************************************************
      GPPU15=XGPROA*XGPROB*XDGP/2.
      GPPV15=GPPU15/XA/XB/XC**3
C*********************************************************************
      GPPU16=XGPROA*XGPROB*(XDUP+XDUBP+XDDP+XDDBP
     &+XDSP*2.+XDCP*2.+XDBP+XDBBP+XDTP+XDTBP)/2.
      GPPV16=GPPU16/XA/XB/XC**3
      ENDIF
      RETURN
      END
