      SUBROUTINE KINEPS(V,W)
      implicit real*8(a-h,o-z)
 	include "vari.f"
c	COMMON/W50511/NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
	include "pdfmod.f"
c      COMMON/PDFMOD/NPTYPRO,NGROPRO,NSETPRO,NPTYTAR,NGROTAR,
c     * NSETTAR
	include "facto.f"
	include "qcd.f"
         include "phoscale.f"
        include "kine.f"
 	include "dreyan.f"
	include "projec.f"
      COMMON/COUPLI/AA
      COMMON/PARTN1/QQ(7)
      COMMON/PARTN2/QQQ(7)
      COMMON/ANOMAL/IANOM,IANORD
      COMMON/X2ANOM/X2FRAG(8),X2G,DLQ,DLQCH
      save
C	forces interstru
c      it=1
c
      ALFMIN=2.*PI/ALFA
c
      V1 = 1.D0-V
      V2 = 2.D0-V
      W1 = 1.D0-W
      VW = V * W
      X1 = 1.D0-V*W
      X2 = 1.D0-V+V*W
      T1 = XX1 /VW
      T2 = XX2 /V1
      SHAT = PT2/(V1*V*W)
      LDQ  =0.D0
      LLQ= LOG(SHAT/0.04D0)
      IF (W.NE.1.D0) THEN
c           LDQ=-2.D0*LX2-LOG(1.-X2)
      ENDIF
      IF(QDS.LT.1.D-9)  QDS  = SHAT
      IF(QDS.GE.1.D-9)  LDQ  = LOG(SHAT/QDS)
          LLQ= LLQ - LDQ
          LLQCH=MAX(LOG(SHAT/2.0D0)-LDQ,0.D0)
c
c     define fragmentation logs for full inclusive cross sections
      DLQ=LLQ
      DLQCH=LLQCH
c
      IF(ISOL.EQ.0) THEN
          IF(IDY.EQ.1.AND.W.NE.1.D0) THEN
              LLQ=LOG(SHAT/(DYMIN**2*(1.D0-X2)))+DYNOR2/DYNOR1
              LLQ= LLQ - LDQ
              LLQCH=MIN(LLQCH,LLQ)
          ENDIF
      ENDIF
c   define the log. factor associated to the accompanied photon
      IF(ISOL.EQ.1.AND.W.NE.1.D0) THEN
          YSTAR=0.5D0*LOG(T1/T2)
          FAC=SHAT*(X2*(1.D0-X2)/(2.D0*COSH(YSTAR)/R0))**2
          LLQ=LOG(FAC/.04D0)
          LLQCH=MAX(LOG(FAC/2.0D0),0.D0)
      ENDIF
c
c     call photon anomalous functions in case ianom=2
      IF(IANOM.EQ.2) THEN
           QP2=QDS
           XD=X2
           IPI=0
c
c  new fragmentation routines corrected by luc in february 1997
c
c           CALL  FONFRA(XD,IPI,QP2,XDGP,XDUP,XDDP,XDSP,XDCP,XDBP) 
c             
           CALL  FONFRA(XD,IPI,QP2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     #                 ,XDBP,XDBBP,XDTP,XDTBP,XDGP)
           X2FRAG(1)=XDUP*ALFMIN
           X2FRAG(2)=XDDP*ALFMIN
           X2FRAG(3)=XDSP*ALFMIN
           X2FRAG(7)=XDCP*ALFMIN
           X2G=XDGP
           X2FRAG(4)=X2FRAG(1)
           X2FRAG(5)=X2FRAG(2)
           X2FRAG(6)=X2FRAG(3)
           X2FRAG(8)=X2FRAG(7)
      ENDIF
c     end of photon anomalous functions in case ianom=2
c
      ALFAS=AA*PI
      ALFPI=AA*0.5D0     
c      if(it.eq.1) then
      call interstru
      return
c	endif
      END
