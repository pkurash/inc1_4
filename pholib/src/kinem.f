      SUBROUTINE KINEM(V,W)
      implicit real*8(a-h,o-z)
	include "vari.f"
c      COMMON/W50511/NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
      include "pdfmod.f"
c      COMMON/PDFMOD/NPTYPRO,NGROPRO,NSETPRO,NPTYTAR,NGROTAR,
c     *  NSETTAR
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
c forces interstru
c      it=1
c      write(65,*) 'in kinem v',v,'w',w
C
      ALFMIN=2.*PI/ALFA
      T1LAST=0.D0
C
      V1 = 1.D0-V
      V2 = 2.D0-V
      W1 = 1.D0-W
      VW = V * W
      X1 = 1.D0-V*W
      X2 = 1.D0-V+V*W
      LV = LOG(V)
      LV1= LOG(V1)
      LW = LOG(W)
      LX1= LOG( X1 )
      LX2= LOG( X2 )
      LR = LOG( X1/V1 )
      LW1 = 0.D0
      LDQ  =0.D0
      IF (W.NE.1.D0) THEN
           LW1= LOG(1.D0-W)
c            LDQ=-2.D0*LX2-LOG(1.D0-X2)
      ENDIF
      T1 = XX1 /VW
      T2 = XX2 /V1
      SHAT = PT2/(V1*V*W)
      LSQ  = LOG(SHAT/Q1S)
      LQC  = LOG(SHAT/Q2S)
      LLQ  = LOG(SHAT/0.04D0)
      IF(QDS.LT.1.D-9) QDS  = SHAT
      IF(QDS.GE.1.D-9) LDQ  = LOG(SHAT/QDS)
      LLQ= LLQ - LDQ
      LLQCH=MAX(LOG(SHAT/2.0D0)-LDQ,0.D0)
c
c     define fragmentation logs for full inclusive cross sections
c     for use with ianom=2
      DLQ=LLQ
      DLQCH=LLQCH
c
c      IF(IDY.EQ.1.AND.W.NE.1.D0) THEN
c          LLQ=LOG(SHAT/(DYMIN**2*(1.D0-X2)))+DYNOR2/DYNOR1
c          LLQCH=MIN(LLQCH,LLQ)
c      ENDIF
      IF(ISOL.EQ.1.AND.W.NE.1.D0) THEN
          YSTAR=0.5D0*LOG(T1/T2)
c          LLQ=LOG((1.D0-X2*SIDEL4)/((1.D0-X2)**2*SIDEL4))
c          LLQCH=LLQ
          LLQ=LOG( (2.D0*COSH(YSTAR)/R0)**2 / (1.D0-X2)**2 )
          LLQCH=MIN(LLQCH,LLQ)
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
c           write(65,*) 'in kinem before fonfra old call'
c          CALL  FONFRA(XD,IPI,QP2,XDGP,XDUP,XDDP,XDSP,XDCP,XDBP) 
c	   write(65,*) 'in kinem after fonfra old call'
c          write(65,*) 'in kinem before fonfra '   
           CALL  FONFRA(XD,IPI,QP2,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDCP
     #                 ,XDBP,XDBBP,XDTP,XDTBP,XDGP)
c	    write(65,*) 'in kinem after fonfra '
           X2FRAG(1)=XDUP*ALFMIN
           X2FRAG(2)=XDDP*ALFMIN
           X2FRAG(3)=XDSP*ALFMIN
           X2FRAG(7)=XDCP*ALFMIN
           X2G=XDGP*ALFMIN
           X2FRAG(4)=X2FRAG(1)
           X2FRAG(5)=X2FRAG(2)
           X2FRAG(6)=X2FRAG(3)
           X2FRAG(8)=X2FRAG(7)
      ENDIF
c     end of photon anomalous functions in case ianom=2
c
      ALFAS=AA*PI
      ALFPI=AA*0.5D0
c      write(65,*) 'in kinem before interstru'
c      if(it.eq.1) then
      call interstru
c      write(65,*) 'in kinem after interstru'
c      WRITE(65,*) ' X t1=',T1,' Q2=',Q1S,' XG1=',GG1,' XU=',QQ1(1),
c     * ' XD=',QQ1(2)       
c     * ' XD=',QQ1(2),' XSEA=',QQ1(3)     
c      WRITE(65,*) ' X t2=',T2,' Q2=',Q1S,' XG1=',GG2,' XU=',QQ2(1),
c     * ' XD=',QQ2(2)     
c     * ' XD=',QQ2(2),' XSEA=',QQ2(3)
c      return
c      endif
      RETURN
      END
