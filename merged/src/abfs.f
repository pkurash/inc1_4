CDECK  ID>, ABFS.
      SUBROUTINE ABFS(KHIOR,totint,anoint,borint,higint)
c calculates cross section with data input via common
c  KHIOR selects order
c  if =0 only born; if=1 born and anoma; if=2 born, anoma and NLL
c total: totint, anoma:anoint, born:borint; high order:higint    
      implicit real*8(a-h,o-z)
      DIMENSION YGHT(20)
C     DIMENSION AAA(15)
c*** commons to handle both programs
	common/output_units/io,ier,itest
	common/choics/facin
c  facin choice of cross section
       common/scale_choice/cc1,cc2,ccd
c   q1s=cc1*pt2  factorization scale
c   q3s=cc2*pt2	 renormalization scale
c   qds=ccd*pt2  fragmentation scale
 	common/constants/cn
C  CN cross section unit
C ISIGM selects DS/DYDPT2 (1) Ed3s/dp3 (2) DS/DYDPT (3)
 	COMMON/XSECT/ISIGM,ILOOP,IHIOR
C***
	
      COMMON/INTEGPREC/NVINT,NWINT,NY2
      COMMON/DATAEXP/IPROJ(20),SIGEX(20),ERREX(20),CORAP(20)
     *,ERSTA(20),ERSYS(20)
      include "conv.f"
      COMMON/BLANKO/BE1(20,20),BE2(20,20),BE12(20,20),BJ1(20,20)
     1  ,BJ2(20,20),BK12(20,20),BOR1(20,20),BOR2(20,20),BORLL(20,20)
     1  ,BCH1(20,20),BCH2(20,20),BCH3(20,20),BCHLL(20,20)
     1  ,GLUCH(20,20),CHCHB(20,20),BCHA(20,20)
     1  ,GLUQU(20,20),QUQUB(20,20),BORN(20,20),BOR3(20,20)
     1  ,BC1(20,20),BC2(20,20),BC3(20,20),BC12(20,20),BD12(20,20)
     1  ,BA1(20,20),BA2(20,20),BA3(20,20),BA12(20,20),BB1(20,20)
     1  ,BB2(20,20),BB3(20,20),BB12(20,20),BJ12(20,20)
     1  ,BA(20,20),BB(20,20),BC(20,20),BD(20,20),BE(20,20),BJ(20,20)
     1  ,BK(20,20)
     1,VNT(40),VGHT(40),WNT(40),WGHT(40)
     1 ,BS(20,20),BAS(20,20),BBS(20,20),BCS(20,20),BES(20,20),BJS(20,20)
      include "facto.f"
      include "qcd.f"
      include "phoscale.f"
c      COMMON/PHOSCALE/Q1S,Q2S,QDS,Q0S,ALBD,ALBD2,XLG0
	include "kine.f"
c     COMMON/KINE/S,RS,PT,PT2,Y,VV,WW,XX1,XX2,XX1M
      COMMON/CHARGE/Q(8)
	include "projec.f"
	include "dreyan.f"
      include "opt.f"
	include "vari.f"
      COMMON/DBOR/DBOR1,DBOR2,DBOR3,DBCH1,DBCH2,DBCH3,DBR1,DBR2,DBR3
     *,DA2,DA3,DA12,DB2,DB3,DB12,DBC1,DBC2,DBC3
     *,DC2,DC3,DC12,DD12,DE2,DE12,DJ2,DJ12,DK12
      COMMON/COUPLI/AA
      COMMON/OPTIM/AB,AC,BBORN,ANOMA,HIGOR
      DIMENSION C1(15),CD(15)
      EXTERNAL STEVE,OPTI,ALI
      EXTERNAL FGLUQU,FQUQUB
C ***      DATA C1/.08,.1,.25,1.,4.,.4,.5,.75,1.,2.,5.,10.,3*20./
C ***      DATA C1/.06,.15,.25,.35,.4,.5,.6,.8,1.,2.,4.,10.,20.,50.,100./
C ***      DATA C1/.11111111,.25,.35,.5,.6,.75,1.,1.5,2.,4.,5*4./
      DATA C1/.11111111D0,.15D0,.20D0,.25D0
     * ,.35D0,.4D0,.5D0,.6D0,.7D0,.8D0,.9D0
     * ,1.D0,1.2D0,1.5D0
     *,2.D0/
C     DATA C1/.0625,.1,.35,1.,5.,10.,50.,8*100./
C      DATA CD/.04,.10,.30,1.,3.0,10.,25.,8*100./
      DATA CD/.11111111D0,.25D0,.35D0,
     * .5D0,.75D0,1.D0,1.5D0,2.D0,4.D0,6*4.D0/
C     DATA C1/.0625,.25,1.,16.,.35,.4,.5,.6,.25,1.,4.,4.,10.,16.,100./
C     DATA C1/100.,50.,20.,10.,5.,3.,2.,1.,.8,.6,.5,.4,.3,.2,.1/
C *** DATA C1/20.,50.,50.,100.,100.,20.,5.,.6,.1,6*.08/
c      DATA AAA/.06,.07,.08,.09,.10,.11,.12,.13,.14,.15,.16,.17,.18
c     1        ,.19,.20/
c      DATA AAA/.05,.1,.15,.20,.25,.30,.161,.151,.140
c     1        ,.128,.117,4*.0293/
c      DATA AAA/.12,.14,.16,.18,.22,.30,.161,.151,.140
c     1        ,.128,.117,4*.0293/
c
	data ifirst/0/
c	save ifirst
	save
c
c        number of integration points (set in INTEGPREC)
c
	if(ifirst.eq.0) then
      write(io,*) ' ABFS NUMBER OF INTEGRATION POINTS : NVINT=',NVINT
     1,'   NWINT=',NWINT
c      PRINT 101
c101   FORMAT(1X,13D10.3)
c     number of rapidity integration points
	if(integy.eq.1) then
      write(io,*) 'NUMBER OF INTEGRATION POINTS in Y: NY2=',NY2
c        else
c      write(io,*) ' NO Y INTEGRATION '	
	endif
        endif 
c *** 
	i68=0
	i66=0
	i67=0
c *** definition des options  DO NOT CHANGE THEM
	ichi2=0
	iopt=0
	idy=0
	isol=0
	ieps=0
	ixbond=0
	icoupls=1
c       icoupls=10 	
	do i=1,20
	  corap(i)=1.D0
	enddo
	
c
c do loop on the coupling aa=alphas/pi values 

      DO 9999 IJKL=1,icoupls
c
c     select (factorization scale)**2=c1(mn)*pt**2(1-xt) for iopt=1
c     select (factorization scale)**2=c1(mn)*pt**2       for iopt=0
c
      if(ifirst.eq.0) then
        write(io,102)
        IF(IOPT.EQ.0) then
	write(io,*) ' RENOR SCALE**2 = ',cc2,'*PT2 (cc2 is used)'
	write(io,*) ' FACTO SCALE**2 = ',cc1,'*PT2 (cc1 is used)'
	endif
	endif
c      IF(IOPT.EQ.0) write(io,*) ' RENOR = FACTO SCALE**2 = ',cc1,'*PT2'
c      IF(IOPT.EQ.1) write(io,*) ' FACTO SCALE**2 = C1(MN)*PT2*(1-XT)'
c        
c possiblity for (1-xt) factorisation scale even with iopt=0
c      IF(IOPT.EQ.0) write(io,*) ' FACTO SCALE**2 = C1(MN)*PT2*(1-XT)'
c      IF(IOPT.EQ.0) write(io,*) ' RENOR SCALE**2 = C1(MN)*PT2*(1-XT)'
c
c possiblity for pt2 factorisation scale even with iopt=1
c      IF(IOPT.EQ.1) write(io,*) ' FACTO SCALE**2 = C1(MN)*PT2'
c
c      mn1 , mn2 are extreme indices on c1(mn) loop
      MN1=1
      MN2=15
      MCD1=1
C      MCD2=11
      MCD2=1
      IF(IOPT.EQ.0) THEN
                 MN1=1
                 MN2=1
      ENDIF
c
c   select factorization scale for final state fragmentation
c   which enters in unisolated photon case :  qds=ccd*pt**2
c   if ccd = 0 then the scale is shat
c
c
c  loop on the fragmentation scale
c
      DO 900 ICD=MCD1,MCD2
c      FRAGM2=EXP(.7D0+.5D0*(ICD-1))
c      CCD=CD(ICD)
c      CCD=0.E0
c      write(io,102)
c      write(io,102)
c      write(io,*) ' FACTO SCALE**2 IN FRAG =', FRAGM2
c      IF(CCD.LT.1.E-9) THEN
c           IF(ISOL.EQ.0)   write(io,*) ' FACTO SCALE**2 IN FRAG = SHAT'
c           ELSE
c      IF(ISOL.EQ.0) write(io,*) ' FACTO SCALE**2 IN FRAG = ',CCD,'*PT2'
c      ENDIF
c
      ALFA = 1.D0/137
      ALFA2=ALFA**2
c      PI   = 3.14159D0
       PI=4.D0*DATAN(1.D0)
      PI2  = PI**2
c -----     muon mass , proton mass
      AMU=.106D0
      AMP=.93828D0
c    dif. cross-section in PICO - barns/gev**2
c      CN   = 0.38935 D+9
      CF   = 4.D0/3.D0
      NC   = 3.D0
c
c     4 flavors  in alfas
c      NF=4
      B1=102.D0-38.D0*NF/3.D0
      B0=(11.D0*NC-2.D0*NF)/3.D0
      B0INV=4.D0*PI/B0
      AB=0.5D0*B0
      AC=.25D0*B1/B0
c
c     4 flavors in particle production
c      NF = 4
      if(ifirst.eq.0) then
      write(io,102)
      write(io,*) ' NUMBER OF OPEN FLAVORS=',NF
      write(io,102)
      endif
C
      Q(1) = 2.D0/3.D0
      Q(2) = -1.D0/3.D0
      Q(3) = -1.D0/3.D0
      Q(4) = -2.D0/3.D0
      Q(5) = 1.D0/3.D0
      Q(6) = 1.D0/3.D0
      Q(7)=  2.D0/3.D0
      Q(8)= -2.D0/3.D0
C     Q(7)=  0.D0
C     Q(8)= -0.D0
	if(ifirst.eq.0) then
      IF (Q(7).NE.0.D0) THEN
          write(io,*) ' CHARM INCLUDED IN STRUCTURE FUNCTIONS'
          ELSE
          write(io,*) ' CHARM NOT INCLUDED IN STRUCTURE FUNCTIONS'
      ENDIF
        endif
c
c         conventions     snoitnevnoc     conventions
c
      FQQ = 1.D0
      IF(IORD.EQ.1.AND.ICONV.EQ.0) FQQ=0.D0
      FQG = FQQ
      FGQ = 0.D0
C *** FQP = 0.
C *** FGG = 0.
      DPQ = 0.D0
C *** DQQ = 0.
C *** DGQ = 0.
C *** DQG = 0.
C *** DGG = 0.
C     KINEMATICS    SCITAMENIK    KINEMATICS
C
C
      S  = RS**2
      PLAB = (S - AMP**2) / (2.*AMP)
       if(ifirst.eq.0) then
c      write(io,102)
c      write(io,102)
      IF(IOPT.EQ.1) write(io,*) '    OPTIMIZATION OF SIGMA TOTAL'
c      write(io,102)
      write(io,*) '    square root of s=',RS,'GeV        FQQ=',FQQ
     1       ,'        FQG=',FQG,'        DPQ=',DPQ
c      write(io,102)
        endif
C
C START LOOP ON PT AND RAPIDITY
C
C    INITIATIALIZE THE CHI2
      CHI2=0.D0
C
      DO 100 I=IPT1,IPT2
c      write(io,102)
c      IPRO=IPROJ(I)
      PT=XPT(I)
      PT2 = PT**2
      XT=2.*PT/RS
c
c define fragmentation scale here
c
c	CCD=FRAGM2/PT2
         if(ifirst.eq.0) then
	  ifirst=1
      	IF(CCD.LT.1.D-9) THEN
           IF(ISOL.EQ.0)   write(io,*) ' FACTO SCALE**2 IN FRAG = SHAT'
           ELSE
      IF(ISOL.EQ.0) 
     1write(io,*) ' FACTO SCALE**2 IN FRAG = ',CCD,'*PT2 (ccd is used)'
      	ENDIF
	endif
c
c
C ----- NORMALIZATION FACTORS: ANORBO( BORN )   ,   ANORM( H.0. TERMS )
C                              SEE FURTHER CHANGES FOR LLQ IN KINEM
      ANORBO=1.
      ANORM = CN*ALFA/PT2**2

      IF(PT.GT.RS/2.) then
      write(ier,*)' PT ',pt,' is higher than sqrt(s)/2', rs/2.
      stop
c      GO TO 100
      endif
c
c	prepare for rapidity
c
      CHYMAX = 0.5D0*RS/PT
      SHYMAX = SQRT(CHYMAX**2-1.D0)
      YMAX   = LOG(CHYMAX+SHYMAX)
      write(io,*)' in abfs YMAX',ymax
      IF(IXBOND.EQ.1) THEN
          YINF=LOG(XINF/XT+SQRT((XINF/XT)**2+1.D0))
          YSUP=LOG(XSUP/XT+SQRT((XSUP/XT)**2+1.D0))
      ENDIF
      YSUP=DMIN1(YMAX,YSUP)
      YINF=DMAX1(-YMAX,YINF)
c      write(io,*) 'IN SUBR YINF ',YINF,' YSUP ',YSUP
      IF(INTEGY.EQ.1) THEN
          JY1=1
          JY2=NY2
          CALL DGSET(YINF,YSUP,NY2,XY,YGHT)
          TOTINT=0.D0
          BORINT=0.D0
          ANOINT=0.D0
          HIGINT=0.D0
      ENDIF
c         write(io,*)' before rap jy1',jy1,'jy2',jy2
c
c	rapidity loop
c
      DO 200 J=JY1,JY2
c      write(io,102)
      ISTOP=0
      TOTOLD=1.D+10
      TOTMIN=1.D+10
      TOTCHA=0.D0
      Y=XY(J)
      CHY=COSH(Y)
      write(io,*)' in abfs y is ',y
      IF ( ABS(Y).GT.YMAX ) then    
      write(ier,*)' in abfs y ',abs(y),' is greater than YMAX',ymax
      stop
c      GO TO 200
      endif
c
c	loop on factorisation scale
c
      DO 1 MN=MN1,MN2
C *** DO 1 MN=I,I
c      CC1=C1(MN)
c       Q1S=EXP(.7D0+.2*(MN-1))
c       CC1=Q1S/PT2
c      CC1=Q1S/PT2/(1.-2.*PT/RS)
c
C    CHOOSE THE FACTORIZATION SCALES
c possiblity for (1-xt) factorisation scale even with iopt=0
c      Q1S=CC1*PT2*(1.-2.*PT/RS)
      IF(IOPT.EQ.0) Q1S=CC1*PT2
      IF(Q1S.LT.2.) THEN
      write(itest,*)' in abfs Q1S',Q1S,' too small cc1=',cc1 
      Q1S=2.01D0
      if(iopt.eq.0) write(itest,*)'in abfs cc1 changed to ',q1s/pt2
      CC1=Q1S/(PT2*(1.-2.D0*PT/RS))
      IF(IOPT.EQ.0) CC1=Q1S/PT2
      ENDIF
c
c   choose the renormalization scale
      Q3S=cc2*pt2
c      write(io,*)' IN SUBR RENO SCALE Q3S ',Q3S

c
C RAZ GENERALE
      L=LOCF(BJS(20,20))-LOCF(BE1(1,1))+2      
      CALL VZERO(BE1,L)
C
C   INITIALIZE THE COUPLING CONSTANT AA=ALFAS/PI=1
      AA=1.
c
      Q2S=RLAM2
c
c	define fragmentation scale from ccd
c	see at the beginning of do 100 loop on p_t
c
      QDS=CCD*PT2
      IF(QDS.LT.2.D0) QDS=2.01D0
      IF(CCD.LT.1.D-9) QDS=0.D0
c      write(io,*)' IN SUBR FRAG SCALE QDS',QDS
c
c      write(io,102)
c      write(io,*) '   Q1S=',Q1S,' FRAG SCALE QDS=',QDS
c      write(io,102)
c
c
c      BOUNDS OF INTEGRATION
c
      XX1  =  PT*EXP(Y) / RS
      XX1M  =  1.D0-XX1
      XX2  =  PT*EXP(-Y)/RS
      VV   = 1.D0- XX2
      WW   = XX1 / VV
      
c                  DO V-INTEGRATION
c
      CALL DGSET(XX1,VV,NVINT,VNT,VGHT)
c      write(io,*) 'before nvint'
      DO 210 K=1,NVINT
C    INITIALIZE AND DO W-INTEGRATION
C
      L=LOCF(DK12)-LOCF(DBOR1)+2
c      CALL VZERO(DBOR1,L)
      CALL VZERO(BE1,L)
      V = VNT(K)
c      write(io,*) 'before kinem k',k,'v',v
      CALL KINEM(V,1.D0)
c      write(io,*) 'after kinem'
      GLUQU(I,J) =GLUQU(I,J)+VGHT(K)*FGLUQU(V)
      QUQUB(I,J) =QUQUB(I,J)+VGHT(K)*FQUQUB(V)
C
      GLUCH(I,J) =GLUCH(I,J)+VGHT(K)*FGLUCH(V)
      CHCHB(I,J) =CHCHB(I,J)+VGHT(K)*FCHCHB(V)
C
c      if(khior.ge.1) then
      BA10 = FBA1(V)
      BA20 = FBA2(V,1.D0)
      BA30 = FBA3(V,1.D0)
C
      BB10 = FBB1(V)
      BB20 = FBB2(V,1.D0)
      BB30 = FBB3(V,1.D0)
C
      BC10 = FBC1(V)
      BC20 = FBC2(V,1.D0)
      BC30 = FBC3(V,1.D0)
C
      BJ10 = FBJ1(V)
      BJ20 = FBJ2(V,1.D0)
C
      BE10 = FBE1(V)      
      BE20 = FBE2(V,1.D0)
c      BE20 = FBE2(V,1.D0)      
c      endif
      WMIN = XX1/VNT(K)
      ALMIN = LOG(1.D0-WMIN)
      
c
c ---- calculate the anomalous part with energy cut: phad/pgam < eps
c
      IF(IEPS.EQ.1) THEN
         WCUT=DMAX1( WMIN,1.D0-EPS/((1.D0+EPS)*V) )
         CALL DGSET(WCUT,1.D0,NWINT,WNT,WGHT)
         DO 215 L=1,NWINT
         W = WNT(L)
         CALL KINEPS(VNT(K),W)
         DBR1=DBR1 + WGHT(L)*FBOR1(VNT(K),WNT(L))
         DBR2=DBR2 + WGHT(L)*FBOR2(VNT(K),WNT(L))
         DBR3=DBR3 + WGHT(L)*FBOR3(VNT(K),WNT(L))
C
         DBC1=DBC1 + WGHT(L)*FBCH1(VNT(K),WNT(L))
         DBC2=DBC2 + WGHT(L)*FBCH2(VNT(K),WNT(L))
         DBC3=DBC3 + WGHT(L)*FBCH3(VNT(K),WNT(L))
C
215      CONTINUE
      ENDIF
c
c----calculate the Higher order and anomalous without eps cut
c
      CALL DGSET(WMIN,1.D0,NWINT,WNT,WGHT)
c      write(io,*)' before nwint'
      DO 220 L=1,NWINT
      W = WNT(L)
      ALWMIP = LOG(1.D0-W)/(1.D0-W)
c     write(io,*)' before kinem nwint',k
      CALL KINEM(VNT(K),WNT(L))
c      write(io,*)' after kinem nwint',k
c
c----if no angular cut the anomalous part here is not needed if eps=1
c
      IF(ISOL.EQ.0.AND.IEPS.EQ.1) GO TO 216
         DBOR1=DBOR1 + WGHT(L)*FBOR1(VNT(K),WNT(L))
         DBOR2=DBOR2 + WGHT(L)*FBOR2(VNT(K),WNT(L))
         DBOR3=DBOR3 + WGHT(L)*FBOR3(VNT(K),WNT(L))
C
         DBCH1=DBCH1 + WGHT(L)*FBCH1(VNT(K),WNT(L))
         DBCH2=DBCH2 + WGHT(L)*FBCH2(VNT(K),WNT(L))
         DBCH3=DBCH3 + WGHT(L)*FBCH3(VNT(K),WNT(L))
C
216   CONTINUE
	if(KHIOR.eq.2) then
      DA2 = DA2 + WGHT(L) * (FBA2(V,W)-BA20) / (1.-W)
      DA3 = DA3 + WGHT(L) * (FBA3(V,W)-BA30) * ALWMIP
      DA12=DA12 + WGHT(L) * FBA12(V,W)
C
      DB2 = DB2 + WGHT(L) * (FBB2(V,W)-BB20) / (1.-W)
      DB3 = DB3 + WGHT(L) * (FBB3(V,W)-BB30) * ALWMIP
      DB12=DB12 + WGHT(L) * FBB12(V,W)
C
      DC2 = DC2 + WGHT(L) * (FBC2(V,W)-BC20) / (1.-W)
      DC3 = DC3 + WGHT(L) * (FBC3(V,W)-BC30) * ALWMIP
      DC12=DC12 + WGHT(L) * FBC12(V,W)
C
      DD12=DD12 + WGHT(L) * FBD12(V,W)
C
      DE2 = DE2 + WGHT(L) * (FBE2(V,W)-BE20) / (1.-W)
      DE12=DE12 + WGHT(L) * FBE12(V,W)
C
      DJ2 = DJ2 + WGHT(L) * (FBJ2(V,W)-BJ20) / (1.-W)
      DJ12=DJ12 + WGHT(L) * FBJ12(V,W)
C
      DK12=DK12 + WGHT(L) * FBK12(V,W)
      endif
220   CONTINUE
c
c---- add here the various anomalous contributions
c
      if(KHIOR.ge.1) then
      BOR1(I,J) = BOR1(I,J) + VGHT(K)*(DBOR1+DBR1)
      BOR2(I,J) = BOR2(I,J) + VGHT(K)*(DBOR2+DBR2)
      BOR3(I,J) = BOR3(I,J) + VGHT(K)*(DBOR3+DBR3)
C
      BCH1(I,J) = BCH1(I,J) + VGHT(K)*(DBCH1+DBC1)
      BCH2(I,J) = BCH2(I,J) + VGHT(K)*(DBCH2+DBC2)
      BCH3(I,J) = BCH3(I,J) + VGHT(K)*(DBCH3+DBC3)
C
      BA1(I,J) = BA1(I,J) +VGHT(K)* (BA10 +BA20*ALMIN+0.5*BA30*ALMIN**2)
      BA2(I,J) = BA2(I,J) +VGHT(K)* DA2
      BA3(I,J) = BA3(I,J) +VGHT(K)* DA3
      BA12(I,J)= BA12(I,J)+VGHT(K)* DA12
C
      BB1(I,J) = BB1(I,J) +VGHT(K)* (BB10 +BB20*ALMIN+0.5*BB30*ALMIN**2)
      BB2(I,J) = BB2(I,J) +VGHT(K)* DB2
      BB3(I,J) = BB3(I,J) +VGHT(K)* DB3
      BB12(I,J)= BB12(I,J)+VGHT(K)* DB12
C
      BC1(I,J) = BC1(I,J) +VGHT(K)* (BC10 +BC20*ALMIN+0.5*BC30*ALMIN**2)
      BC2(I,J) = BC2(I,J) +VGHT(K)* DC2
      BC3(I,J) = BC3(I,J) +VGHT(K)* DC3
      BC12(I,J)= BC12(I,J)+VGHT(K)* DC12
C
      BD12(I,J)= BD12(I,J)+VGHT(K)* DD12
C
      BE1(I,J) = BE1(I,J) +VGHT(K)* (BE10 + BE20*ALMIN)
      BE2(I,J) = BE2(I,J) +VGHT(K)* DE2
      BE12(I,J)= BE12(I,J)+VGHT(K)* DE12
C
      BJ1(I,J) = BJ1(I,J) +VGHT(K)* (BJ10 + BJ20*ALMIN)
      BJ2(I,J) = BJ2(I,J) +VGHT(K)* DJ2
      BJ12(I,J)= BJ12(I,J)+VGHT(K)* DJ12
C
      BK12(I,J)= BK12(I,J)+VGHT(K)* DK12
      endif
210   CONTINUE
C
C          BORN TERMS        SMRET NROB         SMET NROB
C
       if(KHIOR.ge.0) then
      GLUQU(I,J) = ANORBO*CN * (1./NC) *(ALFA/PT2**2) * GLUQU(I,J)
      QUQUB(I,J) = ANORBO*CN *(2*CF/NC)*(ALFA/PT2**2) * QUQUB(I,J)
      BORN(I,J) = GLUQU(I,J) + QUQUB(I,J)
C
      GLUCH(I,J) = ANORBO*CN * (1./NC) *(ALFA/PT2**2) * GLUCH(I,J)
      CHCHB(I,J) = ANORBO*CN *(2*CF/NC)*(ALFA/PT2**2) * CHCHB(I,J)
      BCHA(I,J) = GLUCH(I,J) + CHCHB(I,J)
      endif
C
      if(KHIOR.ge.1) then
      BOR1(I,J) =   ANORM * BOR1(I,J)
      BOR2(I,J) =   ANORM * BOR2(I,J)
      BOR3(I,J) =   ANORM * (CF/NC) * BOR3(I,J)
      BORLL(I,J) = BOR1(I,J)+BOR2(I,J)+BOR3(I,J)
C     BORLL(I,J) =          +BOR2(I,J)+BOR3(I,J)
C
      BCH1(I,J) =   ANORM * BCH1(I,J)
      BCH2(I,J) =   ANORM * BCH2(I,J)
      BCH3(I,J) =   ANORM * (CF/NC) * BCH3(I,J)
      BCHLL(I,J) = BCH1(I,J)+BCH2(I,J)+BCH3(I,J)
C     BCHLL(I,J) =          +BCH2(I,J)+BCH3(I,J)
C
       endif
       if(KHIOR.ge.2) then
      BA1(I,J) = ANORM * BA1(I,J)
      BA2(I,J) = ANORM * BA2(I,J)
      BA3(I,J) = ANORM * BA3(I,J)
      BA12(I,J)= ANORM * BA12(I,J)
C
      BB1(I,J) = ANORM * BB1(I,J)
      BB2(I,J) = ANORM * BB2(I,J)
      BB3(I,J) = ANORM * BB3(I,J)
      BB12(I,J)= ANORM * BB12(I,J)
C
      BC1(I,J) = ANORM * BC1(I,J)
      BC2(I,J) = ANORM * BC2(I,J)
      BC3(I,J) = ANORM * BC3(I,J)
      BC12(I,J)= ANORM * BC12(I,J)
C
      BD12(I,J)= ANORM * BD12(I,J)
C
      BE1(I,J) = ANORM * BE1(I,J)
      BE2(I,J) = ANORM * BE2(I,J)
      BE12(I,J)= ANORM * BE12(I,J)
C
      BJ1(I,J) = ANORM * BJ1(I,J)
      BJ2(I,J) = ANORM * BJ2(I,J)
      BJ12(I,J)= ANORM * BJ12(I,J)
C
      BK12(I,J)= ANORM * BK12(I,J)
C
      BA(I,J) = BA1(I,J)+BA2(I,J)+BA3(I,J)+BA12(I,J)
      BB(I,J) = BB1(I,J)+BB2(I,J)+BB3(I,J)+BB12(I,J)
      BC(I,J) = BC1(I,J)+BC2(I,J)+BC3(I,J)+BC12(I,J)
      BD(I,J) =                            BD12(I,J)
      BE(I,J) = BE1(I,J)+BE2(I,J)         +BE12(I,J)
      BJ(I,J) = BJ1(I,J)+BJ2(I,J)         +BJ12(I,J)
      BK(I,J) =                            BK12(I,J)
C
      BAS(I,J) = BA1(I,J)+BA2(I,J)+BA3(I,J)
      BBS(I,J) = BB1(I,J)+BB2(I,J)+BB3(I,J)
      BCS(I,J) = BC1(I,J)+BC2(I,J)+BC3(I,J)
      BES(I,J) = BE1(I,J)+BE2(I,J)
      BJS(I,J) = BJ1(I,J)+BJ2(I,J)
      BS(I,J) = BAS(I,J)+BBS(I,J)+BCS(I,J)+BES(I,J)+BJS(I,J)
      endif
c
c  add leading log charm contribution to lead log light quarks
c
      HIGOR=0.D0
      ANOMA=0.D0
      BBORN=BORN(I,J) + BCHA(I,J)
      if(KHIOR.ge.2) then
      HIGOR=BA(I,J)+BB(I,J)+BC(I,J)+BD(I,J)+BE(I,J)+BJ(I,J)+BK(I,J)
      endif
      if(KHIOR.ge.1) then
      ANOMA=BORLL(I,J) + BCHLL(I,J)
      endif
C     IF (INTEGY.NE.1) write(io,107) PT,Y,CC1,Q1S
      write(itest,107) PT,Y,CC1,Q1S
107   FORMAT(' PT=',F9.3,'   Y=',F9.3,'   CC1=',D9.3,
     1'   Q1S=',D9.3)
C
      IF(IOPT.EQ.0) THEN
C -------- INVERT THE STEVENSON EQUATION FOR ALFAS  --------------------
c           ALQ2S=LOG(Q1S/RLAM2)
            ALQ2S=LOG(Q3S/RLAM2)
c!!	    
c	write(io,*) ' Q3S ',Q3S, ' RLAM2 ',RLAM2	    
c           AA=DZEROX(1.D-2,.25D0,.1D-4,1000000,STEVE,2)
            AA=ALFLUC(iloop,Q3S)/PI
c	 write(io,*) ' AA from STEVE',AA 
c	 write(io,*) ' alfas ',AA*PI
c!!	  
c          call rzero(.001,.45,aa,zero,.005,100,steve)
      ELSE
C -------- CALCULATE THE OPTIMAL COUPLING  -----------------------------
           IF (OPTI(1.D-2)*OPTI(.25D0).GT.0.D0) THEN
           write(io,102)
           write(io,*) 'NO OPTIM OF SIG-TOTAL FOR PT=',PT,' Y=',Y
           write(io,102)
           write(io,102)
           GO TO 201
           ENDIF
           AA=DZEROX(1.D-2,.25D0,.1D-4,1000000,OPTI,2)
	   write(io,*)' AA from OPTI ',AA
c           call rzero(.001,.45,aa,zero,.005,100,opti)
           ALQ2S=2.*(1./AA+AC*LOG(.5D0*AA*AB/(1.D0+AA*AC)))/AB
C ----- END OF OPTIMIZATION   -------------------------------------
      ENDIF
c inject 'optimal' coupling in iopt=0 option
c      AA=AAA(IJKL)
      Q2S=RLAM2*EXP(ALQ2S)
c     write(io,*) 'Q2S now is ',Q2S
      TTOTT=AA*((1.D0+.5D0*AA*AB*ALQ2S)*BBORN+AA*(HIGOR+ANOMA))
c      write(io,*)'***TTOTT**',TTOTT
      IF(TTOTT-TOTOLD.GE.0.) ISTOP=ISTOP+1
c
c 	write in file the relevant parameters for optimization
c 	combining with jenfilip program
c
       if(i68.eq.1) then
 	write(68,*)'Pt Y log(QDS) log(Q1S) bborn higor anoma '
      WRITE(68,*) PT,Y,LOG(QDS),LOG(Q1S),BBORN,HIGOR,ANOMA
      endif
c      WRITE(67,*) PT,Y,LOG(QDS),LOG(Q1S),BBORN,HIGOR,ANOMA
c
c
C  light quark and charm normalization contribution
c
      BORANN=AA*QUQUB(I,J)
      BORCOM=AA*GLUQU(I,J)
      ANOMA=AA**2*ANOMA
      if(KHIOR.ge.2) then
      HIGOR=AA*(.5*AA*AB*ALQ2S*BBORN+AA*HIGOR)
      endif
      BBORN=AA*BBORN
      if(KHIOR.eq.0) then
      TTOTT=BBORN
      elseif(KHIOR.eq.1) then
      TTOTT=BBORN+ANOMA
      endif
c
      BBCHA=AA*BCHA(I,J)
      ANOCHA=AA**2*BCHLL(I,J)
      TTCHA=BBCHA+ANOCHA
      BCH1(I,J)=AA**2*BCH1(I,J)
      BCH2(I,J)=AA**2*BCH2(I,J)
      BCH3(I,J)=AA**2*BCH3(I,J)
c
c   select the opt value of ttott
c
      IF(TTOTT.LT.TOTMIN) THEN
             TOTCHA=TTCHA
             TOTMIN=TTOTT
             TOTBOR=BBORN
             TOTANO=ANOMA
             TOTHIG=HIGOR	     
      ENDIF
      TOTOLD=TTOTT
c
      CCORR=TTOTT/BBORN
         IF(ISIGM.eq.1)then
	   FACIN=PI
c	   write(io,*) ' DSIG/DY/D2PT' 
	 ELSEIF(ISIGM.eq.2) then
	   FACIN=1.D0 
c	   write(io,*) ' E DSIG/D3P'
	 ELSEIF(ISIGM.eq.3) then
	   FACIN=PI*2.D0*PT
c           write(io,*) ' DSIG/DY/DPT'
	 ENDIF 
      IF (INTEGY.NE.1) THEN    
        write(itest,104) AA,FACIN*TTOTT,FACIN*BBORN,CCORR,
     1 FACIN*ANOMA,ANOMA/TTOTT
104   FORMAT(' AA=',D9.3,'   SIGTOT=',D9.3,'   SIG-BORN=',D9.3
     1,'   K-FACTOR=',D9.3
     2,'      .ANOMA.=',D9.3,'   ANOMA/TOTAL =',D9.3)
c     2,'   BORN-ANNIHI.=',D9.3,'   BORN-COMPTON=',D9.3)     
c        write(io,104) AA,FACIN*TTOTT,FACIN*BBORN,CCORR,
c     1  FACIN*ANOMA,FACIN*(TTOTT-ANOMA)
C     write(io,104) AA/AA,FACIN*TTOTT/AA**2,FACIN*BBORN/AA,CCORR,
C    1  FACIN*ANOMA/AA**2,FACIN*(TTOTT-ANOMA-BBORN)
c    printing of leading log. charm contributing processes
c        write(io,108) TTCHA,BBCHA,BCH1(I,J),BCH2(I,J),BCH3(I,J)
c108   FORMAT(' THE CHARM :    CHA-LL=',D9.3,'   CHA-BORN=',D9.3
c     1,' ANO GLU-GLU=',D9.3,'  ANO GLUE-QRK=',D9.3
c     2,'  ANO QURK-QURK=',D9.3)
        write(itest,102)
	iioo=io
	if(i66.eq.1) iioo=66
	if(i66.eq.2)
     1WRITE(iioo,*) 'PT    QDS      AA log(Q1S)    ',
     1	'ANOMA      TTOTT-ANOMA '
         if(i66.eq.1) 	
     1WRITE(iioo,*) PT,QDS,AA,LOG(Q1S),ANOMA,TTOTT-ANOMA
      ENDIF
      IF(ISTOP.GE.15) GO TO 201
1     CONTINUE
201   CONTINUE
c
      IF(INTEGY.NE.1) THEN
        TOTINT=TOTMIN
	BORINT=TOTBOR
	ANOINT=TOTANO
	HIGINT=TOTHIG
       if(ichi2.eq.1) then
        TOTMIN=CORAP(I)*TOTMIN
        CHI2P=(TOTMIN-SIGEX(I))**2/ERREX(I)**2
        IF(TOTMIN.GT.1D+5) CHI2P=0.
        CHI2=CHI2+CHI2P
        DADTH=SIGEX(I)/TOTMIN
        ERRTH=ERREX(I)/TOTMIN
       write(io,105) IPRO,PT,Y,TOTMIN,DADTH,ERRTH
       write(io,106) IPRO,PT,XT,CHI2P,CHI2
105   FORMAT(' ###IPRO=',I4,' PT=',F9.3,'  Y=',F9.3,'    SIGTOT=',D11.5,
     1'  DAT/TH=',E11.5,' ERR/TH=',E11.5)
106   FORMAT(' ###IPRO=',I4,' PT=',F9.3,' XT=',F9.3,'  CHI2/PNT=',D11.5,
     1' CHI2TOT=',D11.5)
c        else 
c	write(io,*) ' PT ',PT,' Y ',Y,' SIGTOT ',TOTMIN
       endif
      ELSE
      		IF(IXBOND.EQ.0) THEN
      		TOTINT=TOTINT+TOTMIN*YGHT(J)
      		BORINT=BORINT+TOTBOR*YGHT(J)
      		ANOINT=ANOINT+TOTANO*YGHT(J)
      		HIGINT=HIGINT+TOTHIG*YGHT(J)
c      write(io,*) 'DBUG,IPR,PT,RAP,TTMIN,TTINT',IPRO,PT,Y,TOTMIN,TOTINT
		ELSE
      		TOTINT=TOTINT+TOTMIN*YGHT(J)*XT*CHY
      		BORINT=BORINT+TOTBOR*YGHT(J)*XT*CHY
      		ANOINT=ANOINT+TOTANO*YGHT(J)*XT*CHY
      		HIGINT=HIGINT+TOTHIG*YGHT(J)*XT*CHY
c      write(io,*) 'DBUG,IPR,PT,RAP,TTMIN,TTINT',IPRO,PT,Y,TOTMIN,TOTINT
		ENDIF
      ENDIF
c end loop on rapidity
200   CONTINUE
c
      IF(INTEGY.EQ.1) THEN
      		IF(IXBOND.EQ.0) THEN
	          TOTINT=TOTINT/(YSUP-YINF)
	          BORINT=BORINT/(YSUP-YINF)
	          ANOINT=ANOINT/(YSUP-YINF)
	          HIGINT=HIGINT/(YSUP-YINF)
                  write(itest,102)
	          IF(ISIGM.eq.1) THEN
       write(itest,*) ' PT=',PT,' DSIG/DY/D2PT AVE. IN RAP. X-SECTION='
     1    ,TOTINT*FACIN
                    ELSEIF(ISIGM.eq.2) THEN
      write(itest,*) ' PT=',PT,' E DSIG/D3P AVE. IN RAP. X-SECTION='
     1    ,TOTINT*FACIN	  
	            ELSEIF(ISIGM.eq.3) THEN
      write(itest,*) ' PT=',PT,' DSIG/DPT/DY AVE. IN RAP. X-SECTION='
     1    ,TOTINT*FACIN	  
	          ENDIF
          write(itest,*) ' PT=',PT,' BOR=',BORINT*FACIN,
     1    ' AMOMA=',ANOINT*FACIN
     1    ,' HIG.ORD.=',HIGINT*FACIN
          write(itest,102)
	  if(i67.eq.1) then
          WRITE(67,*) 'PT DSIG/DY/D2PT INTEGR. IN Y TOTINT .1*TOTINT'
	  WRITE(67,*) PT, TOTINT, .1*TOTINT
	  endif
		ELSE
	          TOTINT=TOTINT/(XSUP-XINF)
	          BORINT=BORINT/(XSUP-XINF)
	          ANOINT=ANOINT/(XSUP-XINF)
	          HIGINT=HIGINT/(XSUP-XINF)
          write(itest,102)
          write(itest,*) ' PT=',PT,
     1' DSIG/DY/D2PT INTEGR. IN XF. X-SECTION=',TOTINT
          write(itest,*) ' PT=',PT,' BOR=',BORINT,' AMOMA=',ANOINT
     1    ,' HIG.ORD.=',HIGINT
          write(itest,102)
	  if(i67.eq.1) then
	  WRITE(67,*) 'PT  DSIG/DY/D2PT INTEGR TOTINT .1TOTINT'
	  WRITE(67,*) PT, TOTINT, .1*TOTINT
	  endif
		ENDIF
      ENDIF
c end loop on pt      
100    CONTINUE
c end loop on the fragmentation scale
900    CONTINUE
c end loop on the coupling aa=alphas/pi values 
9999  CONTINUE
102    FORMAT (' ')
      RETURN
      END
c
	double precision function alfluc(iloop,scale2)
	implicit real*8(a-h,o-z)
	alfluc = alfas(iloop,scale2)
	return
	end
