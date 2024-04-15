	SUBROUTINE PARAM
c choice of various parameters 
	IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8 MASCH2,MASBO2,MASTO2,LAMBDA2,N
	CHARACTER*20 PARM(20)
	CHARACTER*128 path_bsfile
       DIMENSION VALUE(20)  
       COMMON/OUTPUT_UNITS/IO,IER,ITEST  
       COMMON/ANOMAL/IANOM,IANORD
       COMMON/phoord/lhior 
       include "conv.f"
       COMMON/EVOPT/PTEVO(19)
       COMMON/CONS/PIC,GS,GV,GW,N,GTR,CFC,PT2,VC
      COMMON/VALU/JMAR,IPT
      COMMON/PTVALUE/PT
      COMMON/HADR/IH1,IH2 
      COMMON/XSECT/ISIGM,ILOOP,IHIOR
      COMMON/CHOI/HC2,ZRAP,RS,ZAL,CM,CMU,CMP
      COMMON/ALFA/MASCH2,MASBO2,MASTO2,LAMBDA2
c *** common blocks for use with pdflib tables
      COMMON/W50511/NPTYPE,NGROUP,NSET,MODE,NFL,LO,TMAS
      COMMON/W50512/QCDL4,QCDL5
      COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX
      include "pdfmod.f"
c      COMMON/PDFMOD/NPTYPRO,NGROPRO,NSETPRO,NPTYTAR,NGROTAR,
c     * NSETTAR
      include "pdfchoix.f"
c      LOGICAL LMRS98,LMRS99,LCTEQ5M,LPDFLIB,LGRID
c      COMMON/PDFCHOIX/LMRS98,LMRS99,LCTEQ5M,LPDFLIB,LGRID
c      COMMON/ABFOW1/IMODE
      COMMON/YRANGE/YMIN,YMAX
c      COMMON/NUCLEI/ZSURA
      COMMON/EXPERIENCE/IEXP 
      common/INTCHOIX/lintegy
      common/BAPRECISION/ACCU,JTMX1,JTMX2,JCALL1,JCALL2,JCALL3
      common/baseio/lu
      COMMON/INTEGPREC/NVINT,NWINT,NY2
c OWHPOPI uses  gaussian integration
c number of integration points for V integration (32 recommended)
c number of integration points for W integration (32 recommended)
c number of integration points for rapidity integration (8 may be enough)
      common/pion/ih3
        include "projec.f"
	include "qcd.f"
	CHARACTER*128 name_experiment
	common/name/name_experiment
	common/lname/ilname
	logical lintegy
c used in fonfra
c	logical lzmean
c        common/zmean/lzmean
c	lzmean=.false.
	save
	lmrs98 = .false.
	lmrs99 = .false.
	lcteq5m = .false.
	lpdflib =.false.
	lgrid =.false.
	in=55
	io=69
	ier=15
	itest=8
	lu=8
C open input  file with user's parameters
	OPEN(UNIT=in,file='param.dat',status='OLD',ERR=999)
c path where the ouput files will be stored 
	read (in,*) path_bsfile
	ilen = icharnum(path_bsfile)
c	OPEN(UNIT=ier,file='error.out',status='unknown')
c	OPEN(UNIT=itest,file='test.out',status='unknown')
c	OPEN(UNIT=io,file='result.out',status='unknown')	
	OPEN(UNIT=ier,file=PATH_BSFILE(1:ILEN)//'error.out',
     #	STATUS='unknown')
	OPEN(UNIT=itest,file=PATH_BSFILE(1:ILEN)//'test.out',
     #	STATUS='unknown')
     	OPEN(UNIT=io,file=PATH_BSFILE(1:ILEN)//'result.out',
     #	STATUS='unknown')      	
c        OPEN(UNIT=6,file=PATH_BSFILE(1:ILEN)//'file.out',
c     #	STATUS='unknown')    
c	write(itest,*) 'ilen=',ilen
c	write(itest,*) 'path_bsfile=',path_bsfile(1:ilen)
c ih1 and ih2 define the initial particle
C ih = 0: proton; ih = 1: antiproton, ih = 2: photon, 
C ih = 3: pion; ih = 4: nuclei
c ih2 is the target particle, ih1 is the incident particle
c Type of projectile?
	read (in,*) ih1
	write(itest,*) 'projectile',ih1
c	if (ih1.eq.3) then
	 read (in,*) ZPRO
	 read (in,*) APRO
	 write(itest,*) ' Z=',ZPRO
	 write(itest,*) ' A=',APRO
c  	endif
c determination (pdflib format) of distribution functions in hadron h1
C NTYPE, NGROUP, NSET  
	read (in,*) NPTYPRO
	write(itest,*) 'NPTYPRO=',NPTYPRO
	PARM(1) = 'NPTYPE'
	VALUE(1) = NPTYPRO
	read (in,*) NGROPRO
	write(itest,*) 'NGROPRO=',NGROPRO
	PARM(2) = 'NGROUP'
	VALUE(2) = NGROPRO
	read (in,*) NSETPRO
	write(itest,*) 'NSETPRO=',NSETPRO
	PARM(3) = 'NSET'
	VALUE(3) = NSETPRO
	if(ngropro.eq.98) then
	  lmrs98=.true.
	 elseif(ngropro.eq.99) then
	  lmrs99=.true.
	 elseif(ngropro.eq.100) then
	  lcteq5m=.true. 
	 elseif(ngropro.eq.1000) then
	  lgrid=.true. 
	 else
	   lpdflib=.true. 
	 endif 
c!!!!	 if(lpdflib) CALL PDFSET(PARM,VALUE)
c Type of target?
	read (in,*) ih2
	write(itest,*) 'target',ih2
c	if (ih2.eq.3) then
	 read (in,*) ZCIB
	 read (in,*) ACIB
	 write(itest,*) ' Z=',ZCIB
	 write(itest,*) ' A=',ACIB
	 ZSURA=ZCIB/ACIB
c  	endif
c determination (pdflib format) of distribution functions in hadron h2
	read (in,*) NPTYTAR
	write(itest,*) 'NPTYTAR=',NPTYTAR
	PARM(1) = 'NPTYPE'
	VALUE(1) = NPTYTAR
	read (in,*) NGROTAR
	PARM(2) = 'NGROUP'
	VALUE(2) = NGROTAR
	write(itest,*) 'NGROTAR=',NGROTAR
	read (in,*) NSETTAR
	write(itest,*) 'NSETTAR=',NSETTAR
	PARM(3) = 'NSET'
	VALUE(3) = NSETTAR
c choose the type of produced particle (H1 + H2 --> H3)
	read (in,*) ih3
	write(itest,*) 'IH3=',ih3
c
	read (in,*) LHIOR
c forces production of hadron at nll if final sate is not a photon	
	if (ih3.ne.0) then
	  lhior = 4
	endif  
	write(itest,*) 'LHIOR=',LHIOR
c LHIOR=0 BORN ONLY
c LHIOR=1 BORN and ANOMA
c LHIOR=2 NLL ABFS
c LHIOR=3 NLL ABFS + NLL FRAGMENTATION	
c LHIOR=4 NLL FRAGMENTATION only	
C IHIOR=1 LOWEST ORDER ONLY, IHIOR=2 HIGHER ORDER  IN PIOINCL
C  iord in structure function
c iord = 0 LL 
c iord = 1 NLL
	read (in,*) IHIOR
	write(itest,*) 'IHIOR=',IHIOR
	iord = 1
	write(itest,*)' iord ',iord
	
c anomalous photon structure function
c	read (in,*) ianom
	ianom = 2
	write(itest,*)' ianom=',ianom
	if(ianom.lt.0 .or. ianom.gt.2) then
	write(ier,*)' anomalous order not implemented ianom ',ianom
	stop
	endif
c	if(ianom.eq.1) ih3=21
	read (in,*) ianord
	write(itest,*)' ianord=',ianord	
c	if(ianom.eq.2) then
c ianord determines the choice of fragmentation function from parton
c to photon: 1 Bourhis set I, 2 Bourhis set II
c or from partons to pion: 1 BKK , 2 KKP, 3 BFGW

c Number of active flavors
	read (in,*) jnf
	write(itest,*) 'jnf=',jnf
	nf=jnf
c  Quark masses (used in alpha_s to determine the number of active flavors
c  charm
	MASCH2 = (1.5D0)**2
c  bottom	
	IF (JNF.EQ.5) THEN
	  MASBO2 = (4.5D0)**2
	ELSE IF (JNF.EQ.4) THEN
	  MASBO2 = (1.D+10)**2
	ENDIF
c  top		
	MASTO2 = (1.D+10)**2
c choice of cross sections
C     ISIGM=1 ==> DSIGMA/DY/DPT2
C     ISIGM=2 ==> E*DSIGMA/D3P
C     ISIGM=3 ==> DSIGMA/DY/DPT
	read (in,*) isigm
	write(itest,*) ' isigm ',isigm
	if(isigm.lt.1 .or. isigm.gt.3) then
	write(ier,*)' not implemented cross section choice ',isigm
	stop
 	endif
c factor (hbar*c)**2
	read (in,*) hc2
	write(itest,*) 'hc2=',hc2
c number of loops in alpha_s AND IN THE EVOLUTION EQUATIONS 
C(DEPENDS ON THE CHOICE OF THE DISTRIBUTION FUNCTIONS)
	read (in,*) iloop
	write(itest,*) 'iloop=',iloop
	if(iloop.gt.2 .or. iloop.lt.1) then
	  write(ier,*)' iloop not implemented illop=',iloop
	  stop
	endif
c ischeme:  factorisation scheme 
c 2: MSBAR; 3:DIS
	read (in,*) jmar
	write(itest,*) 'jmar=',jmar
	if(jmar.gt.3 .or. jmar.lt.2) then
	  write(itest,*)' scheme not implemented jmar=',jmar
	  stop
	endif
c ici : to select between collider and fixed target
c ici = 0 collider; ici = 1 fixed target
	read (in,*) ici
	write(itest,*) 'ici=',ici
c center of mass energy for collider
c EBEAM for fixed target
	read (in,*) ebeam
c max. rapidity of produced photon
	read (in,*) ymax
	write(itest,*) 'ymax=',ymax
c min. rapidity of produced photon
	read (in,*) ymin
	write(itest,*) 'ymin=',ymin
c integrate in rapidity?
	read (in,*) lintegy
	write(itest,*) 'lintegy ',lintegy	
c pt values 
	read (in,*) ipt
	write(itest,*) 'nb of pt points=',ipt
	if(ipt.gt.19) then
	  write(ier,*)' too many pt points ',ipt
	  stop
	endif  
	do i=1,ipt
	 read (in,*) ptevo(i)
	 write(itest,*)' pt value ',ptevo(i)
	 if(ptevo(i).lt.4.D0) then
	 write(ier,*) ' too low pt value ',ptevo(i)
	 endif
c	 if(ptevo(i).gt.(rs/2.) then
c         write(ier,*)' too high pt value ',ptevo(i)
c	 endif
	enddo
	pt=ptevo(1)
c constants for SU(3) color
	N = 3.D0
	CFC = (N*N-1.D0)/(2.D0*N)
	GTR = DFLOAT(JNF)/2.D0
c get square root of s (center of mass energy)
	if (ici.eq.0) then
	  rs = ebeam
	  s = rs**2
	  write(itest,*) 'rs =',rs
	else if (ici.eq.1) then
	  xmp = .93828d0
	  xmpi = .139567d0
	  if (ih1.le.1) then
	    xm_inci = xmp	  
	   elseif (ih1.eq.3) then
	    xm_inci = xmp	    
	    write(itest,*) 'cross section per nucleon'
	  elseif (ih1.eq.4.or.ih1.eq.5) then
	    xm_inci = xmpi
	  else
	    write(ier,*) 'bad choice of incident particule ih1',ih1
	    stop
	  endif
	  if (ih2.le.1) then
	    xm_cible = xmp
	  elseif(ih2.eq.3) then
	    xm_cible = xmp
	    write(itest,*) 'cross section per nucleon'
	  else  
	    write(ier,*) 'bad choice of target particule ih2',ih2
	    stop
	  endif
	  rs = dsqrt(xm_inci**2+xm_cible**2+2.d0*xm_cible*ebeam)
	  s = rs**2
	  write(itest,*) 'ebeam =',ebeam
	  write(itest,*) 'rs =',rs
	endif
c scale choice: 1 c*pt; 2 c*pt*(1.-xt)
  	read (in,*) ichoi_scale
	write(itest,*) 'ichoi_scale=',ichoi_scale
	if(ichoi_scale.ne.1) then
	write(ier,*)' not implemented ichoi_scale ',ichoi_scale
	stop
	endif
c cm,cmu and cmf are the normalisation factors for the scale choice
c respectivly for initial factorisation, renormalisation and final factorisation
	read (in,*) cm
	write(itest,*) 'cm=',cm
	read (in,*) cmu
	write(itest,*) 'cmu=',cmu
	read (in,*) cmp
	write(itest,*) 'cmp=',cmp
C experiment name
	read (in,*) name_experiment
c	ilname = icharnum(name_experiment)	
c	write(8,*) 'ilname =',ilname
c	write(8,*) 
c     #	'name_experiment=',name_experiment(1:ilname)
c	write(*,*) 'ilname =',ilname 
c ZAL parameter not allowed to be changed
          ZAL=1.
c Number of point in the gauss integration for abfs
c It uses DGSET of the CERN library
c for integration on v and w
          NVINT=32
	  NWINT=32
c	nvint=8
c	nwint=8
c for integration on rapidity
	  NY2=8
c Precision for bases integration for pioincl
c Bases is an adaptative Monte-Carlo integrator, the user has to provide:
c the precision in percent
c the number of iterations (for the grid step and for the integration step)
c the number of calls per iteration
c The program ends when either the precision or the maximum number 
c of iteration is reached
c Precision in %
        ACCU = 1.D0
c number of iterations
c for grid
        JTMX1 = 5
c for integration
        JTMX2 = 5
c number of calls per iteration
c for Leading Order
        JCALL1 = 1000
c for Next to Leading Order: part in delta(1-w)
        JCALL2 = 2000
c for Next to Leading Order: other parts
        JCALL3 = 3000
	RETURN
999     write(6,*) 'ERROR OPENING FILE param.dat'
	STOP	
	END
	
