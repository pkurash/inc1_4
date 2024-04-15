CDECK  ID>, INITONE.
      subroutine initone
c interface to OWPHOPI (i.e. abfs.f/pholib) code 
	IMPLICIT REAL*8 (A-H,O-Z)     
c*** commons to handle both programs
	common/output_units/io,ier,itest
      COMMON/CHOI/HC2,ZRAP,ZRS,ZAL,CM,CMU,CMP
c HC2 for unit pb /GeV**2 
C ZRAP rapidity of produced particle, 
c ZRS center of mass energy, already set 
C ZAL dummy here  
C CM,CMU,CMP -> cc1,cc2,ccd see below scale_choice  
      COMMON/HADR/IH1,IH2
      COMMON/XSECT/ISIGM,ILOOP,IHIOR
      COMMON/YRANGE/YMIN,YMAX
      COMMON/INTCHOIX/lintegy
      logical lintegy
      COMMON/EVOPT/PTEVO(19)
      COMMON/PTVALUE/PTV
      COMMON/YVALUE/RAPIDITE
      COMMON/VALU/JMAR,IPT 
c      COMMON/NUCLEI/ZSURA
      include "pdfchoix.f"
c for one
	common/constants/cn
	common/scale_choice/cc1,cc2,ccd
c   q1s=cc1*pt2  factorization scale
c   q3s=cc2*pt2	 renormalization scale
c   qds=ccd*pt2  fragmentation scale
       include "projec.f"
c      COMMON/PROJEC/ZCIB,ACIB,ZPRO,APRO
c     *,XPT(20),XY(20)
c     *,YINF,YSUP,XINF,XSUP  
c      COMMON/PROJECI/IPRO,ITAR,ISTRUC
c     *,IPT1,IPT2,JY1,JY2,INTEGY,IXBOND
c projectile choice IPRO: 1=protons,2=antiprotons,3=pis +,4=pis - 6=nucleus
c target choice ITAR =0 (proton), =3 (nucleus) ZCIB,ACIB : Z and  A of the target
c ZPRO,APRO when ipro=6, Z and A of projectile
c structure function choice ISTRUC 
c    ISTRUC=1 MRS98,99 , =2 CTEQ5M; ISTRUC=3 or 4 pdflib, ISTRUC=5 tables
c    futher input for structure function via PDFMOD
c pt points XPT
c rapidity points XY
c pt points from XPT(IPT1) to XPT(IPT2)
c rapidity points from XY(JY1) to XY(JY2)
c INTEGY=0 calculates at given y point
C INTEGY=1 integrates in rapidity  between YINF and YSUP
c IXBOND=1 integrates instead in xf between XINF and XSUP
        include "pdfmod.f"
c       COMMON/PDFMOD/NPTYPRO,NGROPRO,NSETPRO,NPTYTAR,NGROTAR,
c     * NSETTAR
c  for pdflib, give nptype,nset,ngroup
c  nptypro,nsetpro,ngropro  for projectile
c  nptypar,ngrotar,nsettar  for target 
c	COMMON/ANOMAL/IANOM,IANORD
c ----- ianom=0 no qcd corrections ; ianom=1 leading log qcd approx
c ----- ianom=2, ianord=1 ll tables from chiappetta and guillet et al.
c ----- ianom=2, ianord=2 ntl tables from chiappetta and guillet et al	
        include "conv.f"
c	COMMON/CONVI/ IORD,ICONV
c       COMMON/CONV/OWLAM,OWLAM2,RLAM,RLAM2,Q02PR,Q02PI
c   structure functions conventions
c   iord=0                     leading order
c   iord=1             next to leading order
c   iconv=0             universal convention
c   iconv=1          quark distribution = f2
c   lambda values
c    OWLAM : lambda in SF
c    OWLAM2=OWLAM**2
C    RLAM : lambda in alpha_s
c    RLAM2=RLAM**2 
c    Q02PR : Q0**2 min for proton SF
c    Q02PI : Q0**2 min for pion SF  
	include "kine.f"   
c      COMMON/KINE/S,RS,PT,PT2,Y,VV,WW,XX1,XX2,XX1M
c S RS**2
c RS center of mass energy in GeV 
c PT pt in GeV/c
c PT2 pt**2
c Y rapidity 
c VV,WW,XX1,XX2,XX1MIN internal kinematics variable
c only RS need to be initialized here
      COMMON/INTEGPREC/NVINT,NWINT,NY2
c this program will use gaussian integration
c number of integration points for V integration (32 recommended)
c number of integration points for W integration (32 recommended)
c number of integration points for rapidity integration (8 may be enough)     
c --- choice of accuracy (leading log, next to leading)
c	iord=0/1 LL/NLL; IHIOR=0/1/2 born/LL/NLL
c         iord=IHIOR-1
c	iconv=0/1 universal/quark dis = f2 
C JMAR FACTORIZATION SCHEME DEPENDENCE:
C 1= FACTORIZATION SCHEME CQ=1 (C.F. NUCL PHYS. B327,105
C 2=MSBAR FACTORIZATION SCHEME (TO USE WITH HMRS, MT.....)
C 3=DIS FACTORIZATION SCHEME (TO USE WITH DFLM, MT, ....)
	 save
          iconv=JMAR-2	  
c --- choice of projectile
        if(ih1.eq.0) then
	    ipro=1
	elseif(ih1.eq.1) then
	    ipro=2
	elseif(ih1.eq.2) then
	     write(ier,*)' photon not implemented as incident particle'
	     stop
	elseif(ih1.eq.4) then
	     ipro=3
	     write(ier,*)' pi+'
	elseif(ih1.eq.5) then
	     ipro=4  	     
	     write(ier,*)' pi-'   
	elseif(ih1.eq.3) then
	     ipro=6
	     write(itest,*)
     #' nucleon as projectile with Z',ZPRO,' and A',APRO	     
c	     write(io,*)' heavy nucleon as projectile Z',ZPRO,' and A',APRO
	endif                       
c ---- choice of target
c ----- ZCIB,ACIB : Z AND A OF THE TARGET
        if(ih2.eq.0) then
       	  ACIB=1.D0
	  ZCIB=1.D0
	elseif(ih2.eq.1) then
	  write(ier,*)' antiproton not implemented as target'
	  stop
	elseif(ih2.eq.2) then  
	  write(ier,*)' photon not implemented as target particle'
	  stop
        elseif(ih2.eq.4.or.ih2.eq.5) then
	  write(ier,*)' pion not implemented as target particle'
	  stop
	elseif(ih2.eq.3) then
	   write(itest,*)' nucleon as target with Z',ZCIB,' and A',ACIB
	endif
	itar=ih2
c --- choice of structure function
       if(lmrs98) then
           ISTRUC=1
       elseif(lmrs99) then
           ISTRUC=2
       elseif(lpdflib) then
           ISTRUC=4
       elseif(lcteq5m) then
            ISTRUC=3 
       elseif(lgrid) then
            ISTRUC=5	    	   	  	
	endif
c fill PROJEC
	ipt1=1
	ipt2=1
	call ucopy(ptevo,xpt,19*2)
      	YINF=YMIN
      	YSUP=YMAX
	jy1=1
	jy2=1
c	xy(jy1)=rapidite
	XY(JY1)=(YMAX+YMIN)/2.D0
	rapidite=XY(JY1)
c	write(io,*)' INITONE rapidity',rapidite,' xy(1)',xy(1)
	if(lintegy) then
	   integy=1
	 else
	   integy=0
	 endif    	  
c scale_choice
	cc1=CM**2
      	cc2=CMU**2
      	ccd=CMP**2 
c	write(io,*)' IN INITONE cc1=',cc1,' cc2=',cc2,' ccd=',ccd    
c need to initialize only RS = sqrt(s) here
        RS=ZRS
c	WRITE(io,*)' IN INITONE RS is',RS
c full sf initialization	
        call onestru
c choice of fragmentation functions
	  call onefra
c choice of cross-section 
	cn=hc2
C following done elsewhere	
c	IF(ISIGM.EQ.1) THEN
C  DSIGMA/DY/DPT2
c          FACIN=PI
c        ELSE IF (ISIGM.EQ.2) THEN
C E*DSIGMA/D3P	
c          FACIN=1.D0
c        ELSE IF (ISIGM.EQ.3) THEN
C  DSIGMA/DPT/DY
c          FACIN=PI*2.D0*PTV
c        ENDIF
	return
	end
