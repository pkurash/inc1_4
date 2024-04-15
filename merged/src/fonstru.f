       SUBROUTINE FONSTRU(ITI,X,IH,Q2,XUH,XUBH,XDH,XDBH,XSH,XCH,
     # XBH,XTH,XGPRO)
C returns X times Up, UBAR, Down, DBAR, Strange,Charm,
C                 Bottom, Top ,Gluon  distribution functions
C at scale Q2 in proton (IH=0) or anti-proton(IH=1) 
C as respectively XUH, XUBH, XDH, XDBH, XSH, XBH, XTH, XGPRO
C ITI =1: projectile,  ITI=2, TARGET
C      IMPLICIT REAL*8 (A-H,L-Z)
C warning, the whole pioincl has implicit real*8(L-Z)!!
	IMPLICIT REAL*8 (A-H,O-Z)
c	COMMON/NUCLEI/ZSURA
        CHARACTER*20 PARM(20)
        DIMENSION VALUE(20)
        include "pdfchoix.f"
        include "pdfmod.f"
	include "projec.f"
	save 
	xsave=x
c	write(65,*) 'in fonstru x',x,'q2',q2,'iti',iti,' ih ',ih
c	save    
c       write(65,*) 'in fonstru x',x
c       write(65,*) 'fonstru lmrs98 ',lmrs98,' lmrs99 ',lmrs99
c       write(65,*) 'fonstru lcteq5m ',lcteq5m,' lpdflib ',lpdflib
	QQ=DSQRT(Q2)
C  DISTIBUTION DANS LE HADRON H (IH=0 H=P , IH=1 H=PB)
c	if (lmrs99) then
c          if(ngropro.eq.99) then
c	    imode = nsetpro
c	   elseif(ngrotar.eq.99) THEN
c	     imode = nsettar
c	   endif
c	write(6,*)'MRS99 set',imode
c	call mrs99(x,qq,imode,upv,dnv,usea,dsea,str,chm,bot,gl)
c	 top =0.d0
c	elseif (lmrs98) then
c          if(ngropro.eq.98) then
cc	    imode = nsetpro
c	  elseif(ngrotar.eq.98) THEN
c	    imode = nsettar
c	  endif
c	write(6,*)'MRS98 set',imode
c	call mrs98(x,qq,imode,upv,dnv,usea,dsea,str,chm,bot,gl)
c	 top =0.d0 	       
c	else if (lcteq5m) then
c	 Iset=nsetpro
c   	 call SetCtq5(Iset)
c	 USEA=X*Ctq5Pdf (-1, X, QQ)
c	 DSEA=X*Ctq5Pdf (-2, X, QQ)
c	 UPV=X*Ctq5Pdf (1, X, QQ) -USEA
c	 DNV=X*Ctq5Pdf (2, X, QQ)-DSEA
c	 STR=X*Ctq5Pdf (3, X, QQ)
c	 CHM=X*Ctq5Pdf (4, X, QQ)
c	 BOT=X*Ctq5Pdf (5, X, QQ)
c	 GL=X*Ctq5Pdf (0, X, QQ)
c	 top=0.D0
c	elseif(lpdflib) then
C CHOICE OF DISTRIBUTION FUNCTIONS FOLLOWING PDFLIB
          PARM(1) = 'NPTYPE'
          PARM(2) = 'NGROUP'
          PARM(3) = 'NSET'

c	if (lcteq4m) then
c	  VALUE(1) = 1
c	  VALUE(2) = 4
c	  VALUE(3) = 33 +NsetPro
c         endif      
c	  VALUE(1)=NPTYPRO
c	  if (lmrs98) then
c	    VALUE(2)=3 
c	    VALUE(3)=NSETPRO+66
c	  elseif (lmrs99) then
c	    VALUE(2)=3
c	    VALUE(3)=NSETTAR+88
c	  elseif(lcteq5m) then
c	    VALUE(2)=4
c	    VALUE(2)=46+NSETTAR
c	  else	
	IF (ITI.eq.2) THEN
	    VALUE(1)=NPTYTAR
	    VALUE(2)=NGROTAR
	    VALUE(3)=NSETTAR
	elseif(ITI.eq.1) THEN
	    VALUE(1)=NPTYPRO
	    VALUE(2)=NGROPRO
	    VALUE(3)=NSETPRO  
	endif
c setting of lambda4 for cteq4   as it does not work automatically
	if(value(2).eq.4) then
	  if(value(3).eq.34) then
	  PARM(4)='QCDL4'
	  value(4)=0.296D0
	  elseif(value(3).eq.35) then
	  PARM(4)='QCDL4'
	  value(4)=0.213D0
	  elseif(value(3).eq.36) then
	  PARM(4)='QCDL4'
	  value(4)=0.253D0
	  elseif(value(3).eq.38) then
	  PARM(4)='QCDL4'
	  value(4)=0.344D0
	  elseif(value(3).eq.39) then
	  PARM(4)='QCDL4'
	  value(4)=0.399D0
	  elseif(value(3).eq.40) then
	  PARM(4)='QCDL4'
	  value(4)=0.302D0
	endif
	endif
c**standard	
	if(value(1).ne.4) then
c	write(65,*)'before pdfset parm(1)',parm(1),' value(2)',value(2)
        CALL PDFSET(PARM,VALUE)
c	write(65,*)'after pdfset parm(1)',parm(1),'value(2)',value(2)
c	write(65,*)'before structm x ',x,' qq ',qq,' iti ',iti
	CALL STRUCTM(X,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
c	write(65,*)'after structm x ',x,' qq ',qq
c	write(65,*)' UPV ',upv,'DNV',dnv,' USEA ',usea,' DSEA ',dsea
	endif
      IF (IH.EQ.0.OR.IH.EQ.1) THEN                       
        XUH=UPV*DFLOAT(1-IH)+USEA                                         
        XUBH=UPV*DFLOAT(IH)+USEA                                           
        XDH=DNV*DFLOAT(1-IH)+DSEA                                          
        XDBH=DNV*DFLOAT(IH)+DSEA 
c	write(65,*)' XUH ',xuh,' XDH ',xdh                                           
        XSH=STR
        XCH=CHM
        XBH=BOT
        XTH=TOP
        XGPRO=GL
      ELSE IF (IH.EQ.3) THEN
c** ih=3 nucleon 
	if(VALUE(1).eq.1) then
c**using proton SF
	if(ITI.eq.1) then
	ZSURA=ZPRO/APRO
	 elseif(ITI.eq.2) then
	ZSURA=ZCIB/ACIB 
	 endif 
c	 write(65,*)' ITI ',ITI,' ZSURA ',ZSURA   
        XUH= ZSURA*(UPV+USEA) + (1.D0-ZSURA)*(DNV+DSEA)                                         
        XUBH= ZSURA*USEA + (1.D0-ZSURA)*DSEA                                       
        XDH= ZSURA*(DNV+DSEA) + (1.D0-ZSURA)*(UPV+USEA)                                         
        XDBH= ZSURA*DSEA + (1.D0-ZSURA)*USEA                                              
        XSH=STR
        XCH=CHM
        XBH=BOT
        XTH=TOP
        XGPRO=GL
	elseif(VALUE(1).eq.4) then
c**using nucleus SF
	if(ITI.eq.1) then
	A=APRO
	ZSURA=ZPRO/APRO
	elseif(ITI.eq.2) then
	A=ACIB 
	ZSURA=ZCIB/ACIB 
	 endif
	 CALL STRUCTA(X,QQ,A,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
         XUH= ZSURA*(UPV+USEA) + (1.D0-ZSURA)*(DNV+DSEA)                                         
         XUBH= ZSURA*USEA + (1.D0-ZSURA)*DSEA                                       
         XDH= ZSURA*(DNV+DSEA) + (1.D0-ZSURA)*(UPV+USEA)                                         
         XDBH= ZSURA*DSEA + (1.D0-ZSURA)*USEA                                              
         XSH=STR
         XCH=CHM
         XBH=BOT
         XTH=TOP
         XGPRO=GL	
	endif
      ELSE IF (IH.EQ.4.or.IH.eq.5) THEN 
c** pion                       
        XUH=UPV*DFLOAT(5-IH)+USEA                                          
        XUBH=UPV*DFLOAT(IH-4)+USEA                                           
        XDH=DNV*DFLOAT(IH-4)+DSEA                                          
        XDBH=DNV*DFLOAT(5-IH)+DSEA                                           
        XSH=STR
        XCH=CHM
        XBH=BOT
        XTH=TOP
        XGPRO=GL
      ENDIF
      if(x.ne.xsave) then
      write(8,*)'X changed in fonstru from ',xsave,'to ',x
      endif	
c      write(65,*) 'in fonstru x',x,'q2',q2,'iti',iti,'xgpro',xgpro
c	write(65,*)' x=',x,'q2=',q2,' ih=',ih,' U=',xuh,' D=',xdh
      RETURN
      END
 
 
