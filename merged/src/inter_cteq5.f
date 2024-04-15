c subroutine pdfset pour cteq5
	subroutine pdfset(parm,value)
	implicit real*8 (a-h,l-z)
        character*20 parm(20)
	dimension value(20)
        common/w50513/xmin,xmax,q2min,q2max
	common/w50512/qcdl4,qcdl5
	common/mrss/imode
	logical lo1,lo2
	ival1 = int(value(1)+1.d-3)
	ival2 = int(value(2)+1.d-3)
	ival3 = int(value(3)+1.d-3)
	lo1 = ival1.eq.1.and.ival2.eq.4
	if(lo1) ival3=ival3-47
	imode=ival3
	lo2 = ival3.ge.1.and.ival3.le.7
c	write(6,*)' ival1 ',ival1,' ival2 ',ival2,' ival3 ',ival3
c	write(6,*)' in pdfset imode ',imode
	if (lo1.and.lo2) then
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C ---------------------------------------------------------------------------
	  xmin = 1.d-5
	  xmax = 1.d0
	  q2min = 1.d0
	  q2max = 1.d8
          imode = ival3
	  if (imode.eq.1.or.imode.eq.2.or.imode.eq.4.or.imode.eq.5) then
            qcdl4 = (.326)
            qcdl5 = (.226)
	  else if (imode.eq.3) then
            qcdl4 = (.192)
            qcdl5 = (.146)
	  else if (imode.eq.6) then
            write (*,*) 'the program cannot run with 3 flavours'
            write (*,*) 'choose a set with at least four flavours'
	    stop
	  else if (imode.eq.7) then
            qcdl4 = (.309)
c attention c'est completement faux
	    qcdl5 = qcdl4
	  endif
c
	else
	  write (8,*) 'bad choice for value',value
	endif
	return
	end
c subroutine pftopdg pour cteq5
	subroutine pftopdg(x,m,dxpdf)
	implicit real*8 (a-h,l-z)
	dimension dxpdf(-6:6)
	common/mrss/imode
	if (m*m.le.1.d+8) then
	  q = m
	else
	  q = dsqrt(1.d+8)
	endif
	call SetCtq5(imode)
c
	dxpdf(-6) = 0.d0
	dxpdf(-5) = x*Ctq5Pdf(-5,x,q)
	dxpdf(-4) = x*Ctq5Pdf(-4,x,q)
	dxpdf(-3) = x*Ctq5Pdf(-3,x,q)
	dxpdf(-2) = x*Ctq5Pdf(-1,x,q)
	dxpdf(-1) = x*Ctq5Pdf(-2,x,q)
	dxpdf(0) = x*Ctq5Pdf(0,x,q)
	dxpdf(1) = x*Ctq5Pdf(2,x,q)-dxpdf(-1)
	dxpdf(2) = x*Ctq5Pdf(1,x,q)-dxpdf(-2)
	dxpdf(3) = x*Ctq5Pdf(3,x,q)
	dxpdf(4) = x*Ctq5Pdf(4,x,q)
	dxpdf(5) = x*Ctq5Pdf(5,x,q)
	dxpdf(6) = 0.d0	
	return
	end
c subroutine structm pour cteq5
	SUBROUTINE STRUCTM(X,M,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
	implicit real*8 (a-h,l-z)
c	dimension dxpdf(-6:6)
	common/mrss/imode
	if (m*m.le.1.d+8) then
	  q = m
	else
	  q = dsqrt(1.d+8)
	endif
	call SetCtq5(imode)
c
c	dxpdf(-6) = 0.d0
c	dxpdf(-5) = x*Ctq5Pdf(-5,x,q)
c	dxpdf(-4) = x*Ctq5Pdf(-4,x,q)
c	dxpdf(-3) = x*Ctq5Pdf(-3,x,q)
	USEA = x*Ctq5Pdf(-1,x,q)
	DSEA = x*Ctq5Pdf(-2,x,q)
	GL = x*Ctq5Pdf(0,x,q)
	DNV = x*Ctq5Pdf(2,x,q)-DSEA
	UPV = x*Ctq5Pdf(1,x,q)-USEA
	STR = x*Ctq5Pdf(3,x,q)
	CHM = x*Ctq5Pdf(4,x,q)
	BOT = x*Ctq5Pdf(5,x,q)
	TOP = 0.d0	
	return
	end
c dummy	
	SUBROUTINE STRUCTA(X,M,A,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)	
	implicit real*8 (a-h,l-z)
	return
	end
c dummy
	SUBROUTINE PDFSTA
	return
	end
