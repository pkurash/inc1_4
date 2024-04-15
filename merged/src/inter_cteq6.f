
c subroutine pdfset pour cteq6
	subroutine pdfset(parm,value)
	implicit real*8 (a-h,l-z)
        character*20 parm(20)
	dimension value(20)
        common/w50513/xmin,xmax,q2min,q2max
	common/w50512/qcdl4,qcdl5
	common/mrss/imode
	logical lo1,lo2,lot
	ival1 = int(value(1)+1.d-3)
	ival2 = int(value(2)+1.d-3)
	ival3 = int(value(3)+1.d-3)
	lo1 = ival1.eq.1.and.ival2.eq.4
	lo2 = ival3.ge.1.and.ival3.le.3
	lot = ival3.ge.101.and.ival3.le.140
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
C     ------------------------------
C   1xx  CTEQ6M1xx  +/- w.r.t. CTEQ6M     0.118     326   226    cteq6m1xx.tbl
C    (where xx=01--40)	
	if (lo1.and.lo2) then
	  xmin = 1.d-6
	  xmax = 1.d0
	  q2min = 1.3d0**2
	  q2max = 1.d8
          imode = ival3
	  if (imode.eq.1.or.imode.eq.2.or.imode.eq.3) then
            qcdl4 = (.326)
            qcdl5 = (.226)
	  endif
c
	elseif(lo1.and.lot) then
	  xmin = 1.d-6
	  xmax = 1.d0
	  q2min = 1.3d0**2
	  q2max = 1.d8
          imode = ival3 
	  qcdl4 = (.326)
          qcdl5 = (.226)
	  
	 else 
	  write (8,*) 'bad choice for nset in cteq6',value
	endif
	return
	end
c subroutine pftopdg pour cteq6
	subroutine pftopdg(x,m,dxpdf)
	implicit real*8 (a-h,l-z)
	dimension dxpdf(-6:6)
	common/mrss/imode
	if (m*m.le.1.d+8) then
	  q = m
	else
	  q = dsqrt(1.d+8)
	endif
	call SetCtq6(imode)
c
	dxpdf(-6) = 0.d0
	dxpdf(-5) = x*Ctq6Pdf(-5,x,q)
	dxpdf(-4) = x*Ctq6Pdf(-4,x,q)
	dxpdf(-3) = x*Ctq6Pdf(-3,x,q)
	dxpdf(-2) = x*Ctq6Pdf(-1,x,q)
	dxpdf(-1) = x*Ctq6Pdf(-2,x,q)
	dxpdf(0) = x*Ctq6Pdf(0,x,q)
	dxpdf(1) = x*Ctq6Pdf(2,x,q)-dxpdf(-1)
	dxpdf(2) = x*Ctq6Pdf(1,x,q)-dxpdf(-2)
	dxpdf(3) = x*Ctq6Pdf(3,x,q)
	dxpdf(4) = x*Ctq6Pdf(4,x,q)
	dxpdf(5) = x*Ctq6Pdf(5,x,q)
	dxpdf(6) = 0.d0	
	return
	end
c subroutine structm pour cteq6
	SUBROUTINE STRUCTM(X,M,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
	implicit real*8 (a-h,l-z)
	common/mrss/imode
	if (m*m.le.1.d+8) then
	  q = m
	else
	  q = dsqrt(1.d+8)
	endif
	call SetCtq6(imode)

	USEA = x*Ctq6Pdf(-1,x,q)
	DSEA = x*Ctq6Pdf(-2,x,q)
	GL = x*Ctq6Pdf(0,x,q)
	DNV = x*Ctq6Pdf(2,x,q)-DSEA
	UPV = x*Ctq6Pdf(1,x,q)-USEA
	STR = x*Ctq6Pdf(3,x,q)
	CHM = x*Ctq6Pdf(4,x,q)
	BOT = x*Ctq6Pdf(5,x,q)
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
