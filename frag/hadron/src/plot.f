	implicit real*8 (a-h,l-z)
	OPEN(UNIT=9,
     #	FILE='test_ipi3_10000.dat',
     #	STATUS='unknown')
	qp2 = 10000.d0
	ixmax = 21
	ipi = 3
	do i=1,ixmax
	  xc = 0.1d0+dfloat(i-1)*0.04d0
	  call fonfra(xc,ipi,qp2,xdup,xdubp,xddp,xddbp,xdsp,xdcp
     # ,xdbp,xdbbp,xdgp)
	  write (9,*) xc,xdup,xddp,xdsp,xdcp,xdbp,xdgp
	enddo
	end
