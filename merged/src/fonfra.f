      subroutine fonfra(x,ih3,qp2,xdup,xdubp,xddp,xddbp,xdsp,xdcp
     # ,xdbp,xdbbp,xdtp,xdtbp,xdgp)
      implicit real*8 (a-h,l-z)
      dimension dh(0:10)
      COMMON/ANOMAL/IANOM,IANORD
      save
      qs = dsqrt(qp2)
      qp2save=qp2
      if (ih3.eq.0) then
        CALL distributionPert(x,QP2,'UP',up_Pert)
        CALL distributionPert(x,QP2,'DO',down_Pert)
        CALL distributionPert(x,QP2,'SB',strange_Pert)
        CALL distributionPert(x,QP2,'CB',charm_Pert)
        CALL distributionPert(x,QP2,'BB',bottom_Pert)
        CALL distributionPert(x,QP2,'GL',gluon_Pert)
	
	if(ianord .eq. 1) then
 	     CALL distributionNonPert_setI(x,QP2,'UP',up_VDM)
 	     CALL distributionNonPert_setI(x,QP2,'DO',down_VDM)
 	     CALL distributionNonPert_setI(x,QP2,'SB',strange_VDM)
 	     CALL distributionNonPert_setI(x,QP2,'CB',charm_VDM)
 	     CALL distributionNonPert_setI(x,QP2,'BB',bottom_VDM)
 	     CALL distributionNonPert_setI(x,QP2,'GL',gluon_VDM)
        elseif(ianord .eq. 2) then
 	     CALL distributionNonPert_setII(x,QP2,'UP',up_VDM)
 	     CALL distributionNonPert_setII(x,QP2,'DO',down_VDM)
 	     CALL distributionNonPert_setII(x,QP2,'SB',strange_VDM)
 	     CALL distributionNonPert_setII(x,QP2,'CB',charm_VDM)
 	     CALL distributionNonPert_setII(x,QP2,'BB',bottom_VDM)
 	     CALL distributionNonPert_setII(x,QP2,'GL',gluon_VDM)
        else
 	     write(*,*) 'unknown set'
 	     stop
        endif
	if(qp2save.ne.qp2) then
	 write(8,*)'!!!scale changed in bfg!!!from',qp2save,'to',qp2
	endif
	up = up_Pert + up_VDM
	down = down_Pert + down_VDM
	strange = strange_Pert + strange_VDM
	charm = charm_Pert + charm_VDM
	bottom = bottom_Pert + bottom_VDM
	gluon = gluon_Pert + gluon_VDM

        xup=up*x
        xdown=down*x
        xstrange=strange*x
        xcharm=charm*x
        xbottom=bottom*x
        xglue=gluon*x
c rajoute pour interface avec hadlib.f
        xdup = xup
        xdubp = xup
        xddp = xdown
        xddbp = xdown
        xdsp = xstrange
        xdcp = xcharm
        xdbp = xbottom
        xdbbp = xbottom
        xdtp = 0.d0
        xdtbp = 0.d0
        xdgp = xglue
      else
c iset =1 kkp fragmentation function at NTL
	iset = 1
	if(ianord.eq.1) then
	  qsave=qs
          call bkk(ih3,iset,x,qs,dh)
	  if(qs.ne.qsave) then
	   write(8,*)' !!! scale changed in bkk from',qsave,'to',qs
	  endif       
	xdup=dh(1)*x
        xdubp=dh(2)*x
        xddp=dh(3)*x
        xddbp=dh(4)*x
        xdsp=dh(5)*x
        xdcp=dh(7)*x
        xdbp=dh(9)*x
        xdbbp=dh(10)*x
        xdtp=0.d0
        xdtbp=0.d0
        xdgp=dh(0)*x
	elseif(ianord.eq.2) then 
	  qsave=qs 
          call kkp(ih3,iset,x,qs,dh)	
        xdup=dh(1)*x
        xdubp=dh(2)*x
        xddp=dh(3)*x
        xddbp=dh(4)*x
        xdsp=dh(5)*x
        xdcp=dh(7)*x
        xdbp=dh(9)*x
        xdbbp=dh(10)*x
        xdtp=0.d0
        xdtbp=0.d0
        xdgp=dh(0)*x
	  if(qs.ne.qsave) then
	   write(8,*)'!!!scale changed in kkp from',qsave,'to',qs
	  endif
	elseif(ianord.eq.3) then 
	 ipi=1
	   call fonfrac(x,ipi,qp2,xdup,xdubp,xddp,xddbp,xdsp,xdcp
     # ,xdbp,xdbbp,xdgp)
	else
	  write(*,*) 'set hadron-fragmentation not implemented'
	  stop 
	endif 

      endif
      return
      end
