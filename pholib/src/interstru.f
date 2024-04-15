	subroutine interstru
C Interface between Aurenche et al and Guillet et al SF	
	implicit real*8(a-h,o-z)     
c	 REAL*8 NC
	include "vari.f"	
     	include "phoscale.f"  
	include "projec.f"
	save
		
		if(IPRO.eq.1) ih=0
		if(IPRO.eq.2) ih=1
		if(IPRO.eq.3) ih=4
		if(IPRO.eq.4) ih=5
		if(IPRO.eq.6) ih=3
c			
c	write(65,*) 'in interstru before fonstru projectile ih',ih	
c	projectile
	ITI=1
      call fonstru(ITI,T1,ih,Q1S,xuh,xubh,xdh,xdbh,xsh,xch,
     * xbh,xth,xgpro)
c       write(65,*) 'in interstru after fonstru ih',ih
      GG1=xgpro
      QQ1(1)=xuh
      QQ1(2)=xdh
      QQ1(3)=xsh
      QQ1(4)=xubh
      QQ1(5)=xdbh
      QQ1(6)=xsh
      QQ1(7)=xch
      QQ1(8)=qq1(7)
c     qq1(9)=xbh
c     qq1(10)=qq1(9)
	ITI=2
	ih=itar
c	write(65,*) 'in interstru before fonstru target ih',ih	
      call fonstru(ITI,T2,ih,Q1S,xuh,xubh,xdh,xdbh,xsh,xch,
     * xbh,xth,xgpro)
c      write(65,*) 'in interstru before fonstru ih',ih
c	target
      GG2=xgpro
      QQ2(1)=xuh
      QQ2(2)=xdh
      QQ2(3)=xsh
      QQ2(4)=xubh
      QQ2(5)=xdbh
      QQ2(6)=xsh
      QQ2(7)=xch
      QQ2(8)=qq2(7)
c     qq2(9)=xbh
c     qq2(10)=qq1(9) 
      return      
	end
