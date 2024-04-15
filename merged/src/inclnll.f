CDECK  ID>, PHONLL.
C$INIT on
	PROGRAM INCLNLL
	implicit real*8(a-h,o-z)
	COMMON/PTVALUE/PT
	COMMON/VALU/JMAR,IPT
	COMMON/EVOPT/PTEVO(19) 
	include "projec.f"
	common/choics/facin
	common/constants/cn 
	common/phoord/lhior
	COMMON/XSECT/ISIGM,ILOOP,IHIOR
	common/intchoix/lintegy
	common/yrange/ymin,ymax
	common/yvalue/rapidite
	common/output_units/io,ier,itest
        common/pion/ih3
	logical lintegy
	save
C read parameters for both programs
	call param
C initialisations for program one
 	call initone
c write units used for cross section
	write(io,*)'  '
	write(io,*)' Units for cross section (hbar*c)^2=',cn
	write(io,*)' pb/GeV^2 when (hbar*c)^2 =0.38935D+09 '	
	write(io,*)' mb/GeV^2 when (hbar*c)^2 =0.38935D0 '	
	IF(ISIGM.eq.1)then
	   write(io,*) '    cross section   DSIG/DY/D2PT' 
	 ELSEIF(ISIGM.eq.2) then
	   write(io,*) '    cross section    E DSIG/D3P'
	 ELSEIF(ISIGM.eq.3) then
           write(io,*) '    cross section   DSIG/DY/DPT'
	 ENDIF
	 if(lintegy) then
	 write(io,*)' averaged in rapidity from ',ymin,' to ',ymax
	 else
	 write(io,*)' at rapidity=',rapidite
	 endif 
	 if (ih3.eq.0) then
	   write(io,*) ' Inclusive production of photon'
	 else if (ih3.eq.1) then
	   write(io,*) ' Inclusive production of (pi^+ + pi^-)/2'
	 else if (ih3.eq.2) then
	   write(io,*) ' Inclusive production of (K^+ + K^-)/2'
	 else if (ih3.eq.3) then
	   write(io,*) ' Inclusive production of (K^0 + K^0_bar)/2'
	 else if (ih3.eq.4) then
	   write(io,*) ' bad choice for produce particle (ih3)'
	 else if (ih3.eq.5) then
	   write(io,*) ' Inclusive production of pi^0'
	 else if (ih3.eq.6) then
	   write(io,*) ' bad choice for produce particle (ih3)'
	 else if (ih3.eq.6) then
	   write(io,*) ' Inclusive production of h^+ + h^-'
	 endif
c initialisation for program two (done in param + initone)
c loop on pt values
        do i=1,IPT
c call program one ABFS 
C AUTHORS: P.AURENCHE, R.BAIER, M.FONTANNAZ, D.SCHIFF
	PT=PTEVO(I)
	if (lhior.eq.4) then
C call program two fragmentation
C AUTHORS: AVERSA, P.CHIAPPETTA, M.GRECO, J.PH.GUILLET
 
 	  call pioincl(csbor,csect)
	  write(io,*) ' '
	  write(io,*) '  ****** Pt ',pt,' GeV/c  SIGTOT NLL ',csect
	  write(io,*)'    where  BORN =',csbor,
     1    ' and HIGH ORDER =',csect-csbor
	  fk2=csect/csbor
	  write(io,*)' K factor  fragmentation  NLL/BORN ',fk2 
	else
	  ipt1=i
	  ipt2=i
	  khior=lhior
	  if(lhior.eq.3) khior=2       
	  CALL abfs(khior,totint,anoint,borint,higint)
c	write(io,*)' totint ',totint,' anoint ',anoint
c	write(io,*)' BORINT ',BORINT,' HIGINT ',higint
          if(lhior.eq.0) then
	    write(io,*) ' '
	    write(io,*) '  ****** Pt ',pt,' GeV/c  BORN =',borint*facin
	  elseif(lhior.eq.1) then
	    sum=facin*(borint+anoint)
	    write(io,*) ' ' 
	    write(io,*) '  ****** Pt ',pt,' GeV/c  SIGTOT=',sum
	    write(io,*)'  including BORN =',borint*facin,
     1'     and ANOMA =',anoint*facin
          elseif(lhior.eq.2) then
	    sum=totint*facin
	    write(io,*) '  ****** Pt ',pt,' GeV/c  SIGTOT NLL ',sum
	    write(io,*)' including   BORN =',borint*facin,
     1'     ANOMA =',anoint*facin,	 
     1'     and HIGH ORDER =',higint*facin
            fk1=sum/(borint*facin)
	    write(io,*)' K factor ABFS NLL/BORN ',fk1
          elseif(lhior.eq.3) then
	 
C call program two fragmentation
C AUTHORS: AVERSA, P.CHIAPPETTA, M.GRECO, J.PH.GUILLET
 
 	    call pioincl(csbor,csect)
c sum of the two programs
	    diff=anoint*facin-csbor
	    tolerance=1.D-2*csbor
	    if(dabs(diff).gt.tolerance) then
	       write(io,*)' two programs anoma differ by',diff
	       write(io,*)' check number of integration points'
	    endif
	
	     phosum=(totint-anoint)*facin 
	     sum=phosum+csect
	     write(io,*) ' '
	     write(io,*) '  ****** Pt ',pt,' GeV/c  SIGTOT NLL ',sum
	     write(io,*)' including  NLL ABFS',phosum
	     write(io,*)'    where  BORN =',borint*facin,
     1       ' and HIGH ORDER =',higint*facin
	     write(io,*)' and NLL fragmentation',csect
	     write(io,*)'    where BREMSS born =',csbor,
     1       ' and BREMSS high order =',csect-csbor 
	     fk1=phosum/(borint*facin)
	     fk2=csect/csbor
	     fk3=sum/(borint*facin+csbor)
	     write(io,*)' K factor  ABFS  NLL/BORN ',fk1
	     write(io,*)' K factor  fragmentation  NLL/BORN ',fk2 
	     write(io,*)' K factor  NLL/(BORN+BREMSS) ',fk3
	  endif
	endif
	enddo	
	call PDFSTA
      STOP
      END
