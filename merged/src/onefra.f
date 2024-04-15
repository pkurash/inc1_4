	subroutine onefra
c fragmentation function initialization
     	implicit real*8(a-h,o-z)
	COMMON/OUTPUT_UNITS/IO,IER,ITEST
     	COMMON/ANOMAL/IANOM,IANORD
        COMMON/PION/IH3

c ----- choice of anomalous photon structure function
c ----- ianom=0 no qcd corrections ; ianom=1 leading log qcd approx
c ----- ianom=2, ianord=1 ll tables from chiappetta and guillet et al.
c ----- ianom=2, ianord=2 ntl tables from chiappetta and guillet et al.
c
c      standard choice IANOM= 2 and   IANORD= 2
	save
      if (ih3.eq.0) then      
      IF(IANOM.EQ.0) WRITE(io,*) ' LOWEST ORDER ANOMAL EXPRESSION'
      IF(IANOM.EQ.1) WRITE(io,*) ' LEAD LOG ANOMAL EXPR. WITH    GLUON'
      IF(IANOM.EQ.2) THEN
      WRITE(io,*) ' ANOMALOUS PHOTON FRAGMENTATION FROM BOURHIS-GUILLET'
       IF(IANORD.EQ.1) WRITE(io,*) ' SET I, SMALL GLUON'
       IF(IANORD.EQ.2) WRITE(io,*) ' SET II, LARGE GLUON'
      ENDIF 
      else if (ih3.eq.1) then 
        if(ianord.eq.1) then    
        WRITE(io,*) '(PI+ + PI-)/2 FRAGMENTATION FROM BKK'
	elseif(ianord.eq.2) then
	 WRITE(io,*) '(PI+ + PI-)/2 FRAGMENTATION FROM KKP'
	endif 
      else if (ih3.eq.2) then 
         if(ianord.eq.1) then 
        WRITE(io,*) '(K+ + K-)/2 FRAGMENTATION FROM BKK'
	 elseif(ianord.eq.2) then
	WRITE(io,*) '(K+ + K-)/2 FRAGMENTATION FROM KKP'
	 endif
      else if (ih3.eq.3) then   
        if(ianord.eq.1) then  
        WRITE(io,*) '(K0 + K0BAR)/2 FRAGMENTATION FROM BKK'
	elseif(ianord.eq.2) then
	WRITE(io,*) '(K0 + K0BAR)/2 FRAGMENTATION FROM KKP'
	endif
      else if (ih3.eq.5) then
      if(ianord.eq.1) then     
        WRITE(io,*) 'PI0 FRAGMENTATION FROM BKK'
	elseif(ianord.eq.2) then
	WRITE(io,*) 'PI0 FRAGMENTATION FROM KKP'
	endif
      else if (ih3.eq.7) then
      if(ianord.eq.1) then     
        WRITE(io,*) 'H+ + H- FRAGMENTATION FROM BKK'
	elseif(ianord.eq.2) then
	WRITE(io,*) 'H+ + H- FRAGMENTATION FROM KKP'
	elseif(ianord.eq.3) then
	WRITE(io,*) 'H+ + H- FRAGMENTATION FROM BFGW'
	endif
      else
        WRITE(ier,*) 'IH3 has not a good value for fragmentation'
	stop
      endif 
      return
      end
