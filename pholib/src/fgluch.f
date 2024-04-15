CDECK  ID>, FGLUCH.
      DOUBLE PRECISION FUNCTION FGLUCH(V)
      implicit real*8(a-h,o-z)
	include "qcd.f"
	include "vari.f"
      include "phoscale.f"
      include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      TEMP = 0.0D0
      TEMP1 = 0.0D0
      FGLUCH = 0.0D0
      IMI=7
      IMA=8
c      if(nf.eq.5) ima=10
      DO  II=imi,ima
       TEMP = TEMP + Q(II)**2*QQ2(II)*(V+1.D0/V)
       TEMP1= TEMP1+Q(II)**2*QQ1(II)*(V1+1.D0/V1)
      ENDDO
      FGLUCH = ALFAS*V1*V * (GG1*TEMP+GG2*TEMP1)
      RETURN
      END
