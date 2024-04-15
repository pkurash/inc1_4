CDECK  ID>, FGLUQU.
      DOUBLE PRECISION FUNCTION FGLUQU(V)
      implicit real*8(a-h,o-z)
 	include "qcd.f"
	include "vari.f"
      include "phoscale.f"
      include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      TEMP = 0.D0
      TEMP1 = 0.D0
      FGLUQU = 0.D0
      DO II=1,6
       TEMP = TEMP + Q(II)**2*QQ2(II)*(V+1.D0/V)
       TEMP1= TEMP1+Q(II)**2*QQ1(II)*(V1+1.D0/V1)
      ENDDO
      FGLUQU = ALFAS*V1*V * (GG1*TEMP+GG2*TEMP1)
      RETURN
      END
