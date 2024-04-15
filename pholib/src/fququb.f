CDECK  ID>, FQUQUB.
      DOUBLE PRECISION FUNCTION FQUQUB(V)
      implicit real*8(a-h,o-z)
	include "qcd.f"
  	include "vari.f"
      include "phoscale.f"
      include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      TEMP = 0.D0
      FQUQUB = 0.D0
      DO II=1,6
       IF (II.LE.3) IIA=II+3
       IF (II.GT.3) IIA=II-3
       TEMP = Q(II)**2*QQ1(II)*QQ2(IIA)  + TEMP
      ENDDO
      FQUQUB = ALFAS*TEMP*(V**2+V1**2)
      RETURN
      END
