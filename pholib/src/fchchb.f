CDECK  ID>, FCHCHB.
      DOUBLE PRECISION FUNCTION FCHCHB(V)
      implicit real*8(a-h,o-z)
      include "qcd.f"
	include "vari.f"
      include "phoscale.f"
      include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      FCHCHB = 0.D0
      TEMP = Q(7)**2 * (QQ1(7)*QQ2(8)+QQ1(8)*QQ2(7))
c      IF(NF.eq.5) TEMP=TEMP+ Q(9)**2 * (QQ1(9)*QQ2(10)+QQ1(10)*QQ2(9))
      FCHCHB = ALFAS*TEMP*(V**2+V1**2)
      RETURN
      END
