CDECK  ID>, FBCH1.
      DOUBLE PRECISION FUNCTION FBCH1(V,W)
      implicit real*8(a-h,o-z)
	include "qcd.f"
	include "vari.f"
       include "phoscale.f"
      include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      FBCH1 = 0.D0
      TEMP= ANOM(X2,7)
      FBCH1 = ALFAS*ALFPI* LLQCH * V**2*V1*W/X2*
     +TEMP *GG1*GG2 *
     1 ( (V1/VW+VW/V1)/NC-(V1**2+VW**2)/(X2**2*CF) )
      RETURN
      END
