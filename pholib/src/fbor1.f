CDECK  ID>, FBOR1.
      DOUBLE PRECISION FUNCTION FBOR1(V,W)
      implicit real*8(a-h,o-z)
	include "qcd.f"
	include "vari.f"
        include "phoscale.f"
       include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      TEMP = 0.D0
      FBOR1 = 0.D0
      DO  II=1,3
       TEMP = TEMP + ANOM(X2,II)
      ENDDO
      FBOR1 = ALFAS*ALFPI * LLQ * V**2*V1*W/X2 * GG1*GG2 *
     1 ( TEMP*( (V1/VW+VW/V1)/NC-(V1**2+VW**2)/(X2**2*CF) )
     2  +ANOGL(X2)*(9.D0/2.D0)*(3.D0-VW*V1*X2**(-2)
     3  +X2*(V*W/V1**2+V1/(V*W)**2)) )
      RETURN
      END
