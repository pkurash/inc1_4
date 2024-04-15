CDECK  ID>, FBOR2.
      DOUBLE PRECISION FUNCTION FBOR2(V,W)
      implicit real*8(a-h,o-z)
	include "qcd.f"
	include "vari.f"
        include "phoscale.f"
       include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      TEMP = 0.D0
      TEMP1 = 0.D0
      TEMGL = 0.D0
      TEMGL1 = 0.D0
      FBOR2 = 0.D0
      ANOG=ANOGL(X2)
      DO  II=1,6
       ANOI=ANOM(X2,II)
       TEMP = TEMP + ANOI*QQ2(II)
       TEMP1=TEMP1 + ANOI*QQ1(II)
       TEMGL = TEMGL + ANOG*QQ2(II)
       TEMGL1=TEMGL1 + ANOG*QQ1(II)
      ENDDO
      FBOR2 =ALFAS*ALFPI * LLQ * (V**2*V1*W/X2) *
     1      ( (GG1*TEMP+GG2*TEMGL1) *
     2        ((V1/X2+X2/V1)*CF/NC+(V1**2+X2**2)/VW**2) +
     3        (GG2*TEMP1+GG1*TEMGL) *
     4        ((VW/X2+X2/VW)*CF/NC+(X2**2+VW**2)/V1**2) )
      RETURN
      END
