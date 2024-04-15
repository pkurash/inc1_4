CDECK  ID>, FBCH2.
      DOUBLE PRECISION FUNCTION FBCH2(V,W)
      implicit real*8(a-h,o-z)
	include "qcd.f"
	include "vari.f"
        include "phoscale.f"
      include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      FBCH2 = 0.D0
      ANOI=ANOM(X2,7)
      ANOG=ANOGL(X2)
      TEMP=   2.D0*ANOI*QQ2(7)
      TEMP1=  2.D0*ANOI*QQ1(7)
      TEMGL=  2.D0*ANOG*QQ2(7)
      TEMGL1= 2.D0*ANOG*QQ1(7)
      FBCH2 =ALFAS*ALFPI * LLQCH * (V**2*V1*W/X2) *
     1      ( (GG1*TEMP+GG2*TEMGL1) *
     2        ((V1/X2+X2/V1)*CF/NC+(V1**2+X2**2)/VW**2) +
     3        (GG2*TEMP1+GG1*TEMGL) *
     4        ((VW/X2+X2/VW)*CF/NC+(X2**2+VW**2)/V1**2) )
      RETURN
      END
