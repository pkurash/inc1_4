CDECK  ID>, FBCH3.
      DOUBLE PRECISION FUNCTION FBCH3(V,W)
      implicit real*8(a-h,o-z)
	include "qcd.f"
 	include "vari.f"
        include "phoscale.f"
       include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      TEMPE = 0.D0
      TEMPF = 0.D0
      FBCH3 = 0.D0
      ANO7=ANOM(X2,7)
      DO  II=1,6
      IF (II.LT.3) IIA = II+3
      IF (II.GT.3) IIA = II-3
      ANOI=ANOM(X2,II)
      TEMPE = TEMPE + QQ1(II)*QQ2(IIA) *
     1 ( 2.*ANO7 )*(V1**2+VW**2)/X2**2
      TEMPF = TEMPF + QQ1(II)*( (QQ2(7)+QQ2(8)) *
     1(ANOI*(X2**2+VW**2)/V1**2+ANO7*(X2**2+V1**2)/VW**2) )
     1          + 2.* QQ1(7)*( QQ2(II) *
     1(ANO7 *(X2**2+VW**2)/V1**2+ANOI*(X2**2+V1**2)/VW**2) )
      ENDDO
      TEMPJ = 2. *  QQ1(7)*QQ2(8) * ( ANOM(X2,7) *
     1 (2.*(V1**2+VW**2)/X2**2+(X2**2+VW**2)/V1**2+(X2**2+V1**2)/VW**2
     1  +(2./NC)*( VW**2/(X2*V1)+V1**2/(X2*VW)))
     1          +ANOGL(X2) * 2. *
     1  (2.66667*(VW/V1+V1/VW) - 6.*(V1**2+VW**2)/X2**2) )
      TEMPK = 2. *  QQ1(7)*QQ2(7) * ANOM(X2,7) *
     1 ((X2**2+V1**2)/VW**2+(X2**2+VW**2)/V1**2-(2./NC)*X2**2/(V1*VW))
      ANORM = ALFAS*ALFPI * LLQCH * (V**2*V1*W/X2)
      FBCH3 = ANORM * (TEMPE+TEMPF+TEMPJ+TEMPK)
      RETURN
      END
