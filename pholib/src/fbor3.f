CDECK  ID>, FBOR3.
      DOUBLE PRECISION FUNCTION FBOR3(V,W)
      implicit real*8(a-h,o-z)
	include "qcd.f"
	include "vari.f"
        include "phoscale.f"
       include "kine.f"
      COMMON/CHARGE/Q(8)
      save
      TEMPE = 0.D0
      TEMPF = 0.D0
      TEMPJ = 0.D0
      TEMPK = 0.D0
      FBOR3 = 0.D0
      ANOG=ANOGL(X2)
      DO 10 II=1,6
C
      IF (II.LE.3) THEN
       IIA = II+3
       IK  = II+1
       IF(IK.GT.3)  IK=IK-3
       IKA = IK+3
       IL  = II+2
       IF(IL.GT.3)  IL=IL-3
       ILA = IL+3
      ELSEIF(II.GT.3) THEN
       IIA = II-3
       IK  = II+1
       IF (IK.GT.6) IK=IK-3
       IKA = IK-3
       IL=II+2
       IF (IL.GT.6) IL=IL-3
       ILA = IL-3
      ENDIF
C
      ANOI=ANOM(X2,II)
      ANOK=ANOM(X2,IK)
      ANOL=ANOM(X2,IL)
      TEMPE = TEMPE + QQ1(II)*QQ2(IIA) *
     1 ( 2.*ANOK+2.*ANOL )*(V1**2+VW**2)/X2**2
      TEMPF = TEMPF + QQ1(II)*( (QQ2(IK)+QQ2(IKA)) *
     1(ANOI*(X2**2+VW**2)/V1**2+ANOK*(X2**2+V1**2)/VW**2)
     1                            +(QQ2(IL)+QQ2(ILA)) *
     1(ANOI*(X2**2+VW**2)/V1**2+ANOL*(X2**2+V1**2)/VW**2) )
      TEMPJ = TEMPJ + QQ1(II)*QQ2(IIA) *
     1        ( ANOI *
     1 (2.*(V1**2+VW**2)/X2**2+(X2**2+VW**2)/V1**2+(X2**2+V1**2)/VW**2
     1  +(2./NC)*( VW**2/(X2*V1)+V1**2/(X2*VW)))
     1          +ANOG * 2. *
     1  (2.66667*(VW/V1+V1/VW) - 6.*(V1**2+VW**2)/X2**2) )
      TEMPK = TEMPK + QQ1(II)*QQ2(II) * ANOI *
     1 ((X2**2+V1**2)/VW**2+(X2**2+VW**2)/V1**2-(2./NC)*X2**2/(V1*VW))
c
c      TEMPE = TEMPE + QQ1(II)*QQ2(IIA) *
c     1 (Q(IK)**2+Q(IKA)**2+Q(IL)**2+Q(ILA)**2)*(V1**2+VW**2)/X2**2
c      TEMPF = TEMPF + QQ1(II)*( (QQ2(IK)+QQ2(IKA)) *
c     1(Q(II)**2*(X2**2+VW**2)/V1**2+Q(IK)**2*(X2**2+V1**2)/VW**2)
c     1                            +(QQ2(IL)+QQ2(ILA)) *
c     1(Q(II)**2*(X2**2+VW**2)/V1**2+Q(IL)**2*(X2**2+V1**2)/VW**2) )
c      TEMPJ = TEMPJ + QQ1(II)*QQ2(IIA) * Q(II)**2 *
c     1 (2.*(V1**2+VW**2)/X2**2+(X2**2+VW**2)/V1**2+(X2**2+V1**2)/VW**2
c     1  +(2./NC)*( VW**2/(X2*V1)+V1**2/(X2*VW)))
c      TEMPK = TEMPK + QQ1(II)*QQ2(II) * Q(II)**2 *
c     1 ((X2**2+V1**2)/VW**2+(X2**2+VW**2)/V1**2-(2./NC)*X2**2/(V1*VW))
c
10    CONTINUE
      ANORM = ALFAS*ALFPI * LLQ * (V**2*V1*W/X2)
c      ANORM = ALFAS*ALFPI*ANOM(X2)*LLQ * (V**2*V1*W/X2)
      FBOR3 = ANORM * (TEMPE+TEMPF+TEMPJ+TEMPK)
      RETURN
      END
