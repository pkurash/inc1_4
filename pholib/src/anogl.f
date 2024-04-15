CDECK  ID>, ANOGL.
      DOUBLE PRECISION FUNCTION ANOGL(X)
      implicit real*8(a-h,o-z)
      COMMON/CHARGE/Q(8)
      COMMON/ANOMAL/IANOM,IANORD
      COMMON/X2ANOM/X2FRAG(8),X2G,DLQ,DLQCH
      save
      ANOGL=0.D0
      IF (IANOM.EQ.1) THEN
          ANOGL=0.0243*(1.-X) / X**1.97
          ELSE IF(IANOM.EQ.2) THEN
          ANOGL=X2G/X/DLQ
      ENDIF
      RETURN
      END
