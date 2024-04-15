CDECK  ID>, ANOM.
      DOUBLE PRECISION FUNCTION ANOM(X,I)
      implicit real*8(a-h,o-z)
      COMMON/CHARGE/Q(8)
      COMMON/ANOMAL/IANOM,IANORD
      COMMON/X2ANOM/X2FRAG(8),X2G,DLQ,DLQCH
      save
      ANOM=0.D0
      IF (IANOM.EQ.0) THEN
      ANOM=Q(I)**2 * (1.D0+(1.D0-X)**2)/X
          ELSE  IF(IANOM.EQ.1) THEN
          ANOM=Q(I)**2 *
     1    (2.21D0-1.28D0*X+1.29D0*X**2)*X**(.049D0)/
     2    (1.D0-1.63D0*LOG(1.D0-X))/X
     2    +0.0020D0*(1.D0-X)**2 / X**2.54D0
                ELSE IF(IANOM.EQ.2) THEN
                IF(I.LE.6) ANOM=X2FRAG(I)/X/DLQ
                IF(I.GE.7) ANOM=X2FRAG(I)/X/DLQCH
      ENDIF
      RETURN
      END
