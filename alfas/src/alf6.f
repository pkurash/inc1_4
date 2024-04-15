CDECK  ID>, ALF6.
      DOUBLE PRECISION  FUNCTION ALF6(SCALE)
      IMPLICIT REAL*8 (A-H,L-Z)
      COMMON/WW/W03,W04,W05,W06,W13,W14,W15,W16
      COMMON/ALFA/MASCH2,MASBO2,MASTO2,LAMBDA2
      EXTERNAL ALF5
      XLOG=DLOG(SCALE/LAMBDA2)
      XLOGB=DLOG(MASTO2/LAMBDA2)
      ALF6B=-1.D0/ALTAS(6,MASTO2)+ALF5(MASTO2)
      ALF6=1.D0/ALTAS(6,SCALE)+ALF6B
      RETURN
      END