CDECK  ID>, ALTAS.
      DOUBLE PRECISION FUNCTION ALTAS(INF,Q2)
      IMPLICIT REAL*8(A-H,L-Z)
      COMMON/ALFA/MASCH2,MASBO2,MASTO2,LAMBDA2
      COMMON/OPTIMGUI/AB,AC,ALQ2S
      DATA PI/3.141592653589793D0/
      EXTERNAL STEVENSON
      NF=DFLOAT(INF)
      N=3.D0
      B1=102.D0-38.D0*NF/3.D0
      B0=(11.D0*N-2.D0*NF)/3.D0
      B0INV=4.D0*PI/B0
      AB=0.5D0*B0
      AC=.25D0*B1/B0
      ALQ2S=DLOG(Q2/LAMBDA2)
      CALL DZERO(.001D0,.25D0,AA,ZERO,.005D0,100,STEVENSON)
      ALTAS=AA*PI/2.D0/PI
      RETURN
      END
