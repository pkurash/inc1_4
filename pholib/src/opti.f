CDECK  ID>, OPTI.
      DOUBLE PRECISION FUNCTION OPTI(AA)
      implicit real*8(a-h,o-z)
      COMMON/OPTIM/AB,AC,BBORN,ANOMA,HIGOR
	include "dreyan.f"
	save
      CONST=AA*AC/(1.D0+AC*AA)
      ALG=LOG(0.5D0*AB*CONST/AC)
      ANOMM=ANOMA
C     IF(IDY.EQ.1) ANOMM=0.
      OPTI=2.D0*((1.+AA*AC*ALG+.5D0*CONST)*BBORN+AA*(ANOMM+HIGOR))
      RETURN
      END
