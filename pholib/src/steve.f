CDECK  ID>, STEVE.
      DOUBLE PRECISION FUNCTION STEVE(AA)
      implicit real*8(a-h,o-z)
      include "opt.f"
      COMMON/OPTIM/AB,AC,BBORN,ANOMA,HIGOR
      save
      STEVE=ALQ2S-2.D0*(1.D0/AA+AC*
     1      LOG(.5D0*AA*AB/(1.D0+AA*AC)))/AB
COM   WRITE(65,*) ' AAALLAA,'  STEVE=',STEVE,'  ALQ2S=',ALQ2S
      RETURN
      END
