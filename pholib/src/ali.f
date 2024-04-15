CDECK  ID>, ALI.
      DOUBLE PRECISION FUNCTION ALI(X)
      implicit real*8(a-h,o-z)
      save
      ALI=LOG(1.-X)/X
      RETURN
      END
