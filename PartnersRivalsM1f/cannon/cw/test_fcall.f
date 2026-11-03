*     test_fcall.f -- m1cls called from fixed-form Fortran as an implicit external, the way the
*     simulator calls it.  gfortran and ifx (Linux) both mangle this to the symbol m1cls_.
*        gfortran -O2 -o test_fcall.x test_fcall.f m1cls.o
*     Prints four atoms / strata; test_cls.py checks them against the C driver.
      program tfcall
      implicit real*8 (a-h,o-z)
      dimension s(4)
      integer iat,ist,inr
      eps=1.d-04
      om=1.d+00-eps
*     ALLD at the PD (u,v) = (-0.5,0.5): 011 = 3, stratum 0
      s(1)=eps
      s(2)=eps
      s(3)=eps
      s(4)=eps
      u=-0.5d+00
      v=0.5d+00
      call m1cls(s,eps,u,v,iat,ist,inr)
      write (*,'(3i4)') iat,ist,inr
*     WSLS at (-0.1,0.1): 110 = 6, stratum 28
      s(1)=om
      s(4)=om
      u=-0.1d+00
      v=0.1d+00
      call m1cls(s,eps,u,v,iat,ist,inr)
      write (*,'(3i4)') iat,ist,inr
*     an interior point at a wedge-S game (u,v) = (1.5,-0.25)
      s(1)=om
      s(2)=0.7d+00
      s(3)=0.6d+00
      s(4)=0.2d+00
      u=1.5d+00
      v=-0.25d+00
      call m1cls(s,eps,u,v,iat,ist,inr)
      write (*,'(3i4)') iat,ist,inr
*     GTFT-like (1, x, 1-x, 0) at the PD
      s(1)=om
      s(2)=0.75d+00
      s(3)=0.25d+00
      s(4)=eps
      u=-0.5d+00
      v=0.5d+00
      call m1cls(s,eps,u,v,iat,ist,inr)
      write (*,'(3i4)') iat,ist,inr
      end
