! test_fbind.f90 -- m1cls called through an iso_c_binding interface block (symbol `m1cls`),
! the portable route for any Fortran compiler.  Same four cases as test_fcall.f.
!    gfortran -O2 -o test_fbind.x test_fbind.f90 m1cls.o
program tfbind
  use iso_c_binding, only: c_double, c_int
  implicit none
  interface
     subroutine m1cls(s, eps, u, v, atom, stratum, neartie) bind(C, name='m1cls')
       import :: c_double, c_int
       real(c_double), intent(in) :: s(4), eps, u, v
       integer(c_int), intent(out) :: atom, stratum, neartie
     end subroutine m1cls
  end interface
  real(c_double) :: s(4), eps, om, u, v
  integer(c_int) :: iat, ist, inr
  eps = 1.0d-4
  om = 1.0d0 - eps
  s = eps; u = -0.5d0; v = 0.5d0
  call m1cls(s, eps, u, v, iat, ist, inr); write (*, '(3i4)') iat, ist, inr
  s(1) = om; s(4) = om; u = -0.1d0; v = 0.1d0
  call m1cls(s, eps, u, v, iat, ist, inr); write (*, '(3i4)') iat, ist, inr
  s = [om, 0.7d0, 0.6d0, 0.2d0]; u = 1.5d0; v = -0.25d0
  call m1cls(s, eps, u, v, iat, ist, inr); write (*, '(3i4)') iat, ist, inr
  s = [om, 0.75d0, 0.25d0, eps]; u = -0.5d0; v = 0.5d0
  call m1cls(s, eps, u, v, iat, ist, inr); write (*, '(3i4)') iat, ist, inr
end program tfbind
