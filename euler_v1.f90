! euler_v1.f90
!
! <2019-11-09> CodeRoninSY
!
! Compile fortran90 code using f2py.numpy: $> f2py -c -m euler_v1 euler_v1.f90
! Use this Fortran code in "euler_v1.py" via import euler_v1
! Approximates y(t) in y'(t) = f(y, t) with y(a) = y0 and t = a..b and the
! step size h.
subroutine euler(func, y0, a, b, h)
!f2py intent(in) :: y0
!f2py intent(in) :: a
!f2py intent(in) :: b
!f2py intent(in) :: h
  integer, PARAMETER :: DP = selected_real_kind(8)
    external func
    real(kind=DP) :: y0, a, b, h
    real(kind=DP) :: func
    real(kind=DP) :: t, y
    if (a > b) return
    if (h <= 0) stop "negative step size"
    print '("# h = ", F0.3)', h
    print '("Time    ",4x,"Temperature     ")'
    print '("--------",4x,"----------------")'
    y = y0
    t = a
    do
      print '(F8.2,4x,F16.3)', t, y
      t = t + h
      if (t > b) return
      y = y + h * func(y, t)
    end do
end subroutine

! real*8 function func(temp, t) result(dTdt)
! !f2py intent(in) :: temp
! !f2py intent(hide) :: t
!     real*8, intent(in) :: temp
!     real*8 :: t

!     dTdt = -0.07 * (temp - 20)
! end