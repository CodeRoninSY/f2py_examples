! euler_v1.f90
!
! <2019-11-09> CodeRoninSY
!
! Compile fortran90 code using f2py.numpy: $> f2py -c -m euler_v1 euler_v1.f90
! Use this Fortran code in "euler_v1.py" via import euler_v1
! Approximates y(t) in y'(t) = f(y, t) with y(a) = y0 and t = a..b and the
! step size h.

module ode
  implicit none
  integer, parameter :: ou = 11

contains

!! Euler
subroutine euler(y, y0, t, n)
  implicit none
  !f2py intent(callback) f0
  external f0
  real(8) :: f0
  !f2py real(8) ret, arg1, arg2
  !f2py ret = f0(arg1, arg2)

  real(8), intent(in) :: y0
  real(8), intent(out) :: y(n)
  real(8), intent(in) :: t(n)
  integer n, i
  !f2py integer intent(hide),depend(t) :: n=shape(t,0)
  real(8) h

  y(1) = y0
  h = ((t(n) - t(1)) / (n-1))
  write(ou,'(//A)') '......'
  write(ou,*) "y0: ", y0, "h: ", h

  write(ou,'(//A)') '... Euler ...'
  write(ou,'("# h = ", F12.5)') h
  write(ou,'("t       ",4x,"f(t)       ")')
  write(ou,'("--------",4x,"----------------")')

  do i = 1,n-1
    write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)
    y(i+1) = y(i) + h * f0(y(i), t(i))
  end do
  write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine euler


!! Heun
subroutine heun(y, y0, t, n)
    implicit none
    !f2py intent(callback) f1
    external f1
    real(8) :: f1
    !f2py real(8) ret, arg1, arg2
    !f2py ret = f1(arg1, arg2)

    real(8), intent(in) :: y0
    real(8), intent(out) :: y(n)
    real(8), intent(in) :: t(n)
    integer n, i
    !f2py depend(t) n

    real(8) :: h, k1, k2

    y(1) = y0
    h = ((t(n) - t(1)) / (n-1))
    write(ou,'(//A)') '......'
    write(ou,*) "y0: ", y0, "h: ", h

    write(ou,'(//A)') '... Heun ...'
    write(ou,'("# h = ", F12.5)') h
    write(ou,'("t       ",4x,"f(t)     ")')
    write(ou,'("--------",4x,"----------------")')

    do i = 1, n-1
      write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)
      k1 = h * f1(y(i), t(i))
      k2 = h * f1(y(i) + k1, t(i) + h)
      y(i+1) = y(i) + ( k1 + k2) / 2.0
    end do
    write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine heun

!! RK2A
subroutine rk2a(y, y0, t, n)
    implicit none
    !f2py intent(callback) f
    external f
    real(8) :: f
    !f2py real(8) ret, arg1, arg2
    !f2py ret = f(arg1, arg2)

    real(8), intent(in) :: y0
    real(8), intent(out) :: y(n)
    real(8), intent(in) :: t(n)
    integer n, i
    !f2py depend(t) n

    real(8) :: h, k1

    y(1) = y0
    h = ((t(n) - t(1)) / (n-1))

    write(ou,'(//A)') "... RK2A ..."
    write(ou,'("# h = ", F12.5)') h
    write(ou,'("t       ",4x,"f(t)     ")')
    write(ou,'("--------",4x,"----------------")')

    do i = 1, n-1
      write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)
      h = t(i+1) - t(i)
      k1 = h * f(y(i), t(i)) / 2.0
      y(i+1) = y(i) + h * f(y(i) + k1, t(i) + h / 2.0)
    end do
    write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine rk2a

!! RK2B
subroutine rk2b(y, y0, t, n)
    implicit none
    !f2py intent(callback) f6
    external f6
    real(8) :: f6
    !f2py real(8) ret, arg1, arg2
    !f2py ret = f6(arg1, arg2)

    real(8), intent(in) :: y0
    real(8), intent(out) :: y(n)
    real(8), intent(in) :: t(n)
    integer n, i
    !f2py depend(t) n

    real(8) :: h, k1, k2

    y(1) = y0
    h = ((t(n) - t(1)) / (n-1))

    write(ou,'(//A)') "... RK2B ..."
    write(ou,'("# h = ", F12.5)') h
    write(ou,'("t       ",4x,"f(t)     ")')
    write(ou,'("--------",4x,"----------------")')

    do i = 1, n-1
      write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)
      h = t(i+1) - t(i)
      k1 = h * f6(y(i), t(i))
      k2 = h * f6(y(i) + k1, t(i+1))
      y(i+1) = y(i) + (k1 + k2) / 2.0
    end do
    write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine rk2b

!! RK3
subroutine rk3(y, y0, t, n)
    implicit none
    !f2py intent(callback) f7
    external f7
    real(8) :: f7
    !f2py real(8) ret, arg1, arg2
    !f2py ret = f7(arg1, arg2)

    real(8), intent(in) :: y0
    real(8), intent(out) :: y(n)
    real(8), intent(in) :: t(n)
    integer n, i
    !f2py depend(t) n

    real(8) :: h, k1, k2, k3

    y(1) = y0
    h = ((t(n) - t(1)) / (n-1))

    write(ou,'(//A)') "... RK3 ..."
    write(ou,'("# h = ", F12.5)') h
    write(ou,'("t       ",4x,"f(t)     ")')
    write(ou,'("--------",4x,"----------------")')

    do i = 1, n-1
      write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)
      h = t(i+1) - t(i)
      k1 = h * f7(y(i), t(i))
      k2 = h * f7(y(i) + 0.5 * k1, t(i) + 0.5 * h)
      k3 = h * f7(y(i) - k1 + k2 * 2.0, t(i) + h)
      y(i+1) = y(i) + (k1 + 4.0 * k2 + k3) / 6.0
    end do
    write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine rk3

!! RK4_38
subroutine rk4_38(y, y0, t, n)
    implicit none
    !f2py intent(callback) f8
    external f8
    real(8) :: f8
    !f2py real(8) ret, arg1, arg2
    !f2py ret = f8(arg1, arg2)

    real(8), intent(in) :: y0
    real(8), intent(out) :: y(n)
    real(8), intent(in) :: t(n)
    integer n, i
    !f2py depend(t) n

    real(8) :: h, k1, k2, k3, k4

    y(1) = y0
    h = ((t(n) - t(1)) / (n-1))

    write(ou,'(//A)') "... RK4_3/8 ..."
    write(ou,'("# h = ", F12.5)') h
    write(ou,'("t       ",4x,"f(t)     ")')
    write(ou,'("--------",4x,"----------------")')

    do i = 1, n-1
      write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)
      k1 = h * f8(y(i), t(i))
      k2 = h * f8(y(i) + (1.0/3.0) * k1, t(i) + (1.0/3.0) * h)
      k3 = h * f8(y(i) - (1.0/3.0) * k1 + k2, t(i) + (2.0/3.0) * h)
      k4 = h * f8(y(i) + k1 - k2 + k3, t(i+1))
      y(i+1) = y(i) + (k1 + 3.0 * (k2 + k3) + k4) / 8.0
    end do
    write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine rk4_38

!! RK4
subroutine rk4(y, y0, t, n)
    implicit none
    !f2py intent(callback) f2
    external f2
    real(8) :: f2
    !f2py real(8) ret, arg1, arg2
    !f2py ret = f2(arg1, arg2)

    real(8), intent(in) :: y0
    real(8), intent(out) :: y(n)
    real(8), intent(in) :: t(n)
    integer n, i
    !f2py depend(t) n

    real(8) :: h, k1, k2, k3, k4

    y(1) = y0
    h = ((t(n) - t(1)) / (n-1))

    write(ou,'(//A)') "... RK4 ..."
    write(ou,'("# h = ", F12.5)') h
    write(ou,'("t       ",4x,"f(t)     ")')
    write(ou,'("--------",4x,"----------------")')

    do i = 1, n-1
      write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)
      k1 = h * f2(y(i), t(i))
      k2 = h * f2(y(i) + 0.5 * k1, t(i) + 0.5 * h)
      k3 = h * f2(y(i) + 0.5 * k2, t(i) + 0.5 * h)
      k4 = h * f2(y(i) + k3, t(i) + h)
      y(i+1) = y(i) + (k1 + 2.0 * (k2 + k3) + k4) / 6.0
    end do
    write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine rk4

!! 4th order Runge-Kutta method
subroutine rk45(y, e, y0, t, n)
    implicit none
    !f2py intent(callback) f3
    external f3
    real(8) :: f3
    !f2py real(8) ret, arg1, arg2
    !f2py ret = f3(arg1, arg2)

    real(8), intent(in) :: y0
    real(8), intent(out) :: y(n)
    real(8), intent(in) :: t(n)
	real(8), intent(out) :: e(n)
    integer n, i
    !f2py depend(t) n
    real(8) h

  real(8) :: k1, k2, k3, k4, k5, k6
  real(8) :: y5

  !! 4th order Runge-Kutta method constants
  real(8) :: c20 = 2.500000000000000e-01, &
    & c30 = 3.750000000000000e-01, &
    & c40 = 9.230769230769231e-01, &
    & c50 = 1.000000000000000e+00, &
    & c60 = 5.000000000000000e-01, &
    & c21 = 2.500000000000000e-01, &
    & c31 = 9.375000000000000e-02, &
    & c32 = 2.812500000000000e-01, &
    & c41 = 8.793809740555303e-01, &
    & c42 = -3.277196176604461e+00, &
    & c43 = 3.320892125625853e+00, &
    & c51 = 2.032407407407407e+00, &
    & c52 = -8.000000000000000e+00, &
    & c53 = 7.173489278752436e+00, &
    & c54 = -2.058966861598441e-01, &
    & c61 = -2.962962962962963e-01, &
    & c62 = 2.000000000000000e+00, &
    & c63 = -1.381676413255361e+00, &
    & c64 = 4.529727095516569e-01, &
    & c65 = -2.750000000000000e-01, &
    & a1  = 1.157407407407407e-01, &
    & a2  = 0.000000000000000e-00, &
    & a3  = 5.489278752436647e-01, &
    & a4  = 5.353313840155945e-01, &
    & a5  = -2.000000000000000e-01, &
    & b1  = 1.185185185185185e-01, &
    & b2  = 0.000000000000000e-00, &
    & b3  = 5.189863547758284e-01, &
    & b4  = 5.061314903420167e-01, &
    & b5  = -1.800000000000000e-01, &
    & b6  = 3.636363636363636e-02

    y(1) = y0
    h = ((t(n) - t(1)) / (n-1))

  write(ou, '(//A)') "... RK45 ..."
  write(ou,'("# h = ", F12.5)') h
  write(ou,'("t       ",4x,"f(t)     ")')
  write(ou,'("--------",4x,"----------------")')

  do i = 1, n-1
    write(ou,'(F12.5,4x,ES23.16)') t(i), y(i)

    k1 = h * f3(y(i), t(i))
    k2 = h * f3(y(i) + c21 * k1, t(i) + c20 * h )
    k3 = h * f3(y(i) + c31 * k1 + c32 * k2, t(i) + c30 * h )
    k4 = h * f3(y(i) + c41 * k1 + c42 * k2 + c43 * k3, t(i) + c40 * h )
    k5 = h * f3(y(i) + c51 * k1 + c52 * k2 + c53 * k3 + c54 * k4, t(i) + h )
    k6 = h * f3(y(i) + c61 * k1 + c62 * k2 + c63 * k3 + c64 * k4 + c65 * k5, &
    & t(i) + c60 * h )
    y(i+1) = y(i) + a1 * k1 + a3 * k3 + a4 * k4 + a5 * k5
    y5 = y(i) + b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6

    e(i+1) = abs( y5 - y(i+1) )
  end do
  write(ou,'(F12.5,4x,ES23.16)') t(n), y(n)

end subroutine rk45

end module ode