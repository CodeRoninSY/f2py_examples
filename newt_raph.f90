! $UWHPSC/codes/fortran/newton/newton.f90
! https://staff.washington.edu/rjl/classes/am583s2014/notes/fortran_newton.html

module newton

    ! module parameters:
    implicit none
    integer, parameter :: maxiter = 100
    real(8), parameter :: tol = 1.0e-5

contains

subroutine solve(f, fp, x0, x, iters, debug)

    ! Estimate the zero of f(x) using Newton's method.
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!)
    !   the number of iterations iters

    implicit none
    !f2py intent(in) :: x0
    !f2py intent(in) :: debug
    !f2py intent(in,out) :: x
    !f2py intent(in,out) :: iters
    external f
    external fp
    real(8) :: f, fp

    real(8) :: x0
    logical :: debug
    real(8) :: x
    integer :: iters

    ! Declare any local variables:
    real(8) :: deltax, fx, fxprime
    integer :: k

    ! initial guess
    x = x0

    if (debug) then
        print 11, x
 11     format(//'Initial guess: x = ', es23.16)
    endif

    ! Newton iteration to find a zero of f(x)

    do k = 1,maxiter

        ! evaluate function and its derivative:
        fx = f(x)
        fxprime = fp(x)

        if (abs(fx) < tol) then
            exit  ! jump out of do loop
            endif

        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        x = x - deltax

        if (debug) then
            print 12, k,x
 12         format('After', i4, ' iterations, x = ', es23.16)
        endif

    enddo


    if (k > maxiter) then
        ! might not have converged

        fx = f(x)
        if (abs(fx) > tol) then
            print *, '*** Warning: has not yet converged'
            endif
        endif

    ! number of iterations taken:
    iters = k - 1

end subroutine solve

end module newton