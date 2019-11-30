subroutine test(r, x, t, n, m)
! test.f90
! r : real(8) (OUT) -> return value
! x, t: real(8) (IN) -> for func args
! m,n: integer -> dimensions of x & t
! py_func : external f; pass from python to fortran

!f2py intent(callback) f
external f
! The following lines define the signature of py_func for f2py
real(8) f
!f2py real(8) y,arg1,arg2
!f2py y = f(arg1, arg2)

real(8), intent(in) :: x(m)
real(8), intent(in) :: t(n)
integer :: n, m
!f2py depend(t) n
!f2py depend(x) m
real(8), INTENT(OUT) :: r(m, n)

write(6,'(A4,2x,3F16.8)') "x: ", (x(i),i=1,m)
write(6,'(A)') "i       j       t               r"

do j = 1,m
    do i = 1,n
        r(j, i) = f(x(j), t(i))
        write(6,'(I4,2x,I4,2x,F16.8,2x,E35.28)') i,j, t(i), r(j, i)
    end do
end do

end subroutine