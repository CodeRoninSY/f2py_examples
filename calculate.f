      subroutine calculate(x,n)
      implicit none

cf2py intent(callback) func
      external func
      real*8 func
c     The following lines define the signature of func for F2PY:
      real*8 y, r
cf2py y = func(r)
c
cf2py intent(in,out,copy) x
      integer n,i
      real*8 x(n)
      do i=1,n
         x(i) = func(x(i))
      end do
      end