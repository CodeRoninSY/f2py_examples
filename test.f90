subroutine test()

!f2py intent(callback) py_func
external py_func

real*8 py_func
real*8 y,x
!f2py y = py_func(x)

real*8 :: a
real*8 :: b, r

a = 12
write(6,*) a
r = py_func(a)
write(6,*) r

return
end subroutine