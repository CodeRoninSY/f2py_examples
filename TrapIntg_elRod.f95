subroutine TrapezoidalIntegrator(phi_x, phi_res, int_res, L)

implicit none

! IN-/OUTPUT
integer, intent(in)  :: phi_res, int_res
real, intent(in) :: L
real, intent(out) :: phi_x(phi_res)

! OTHER VARs
integer  :: i, j
real :: x, z, dx, dz, intg, chDens

dz = (L - 0.0)/real(int_res)
dx = (L - 0.0)/real(phi_res)
z = 0.0
x = 0.0

do i = 1, phi_res       ! for every x ..
    do j = 1, int_res   ! .. integrate
        chDens   = z
        intg     = chDens / ( sqrt( (x/L - z)**2.0 + 1.0) )
        phi_x(i) = phi_x(i) + intg


        z = z+dz
        chDens   = z
        intg     = chDens / (sqrt(  (x/L - z)**2.0 + 1.0 ) )
        phi_x(i) = phi_x(i) + intg

     end do
     phi_x(i) = phi_x(i) * (0.5*dz)
     x = x+dx
     z = 0.0
end do

return
end subroutine