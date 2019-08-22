module dynamics_lib

contains
! *****************************************************************************

subroutine diffusion(tn, n, cfx, cfy, cfz, wt)

implicit none
integer :: i, j, k, tn, n
real*8  :: dx, dy, dz, msd
real*8, dimension(tn, n) :: cfx, cfy, cfz
real*8, dimension(tn)    :: t, wt

wt = 0.d0

do i = 1, tn - 1
  msd = 0.d0
  do j = 1, tn - i
    do k = 1, n
      dx  = cfx(j + i, k) - cfx(j, k)
      dy  = cfy(j + i, k) - cfy(j, k)
      dz  = cfz(j + i, k) - cfz(j, k)
      msd = msd + dx*dx + dy*dy + dz*dz
    end do
  end do
  wt(i) = wt(i) + msd/(n*(tn - i))
end do

return

end subroutine diffusion

! *****************************************************************************
end module
