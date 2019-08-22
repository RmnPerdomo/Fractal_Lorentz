module   structure_lib

contains

! ------------------------------------------------------------------------------

subroutine gr_3d(n, x, y, z, box, dr, ndiv, gr)

implicit none
integer :: i, j, k, n
integer :: ndiv
real*8  :: xij, yij, zij, r, dr
real*8  :: box

real*8, dimension(0:ndiv)        :: gr
real*8, intent(in), dimension(:) :: x, y, z

do i = 1, n - 1
  do j = i + 1, n

    xij = x(i) - x(j)
    yij = y(i) - y(j)
    zij = z(i) - z(j)

    xij = xij - dnint(xij/box)*box
    yij = yij - dnint(yij/box)*box
    zij = zij - dnint(zij/box)*box

    r = dsqrt(xij*xij + yij*yij + zij*zij)
    k = dnint(r/dr)

    if (k .le. ndiv) then
      gr(k) = gr(k) + 2.d0
    endif

  enddo
enddo

return

end subroutine

! ------------------------------------------------------------------------------

end module
