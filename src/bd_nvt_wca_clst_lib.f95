module bd_nvt_wca_clst_lib

use random_generator

contains

!******************************************************************

subroutine bd_WCA_clst(tn, n, x, y, z, nf, xf, yf, zf, box, sigma, &
                       temp, lr, la, dt, energy, nn_list, nf_list)
implicit none
integer :: i, j, it, tn
integer :: n, nf
integer :: lr, la
real*8  :: dx, dy, dz, dt
real*8  :: box, sigma, temp, energy
real*8  :: sig
real*8, dimension(n)  :: x, y, z, Fx, Fy, Fz
real*8, dimension(nf) :: xf, yf, zf
integer*8, dimension(n, 0:50)  :: nn_list
integer*8, dimension(n, 0:50) :: nf_list

sig = sqrt(2.d0*dt)

do it = 1, tn

  call Force_WCA_clst(n, x, y, z, nf, xf, yf, zf, box, sigma, &
                      temp, lr, la, Fx, Fy, Fz, energy, nn_list, nf_list)

  do i = 1, n
    dx = sig*gasdev()
    dy = sig*gasdev()
    dz = sig*gasdev()

    x(i) = x(i) + Fx(i)*dt + dx
    y(i) = y(i) + Fy(i)*dt + dy
    z(i) = z(i) + Fz(i)*dt + dz
  end do
end do

return

end subroutine bd_WCA_clst

!------------------------------------------------------------------

subroutine Force_WCA_clst(n, x, y, z, nf, xf, yf, zf, box, sigma, &
                          temp, lr, la, Fx, Fy, Fz, energy, nn_list, nf_list)

implicit none
integer :: i, j, k, n, nf
integer :: lr, la
real*8  :: box, sigma, energy, force, temp, A, B
real*8  :: xij, yij, zij, rij, uij, fij
real*8, dimension(n)  :: x, y, z, Fx, Fy, Fz
real*8, dimension(nf) :: xf, yf, zf

integer*8, dimension(n, 0:50) :: nn_list
integer*8, dimension(n, 0:50) :: nf_list

energy = 0.d0
Fx = 0.d0
Fy = 0.d0
Fz = 0.d0

A = 1.d0*lr*(1.d0*lr/la)**la
B = 1.d0*lr/la

do i = 1, n
  do k = 1, nf_list(i, 0)
   
    j = nf_list(i, k)

    xij = x(i) - xf(j)
    yij = y(i) - yf(j)
    zij = z(i) - zf(j)

    xij = xij - box*anint(xij/box)
    yij = yij - box*anint(yij/box)
    zij = zij - box*anint(zij/box)

    rij = sqrt(xij*xij + yij*yij + zij*zij)

    if (rij .lt. sigma*B) then
      uij = A*((sigma/rij)**lr - (sigma/rij)**la) + 1.d0
      fij = A*(lr*(sigma/rij)**(lr + 1.d0) - la*(sigma/rij)**(la + 1.d0))

      energy = energy + uij/temp
      force  = fij/temp

      Fx(i) = Fx(i) + force*(xij/rij)
      Fy(i) = Fy(i) + force*(yij/rij)
      Fz(i) = Fz(i) + force*(zij/rij)
      if (uij .gt. 100000) then
        write(*, *) 'Error in fractal'
        write(*, *) i, j, rij
      end if
    end if
 
  end do
end do

end subroutine Force_WCA_clst

!******************************************************************

subroutine update_verlet(n, x, y, z, box, sigma, lr, la, nn_list, &
                         nf, xf, yf, zf, nf_list)

implicit none
integer :: i, j, k, n, nf, lr, la
real*8  :: box, sigma
real*8  :: r, xij, yij, zij

real*8, dimension(n) :: x, y, z
real*8, dimension(nf) :: xf, yf, zf
integer*8, dimension(n, 0:50) :: nn_list
integer*8, dimension(n, 0:50) :: nf_list

nn_list = 0
nn_list(:, 0) = 0

nf_list = 0
nf_list(:, 0) = 0

do i = 1, n
  do j = 1, nf

    xij = x(i) - xf(j)
    yij = y(i) - yf(j)
    zij = z(i) - zf(j)

    xij = xij - box*anint(xij/box)
    yij = yij - box*anint(yij/box)
    zij = zij - box*anint(zij/box)

    r   = sqrt(xij*xij + yij*yij + zij*zij)

    if (r .le. 1.20*sigma) then
      nf_list(i, 0) = nf_list(i, 0) + 1
      nf_list(i, nf_list(i, 0)) = j
      if (nf_list(i, 0) .gt. 50) then
        write(*, *) 'Error in fractal neighbors'
        stop
      end if
    endif

  end do
end do

end subroutine update_verlet

!******************************************************************

end module

