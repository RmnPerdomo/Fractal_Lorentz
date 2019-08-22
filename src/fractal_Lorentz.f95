program wca_fractal_bd

use random_generator
use bd_nvt_wca_clst_lib
use structure_lib
use dynamics_lib
use cpu_time_lib

! *****************************************************************************

implicit none
integer :: i, j, k, l, m, it, teq, tn, tm, block, n_dec
integer :: n, nf, frame, seed
real*8  :: box, sigma, dt, df
real*8  :: rho, phi, rho_f, phi_f, pi
real*8  :: energy

real*8,    allocatable, dimension(:)       :: x, y, z, xf, yf, zf
real*8,    allocatable, dimension(:, :)    :: t, wt
real*8,    allocatable, dimension(:, :, :) :: cfx, cfy, cfz
integer*8, allocatable, dimension(:, :)    :: nn_list, nf_list
integer*8, allocatable, dimension(:)       :: nblock, mblock

real*8 :: dx, dy, dz

character(len=2) :: id

integer :: ndiv
real*8  :: r, dr
real*8, allocatable, dimension(:) :: gr

real :: time, time0

! *****************************************************************************
! Read files

open(100, file = 'conf.xyz', status = 'unknown')
read(100, *) n, box
read(100, *) frame
allocate(x(n), y(n), z(n))
do i = 1, n
  read(100, *) id, x(i), y(i), z(i)
end do
close(100)

open(101, file = 'cluster.xyz', status = 'unknown')
read(101, *) nf, box
read(101, *) frame
allocate(xf(nf), yf(nf), zf(nf))
do i = 1, nf
  read(101, *) id, xf(i), yf(i), zf(i)
end do
close(101)

open(102, file = 'input_bd', status = 'unknown')
read(102, *) dt
read(102, *) teq
read(102, *) tm
read(102, *) tn
read(102, *) dr
read(102, *) seed
close(102)

! -----------------------------------------------------------------------------
! Parameters

sigma = 1.d0
pi    = 4.d0*datan(1.d0)
rho   = n/box**3.d0
rho_f = nf/box**3.d0
phi   = pi*rho/6.d0
phi_f = pi*rho_f/6.d0
ndiv  = int(0.5d0*box/dr)

call sgrnd(seed)

! -----------------------------------------------------------------------------
! Print information

write(*, '(/, a)')     'Mobile  particles: '
write(*, '(a, i10)')   'Particles     = ', n
write(*, '(a, f14.6)') 'Density       = ', rho
write(*, '(a, f14.6)') 'Vol. fraction = ', phi
write(*, '(a, f14.6)') 'Diameter      = ', sigma

write(*, '(/, a)')     'Fractal  particles: '
write(*, '(a, i10)')   'Particles         = ', nf
write(*, '(a, f14.6)') 'Density           = ', rho_f
write(*, '(a, f14.6)') 'Vol. fraction     = ', phi_f
write(*, '(a, f14.6)') 'Diameter          = ', sigma
write(*, '(a, f14.6)') 'Fractal dimension = ', df

write(*, '(/, a)')     'Box: '
write(*, '(a, f14.6)') 'Box lenght        = ', box

! -----------------------------------------------------------------------------

! Initial energy and forces

allocate(nn_list(n, 0:50))
allocate(nf_list(n, 0:50))
call update_verlet(n, x, y, z, box, sigma, 50, 49, nn_list, &
                   nf, xf, yf, zf, nf_list)

write(*, '(/, a)') 'Thermalization'

call set_zerotime(time)

do it = 1, teq

  if (mod(it, 20) .eq. 0) then
    call update_verlet(n, x, y, z, box, sigma, 50, 49, nn_list, &
                   nf, xf, yf, zf, nf_list)
  end if

  call bd_WCA_clst(tm, n, x, y, z, nf, xf, yf, zf, box, sigma, &
                    1.4737d0, 50, 49, dt, energy, nn_list, nf_list)

  if (mod(it, 200) .eq. 0) then
    call progress_bar_time(time, it, teq)
  end if

end do

! -----------------------------------------------------------------------------

n_dec = int(log10(tn*1.d0))

write(*, *) n_dec

allocate(gr(0: ndiv))
allocate(cfx(0:n_dec, 10, n), cfy(0:n_dec, 10, n), cfz(0:n_dec, 10, n))
allocate(t(0:n_dec, 10), wt(0:n_dec, 10))
allocate(nblock(0:n_dec))
allocate(mblock(0:n_dec))

gr  = 0.d0
cfx = 0.d0
cfy = 0.d0
cfz = 0.d0
t   = 0.d0
wt  = 0.d0

nblock = 0
mblock = 0

call set_zerotime(time)

do it = 1, tn

  if (mod(it, 20) .eq. 0) then
    call update_verlet(n, x, y, z, box, sigma, 50, 49, nn_list, &
                   nf, xf, yf, zf, nf_list)
  end if

  call bd_WCA_clst(tm, n, x, y, z, nf, xf, yf, zf, box, sigma, &
                    1.4737d0, 50, 49, dt, energy, nn_list, nf_list)

!  call gr_3d(n, x, y, z, box, dr, ndiv, gr)

  do block = 0, n_dec
    if (mod(it, 10**block) .eq. 0) then
      nblock(block) = nblock(block) + 1
      k = nblock(block)
      cfx(block, k, :) = x
      cfy(block, k, :) = y
      cfz(block, k, :) = z
      if (mod(k, 10) .eq. 0) then
        mblock(block) = mblock(block) + 1
        do l = 1, 9
          do m = 1, 10 - l
            do i = 1, n
              dx = cfx(block, l + m, i) - cfx(block, l, i)
              dy = cfy(block, l + m, i) - cfy(block, l, i)
              dz = cfz(block, l + m, i) - cfz(block, l, i)
              wt(block, m) = wt(block, m) + dx*dx + dy*dy + dz*dz
            end do
          end do
        end do
        nblock(block) = 0
      end if
    end if
  end do

  if (mod(it, 200) .eq. 0) then
    call progress_bar_time(time, it, tn)
  end if

end do

! -----------------------------------------------------------------------------
! Saving data

open(300, file = 'gr_bd.dat', status = 'unknown')
gr = gr/(4.d0*pi*rho*tn*n)
do i = 1, ndiv
  r  = i*dr
  write(300, '(2f14.6)') r, gr(i)/(r**2*dr)
end do
close(300)

open(302, file = 'msd.dat', status = 'unknown')
do block = 0, n_dec
  do i = 1, 9
    write(302, *) i*dt*10**block, wt(block, i)/(mblock(block)*n*(10 - i))
  end do
end do
close(302)

! -----------------------------------------------------------------------------
! Movie Marker

!open(400, file = 'coordinates.xyz', status = 'unknown')
!do it = 1, 200000
!
!  if (mod(it, 20) .eq. 0) then
!    call update_verlet(n, x, y, z, box, sigma, 50, 49, nn_list, &
!                   nf, xf, yf, zf, nf_list)
!  end if
!
!  call bd_WCA_clst(tm, n, x, y, z, nf, xf, yf, zf, box, sigma, &
!                    1.4737d0, 50, 49, dt, energy, nn_list, nf_list)
!
!  if (mod(it, 1000) .eq. 0) then
!    write(400, '(i10, f14.6, /)') n, box
!    x(1) = x(1) - dnint(x(1)/box)*box
!    y(1) = y(1) - dnint(y(1)/box)*box
!    z(1) = z(1) - dnint(z(1)/box)*box
!    do i = 1, n
!      x(i) = x(i) - dnint(x(i)/box)*box
!      y(i) = y(i) - dnint(y(i)/box)*box
!      z(i) = z(i) - dnint(z(i)/box)*box
!      write(400, *) 'H', x(i), y(i), z(i)
!    end do
!  end if
!end do
!close(400)

! *****************************************************************************

end program
