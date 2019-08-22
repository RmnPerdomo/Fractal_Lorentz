program dynamical_exponent

implicit none
integer :: i, j
real*8, allocatable, dimension(:) :: t, msd

allocate(t(63), msd(63))

open(100, file = 'msd.dat', status = 'old')
do i = 1, 63
  read(100, *) t(i), msd(i)
end do
close(100)

open(200, file = 'log10msd.dat', status = 'unknown')
do i = 1, 63
  write(200, *) log10(t(i)), log10(msd(i))
end do
close(200)

t = log10(t)
msd = log10(msd)

open(300, file = 'exponent.dat', status = 'unknown')
do i = 2, 62
  write(300, *) t(i), (msd(i + 1) - msd(i - 1))/(t(i + 1) - t(i - 1))
end do
close(300)


! Agregado a GitHub el 22 de agosto del 2019

end program
