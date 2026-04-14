program findinv
implicit none
real::j,x
write(*,*)'Enter the integer whose inverse you wish to find.'
read(*,*)j
x=-j
write(*,*)'The inverse of ',j,'is ',x
stop
end program findinv