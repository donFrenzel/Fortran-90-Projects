program roots1
implicit none
integer::i,n
real::rooti
!Now put the do loop.  Always
write(*,*)'Enter an integer'
read(*,*)n
do i=2*n,2,-2 !Start increment variable condition, end condition, and increment.  
    rooti=sqrt(real(i))
    write(*,*)i,rooti
end do
stop
end program roots1
