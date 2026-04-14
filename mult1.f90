program mult1
implicit none
integer:: i,j,k
!
! This simple FORTRAN program multiplies two integers
! It then displays the two integers and their product
!
i=5
j=8
k=i*j
write(*,*)i,'*',j,'=',k
stop 
end program mult1