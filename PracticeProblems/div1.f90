program div1
implicit none
integer::i,j
real::k
write(*,*)'Please enter a divisior, then a dividend:'
read(*,*)i,j
k=real(i)/real(j)
write(*,*)'The result is: ',k
stop
end program div1
