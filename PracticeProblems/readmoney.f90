program readMoney
implicit none
integer::i
real::rate,total
open(10,file='money.txt')
do i=1,4,1
read(10,*)rate,total
write(*,*)rate,total
end do
stop
end program readMoney
