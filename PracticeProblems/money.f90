program moneyMaker
implicit none
!!Formula for interest rate is final = start(1+(rate/100))^N years
integer::final,start,rate,years
real::total
start=1000
years=5
open(15,file='money.txt') !Supports file creation on the fly and overwriting should the opened file already exist.  No specification for r/w
!!Now for the calculations
do rate = 2,8,2 !!do loop.  Args: i = starting, end condition, increment)
    total= start*(1+(real(rate)/real(100)))**years
    write(15,*)rate,total
end do
close(15)
stop
end program moneyMaker
