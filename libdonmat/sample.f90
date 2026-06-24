program libdonmatsample
use libdonmat
implicit none

!Declare values below
integer::n,m
real,dimension(:,:),allocatable::mymatrix

!Get user input
write(*,*)'Please enter the number of rows and then the number of columns'
read(*,*)n,m

!Allocate to fit that input
allocate(mymatrix(n,m))

!Function Call
mymatrix=readmatrix('matrix2.txt',n,m)

!Subroutine Call
call printMatrix(mymatrix)
end program libdonmatsample
