program readFile
integer::i,j,k
open(11,file='info.txt')
read(11,*)k,l,m !takes data line by line and assigns them automatically to variables as per the standard denoted by *.
!write(*,*)k,l,m
write(11,*)'Hello, I am on line 4'
i=5
j=1000
k=52345617
11 format(i2,i4,i10)
write(11,*)i,j,k
close(11)
stop
end program readFile