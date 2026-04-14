program arrays
!program which adds two vectors
implicit none
real,dimension(3):: vect1, vect2, vectsum
!!Read in the elements of the two vectors
write(*,*)'Enter the three components of vector 1'
read(*,*)vect1
write(*,*)'Enter the three components of vector 2'
read(*,*)vect2
!!now add the vectors together by adding their components
vectsum = vect1 + vect2
write(*,10)vectsum
10 format('The sum of the two vectors is:'3f6.2) !!Not sure what the 3 means 3 potential outputs, f means format, 6 means spaces from the str finish
!!and the .2 means Two decimal places
close(10)
write(*,11)maxval(vect1)+maxval(vect2) !!Allows maxval functions.  Very useful.  
11 format('The sum of the maximum values of the arrays are:'1f6.2)
close(11)
write(*,12)sqrt(sum(vectsum**2))
12 format('The normal (length) of vectsum is:'1f6.2)
close(12)
end program arrays