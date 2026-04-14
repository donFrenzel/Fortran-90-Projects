program trig
implicit none
real::inputdeg,sin1,cos1,tan1
write(*,*)'Please input an angle in degrees'
read(*,*)inputdeg
sin1=sin(inputdeg)
cos1=cos(inputdeg)
tan1=tan(inputdeg)
write(*,*)sin1,cos1,tan1

!character*5,dimension(5):: students=(Fred,Susie,Tom,Anita,Peter)
!write(*,*)students

stop
end program trig
