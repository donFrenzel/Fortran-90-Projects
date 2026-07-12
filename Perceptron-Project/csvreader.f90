module fortranCSVRead
implicit none
contains
!fill here for module.

subroutine readCSVPerceptron(filename) 
implicit none 
character(len=*), intent(in)::filename
integer::label,x1,x2,eof,numIters
eof=0
!open the file
open(15,file=filename,status='old')
numIters=0
do 
    !handle case for the labels/headers since the first iteration will always have them
    if(numIters==0) then 
        read(15,*)
        numIters=numIters+1
        cycle
    end if
    read(15,*,iostat=eof)label,x1,x2 !just use the standard form bc its already comma delimited, just as the needs of the csv are. 
    !handle case where end of file is met
    if(eof<0) then 
        exit
    end if
    print *,'label',label,'X1: ',x1,'X2: ',x2
    numIters=numIters+1
end do
close(15)
end subroutine readCSVPerceptron

end module fortranCSVRead
!declare program header here.  
program csvreader
use fortranCSVRead
implicit none
call readCSVPerceptron('sampleData.csv')
end program csvreader
