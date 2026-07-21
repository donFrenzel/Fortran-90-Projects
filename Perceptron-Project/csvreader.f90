module fortranCSVRead
implicit none
contains
!fill here for module.
subroutine csv_to_txt(filename,filenameC,rows,cols) 
implicit none 
character(len=*), intent(in)::filename
integer::rows,cols
integer::eof,suffixPos,i,commaPosition,startPos
integer::numCols,numRows
character(len=4)::oldSuf,newSuf
character(len=255)::filenameC,field
character(len=2048)::firstLine,curLine !for first line of the file and current file line (for dividing up and entering as data into txt.  )
!define old and new suffixes for transfusion
suffixPos=0
oldSuf=".csv"
newSuf=".txt"
filenameC=filename
suffixPos=index(filenameC,oldSuf)!find position of the old .csv suffix AND replace with new suffix, new filename to open up. 
if (suffixPos.gt.0) then 
    filenameC=filenameC(:suffixPos-1)//newSuf//filenameC(suffixPos+len(oldSuf):)
end if
eof=0
numCols=1 !will always have at least 1 data column for the csv not to be empty.  If you feed it an empty CSV its your fault.  
numRows=0
open(15,file=filename,status='old')
open(20,file=filenameC,status='unknown',action="write")
read(15, '(A)', IOSTAT=eof) firstLine !read to the first line.  
!NOTE: len_trim() intrinsic is the trimmed (no nullspace) length of the character array being entered.  Loop thru and get comma delimit count
do i=1,len_trim(firstLine)
    if(firstLine(i:i)==',') then 
        numCols=numCols+1
    end if
end do
do 
    read(15,'(A)',iostat=eof)curLine !now we need to take the number of columns and separate the data into those x columns PER UNIT
    if(eof<0) then 
        exit
    end if
    startPos=1
    commaPosition=1
    !now that it is curLine holding the line data, we take the numCols and divide the data up into them by finding the index
    do
        !get the index of the comma position and check if there is no next comma, in which case just grab the data until the end. 
        commaPosition=index(curLine(startPos:),',')
        if(commaPosition==0) then 
            field=curLine(startPos:)
            write(20,*)trim(adjustl(field)) !write to the txt the last bit (since that's what we're at in this if statement)
            exit
        else
            !now basically we're gonna need to get the comma position and read everything from that until the next comma.  
            commaPosition= startPos+commaPosition-1 !basically here because you're slicing a small portion with INDEX() above, must be added to 
            !total position.  Caused many problems.  Not any more.  
            field=curLine(startPos:commaPosition-1) !grab between the position of the start and the comma (-1 necessary because it does index that value.  )
            write(20,*)trim(adjustl(field)) 
            startPos=commaPosition+1
        end if
    end do
    numRows=numRows+1
end do
rows=numRows
cols=numCols
!close out the connections.
close(15)
close(20)
end subroutine csv_to_txt

end module fortranCSVRead
!declare program header here.  
program csvreader
use fortranCSVRead
use libdonmat
implicit none
integer::numRows,numCols
character(len=255)::txtname
real,dimension(:,:),allocatable::mymatrix
call csv_to_txt('sampleData.csv',txtname,numRows,numCols)
write(*,*)'Num Rows: ',numRows
write(*,*)'Num Cols: ',numCols
write(*,*)'Txt file name: ',txtname
!now that we've made a txt file we can use our libdonmat functions to send it!!
!TO COMPILE USE: gfortran -I. -o csvreader.exe csvreader.f90 .\libdonmat.a
!TO RUN USE: ./csvreader.exe
allocate(mymatrix(numRows,numCols))
mymatrix=readmatrix(txtname,numRows,numCols)
call printMatrix(mymatrix)
deallocate(mymatrix)
end program csvreader
