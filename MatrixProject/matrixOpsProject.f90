program matrices2
implicit none
!! Takes file as input for matrix values, has output file as well. Creates an allocatable matrix determined by user input.  
!! Should determine Reduced Row Echelon Form using a radix sort of values
!! NOTE: MAKE SURE THAT THE FILE IS WRITTEN IN ORDER: FIRST ROW->NthRow n=noRows, m = noCols

integer::n,m,i,j !Assigns values of row length and column height
real,dimension(:,:),allocatable::mymatrix,retMatrix

write(*,*)'Please input the n and m of your matrix.  FORMAT: Num Rows [ENTER] Num Cols [ENTER]' !Takes input and allocates values of the matrix to them.  
read(*,*)n,m
allocate(mymatrix(n,m)) !!Should allocate the necessary memory for an array of a given size.  Need to review booleans.  
allocate(retMatrix(n,m))

!remember that n = nRows, m = nCols
!! Insert file values into the matrix from input 10.
open(20,file='sample_matrix_input.txt', status='old')
!! Create file-read method for this.  Basically just from one input to another, or allow multi-input?  File works better I think.  
do i=1,n,1
    do j=1,m,1
        !!!loop through and fill matrix values in.  
        read(20,*)mymatrix(i,j)
    end do
end do
!! From there present the options of what can be done, such as finding the span, inverse, etc.  


!! Step 1: Get RREF & Also Determinant Check - Write subroutines for these.  
!! Step 2: Span, Inverse, & Identity Matrix.  
!! Step 3: 

!prints matrix first by rows, then by columns
!FOR ENTRY OF MATRIX VIA TXT: Format is: Each subsequent value is a column.  It goes first line, second line third line; those are
!three columns of the first row.  Then third line fourth line fifth line will be columns for the second row.  
!It must be entered as 2 (num rows) and then 3 (num columns).  

!mymatrix = transpose(mymatrix)
call printMatrix(mymatrix,n,m)
call GaussJordan(mymatrix,retMatrix,n,m) !Remember, subroutines modify the values in-place.  So RREF can be called using the result.  
mymatrix = retMatrix

call printMatrix(mymatrix,n,m)
!! Conducts RREF through loops.  1 in position is the last one.  Figure out how to tie it to another one or divide it such that it is 1.
!! Write them to the file as each row is finished.  
close(20)
stop
!primary program is done!
!Define print matrix function; inputs are the matrix and number of rows.
contains 
subroutine printMatrix(matrix,n,m)
implicit none
integer::n,m,k
real,dimension(n,m)::matrix
do k = 1, n
    write(*,*)matrix(k,:)
end do
write(*,*)
return 
end subroutine printMatrix

!!!Create GJ Subroutine here:
subroutine GaussJordan(matrix,retMatrix,n,m) 
implicit none
integer::n,m,i,k
real::currVal,nextVal
real,dimension(n,m)::matrix, retMatrix
real,dimension(m)::selectedRow
retMatrix = matrix !assign matrix value to keeper, avoids conflict.  
!! Create a way to find the upper rightmost corner.  
!! Divide the first value (1,1) for upper corner by itself 
!! Then select that row and iterate through all rows below it until it reaches the last one n-1 necessary to avoid bogus values to use as scalar.  
!! Then take scalar and multiply by the selected row and subtract that combination from the rows beneath the selected row. 
do i=1, m-1, 1
    currVal=retMatrix(i,i)
    retMatrix(i,:)=retMatrix(i,:)/currVal !Reduces the row to it's leading 1 form.  
    selectedRow = retMatrix(i,:)
    !find succeeding values in the next row (needs a ceiling)
    do k=i, n-1, 1
        nextVal = retMatrix(k+1,i) 
        !varRow = retMatrix(k+1,:)
        retMatrix(k+1,:) = retMatrix(k+1,:)-(selectedRow*nextVal)
    end do
end do
return
end subroutine GaussJordan

!Remember, a subroutine returns the value in place.  Return var specified at the top. 
subroutine RREF(matrix, retMatrix, n, m)
implicit none
integer::n,m
real,dimension(n,m)::matrix, retMatrix
retMatrix = matrix !swaps returnMatrix for the original input

end subroutine RREF
end program matrices2
