program matrixOperations
implicit none
!! Takes file as input for matrix values, has output file as well. Creates an allocatable matrix determined by user input.  
!! Should determine Reduced Row Echelon Form using a radix sort of values
!! NOTE: MAKE SURE THAT THE FILE IS WRITTEN IN ORDER: FIRST ROW->NthRow n=noRows, m = noCols

integer::n,m,i,j !Assigns values of row length and column height
real,dimension(:,:),allocatable::mymatrix,retMatrix
real::determinant

write(*,*)'Please input the n and m of your matrix.  FORMAT: Num Rows [ENTER] Num Cols [ENTER]' !Takes input and allocates values of the matrix to them.  
read(*,*)n,m
allocate(mymatrix(n,m)) !!Should allocate the necessary memory for an array of a given size.  Need to review booleans.  
allocate(retMatrix(n,m))

!remember that n = nRows, m = nCols
!! Insert file values into the matrix from input 10.
open(20,file='matrix2.txt', status='old')
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
!call GaussJordan(mymatrix,retMatrix,n,m) !Remember, subroutines modify the values in-place.  So RREF can be called using the result.  
!call RREF(mymatrix,retMatrix,n,m)
determinant = det(mymatrix,n,m)
write(*,*)'Determinant:', determinant !!Unpause this
!mymatrix = retMatrix

!call printMatrix(mymatrix,n,m)

!! Conducts RREF through loops.  1 in position is the last one.  Figure out how to tie it to another one or divide it such that it is 1.
!! Write them to the file as each row is finished.  
close(20)
stop
!primary program is done!

contains 
!Define print matrix subroutine; inputs are the matrix and number of rows.  Rets matrix as-is.  
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
integer::n,m,i,k
real::currVal,nextVal
real,dimension(n,m)::matrix, retMatrix
real,dimension(m)::selectedRow
retMatrix = matrix !swaps returnMatrix for the original input
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
!!Gauss Jordan portion is done.  Now bascially loop backwards in an upper triangular manner and when doing so, grab the values appropriately
!!of the scalars necessary to zero them out and then 
do k=m-1, 1, -1
    do i=k-1, 1, -1
        !write(*,*)retMatrix(i,k)
        selectedRow=retMatrix(k,:)
        nextVal = retMatrix(i,k)
        retMatrix(i,:)=retMatrix(i,:)-(selectedRow*nextVal)
    end do
end do
!Clean up negative zero values; could be a problem later on.  s
do i=1, n, 1
    do k=1, m ,1
        if (abs(retMatrix(i,k))==0) then
            retMatrix(i,k)=0
        end if
    end do
end do
end subroutine RREF

!!!Function for determinant
real function det(matrix,n,m) RESULT(r)
implicit none
integer, intent(in)::n,m
real,dimension(n,m)::matrix,workMat
real,dimension(:,:),allocatable::LMatrix
real,dimension(m)::selectedRow
real::a,b,c,d,detL,detU, currVal, nextVal, inputVal  !return for determinant
integer::i,j

if(n/=m) then
    write(*,*)'Cannot take the determinant of a nonsquare matrix.'
    return
end if
workMat = matrix
!!!Write the main gaussian upper trianguarization.  
!check if 2x2matrix
if(n==2) then
    a=workMat(1,1)
    b=workMat(1,2)
    c=workMat(2,1)
    d=workMat(2,2)
    r = ((a*d)-(c*b)) !Calculates determinant of 2x2 matrix specifically; regular algorithm will not work.  
    return
end if

!
if(n>2) then
    !Go through and construct upper triangular matrix while keeping track of diagonals as entries.  
    !Do Gauss-Jordan Elimination Method but only to the second-to-last column for U.  While doing this, keep track of the multiplicants and 
    !input those values into the L matrix in their respective positions (which is the identity matrix to that point.)
    !Then take the multiplication of the diagonals of both (i.e. the determinant) and multiply them together to get the diagonal of the whole matrix.  
    !allocate identity matrix next.  
    allocate(LMatrix(n,m), source=0.0) !allocates the memory to it and fills all values with zeroes.    
    
    !first fill L matrix with 1's along the diagonal. 
    do i=1,n,1
        LMatrix(i,i)=1.0
    end do
    !Now fill out the U (upper triangular) matrix.  
    do i=1, m-1, 1
        selectedRow = workMat(i,:)
        currVal = workMat(i,i)
        do j=i, n-1, 1
            nextVal = workMat(j+1,i) !!Nextval will be input into its place in the identity matrix.  
            inputVal=(nextVal/currVal) !inputVal is to be entered into the precise place in the iMatrix
            LMatrix(j+1,i)=inputVal
            workMat(j+1,:)=workMat(j+1,:)-(selectedRow*inputVal)
        end do
    end do
!call printMatrix(LMatrix,n,m) !L Matrix
!call printMatrix(workMat,n,m) !U Matrix
!Get the determinant of the U matrix
do i=1,n,1
    currVal=workMat(i,i)
    if(i==1) then
        detU=currVal
    else
        detU=detU*currVal
    end if
end do
!Get the determinant of the L matrix
do i=1,n,1
    currVal=LMatrix(i,i)
    if(i==1) then
        detL=currVal
    else
        detL=detL*currVal
    end if
end do
r = detU*detL

end if 

 
end function det


end program matrixOperations

!!Create a function for the determinant.  
