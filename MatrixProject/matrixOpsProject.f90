program matrixOperations
implicit none
!! Takes file as input for matrix values, has output file as well. Creates an allocatable matrix determined by user input.  
!! Should determine Reduced Row Echelon Form using a radix sort of values
!! NOTE: MAKE SURE THAT THE FILE IS WRITTEN IN ORDER: FIRST ROW->NthRow n=noRows, m = noCols
integer::n,m,i,j !Assigns values of row length and column height
real,dimension(:,:),allocatable::mymatrix,retMatrix,inverse
real::determinant

write(*,*)'Please input the n and m of your matrix.  FORMAT: Num Rows [ENTER] Num Cols [ENTER]' !Takes input and allocates values of the matrix to them.  
read(*,*)n,m
allocate(mymatrix(n,m)) !!Should allocate the necessary memory for an array of a given size.  Need to review booleans.  
allocate(retMatrix(n,m))
allocate(inverse(n,m))

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

!! Step 1: Get Gauss-Jordan Elim, RREF, Determinant, Identity - Write subroutines or functions for these.  
!! Step 2: Span, Rank, Inverse.   
!! Step 3: Eigenvalues/Eigenmatrices

!prints matrix first by rows, then by columns
!FOR ENTRY OF MATRIX VIA TXT: Format is: Each subsequent value is a column.  It goes first line, second line third line; those are
!three columns of the first row.  Then third line fourth line fifth line will be columns for the second row.  
!It must be entered as 2 (num rows) and then 3 (num columns).  

!mymatrix = transpose(mymatrix)
write(*,*)'Input Matrix:'
call printMatrix(mymatrix,n,m)
!call GaussJordan(mymatrix,retMatrix,n,m) !Remember, subroutines modify the values in-place.  So RREF can be called using the result.  
!call RREF(mymatrix,retMatrix,n,m)

determinant = det(mymatrix,n,m)
write(*,*)'The Determinant of the Matrix is:', determinant
write(*,*)
inverse = inv(mymatrix,n,m)
write(*,*)'The Inverse of the Matrix is:'
call printMatrix(inverse,n,m)
!mymatrix = retMatrix

!call printMatrix(mymatrix,n,m)
close(20)
deallocate(mymatrix)
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
real,dimension(n,m), intent(in)::matrix
real, dimension(n,m)::workMat
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

!!Find the inverse of a matrix using a function to return the pure inverse.  This is done using LU Decomposition
function inv(matrix,n,m) RESULT(invMat)
implicit none

integer, intent(in)::n,m
real,dimension(n,m), intent(in)::matrix
real,dimension(n,m)::UMatrix,invMat,iMatrix,iMatrix2 !invMat is return value
real,dimension(:,:),allocatable::LMatrix
real,dimension(m)::selectedRow,selectedRow2 !second specifically for U inverse.  
real::currVal,nextVal,inputVal

UMatrix = matrix !Assigns UMat
allocate(LMatrix(n,m), source=0.0) !allocates the memory to it and fills all values with zeroes.  
do i=1,n,1
        LMatrix(i,i)=1.0
end do
iMatrix = LMatrix !Assigns identity matrix for later use.  
iMatrix2 = LMatrix
!Now fill out the U (upper triangular) matrix.  
do i=1, m-1, 1
    selectedRow = UMatrix(i,:)
    currVal = UMatrix(i,i)
    do j=i, n-1, 1
        nextVal = UMatrix(j+1,i) !!Nextval will be input into its place in the identity matrix.  
        inputVal=(nextVal/currVal) !inputVal is to be entered into the precise place in the iMatrix
        LMatrix(j+1,i)=inputVal
        UMatrix(j+1,:)=UMatrix(j+1,:)-(selectedRow*inputVal)
    end do
end do
!!Once U and L Matrices have been gotten, and now they should be, take the inverse of the L Matrix and U Matrix.  
!Inverse L (L will always have 1's on the diagonal and so inverse will just be negative of lower triangular values)
!use iMatrix2
do i=1,n,1
    selectedRow=LMatrix(i,:)
    selectedRow2=iMatrix2(i,:)

    currVal = LMatrix(i,i)
    do j=i+1,n,1
        !write(*,*)LMatrix(j,i) !grabs correct values
        nextVal = LMatrix(j,i) 
        LMatrix(j,:)=LMatrix(j,:)-selectedRow*nextVal
        iMatrix2(j,:)=iMatrix2(j,:)-selectedRow2*nextVal
    end do
end do
LMatrix = iMatrix2
!Inverse U, use iMat to set up the whole thing.  Reduce one to all 1's while doing the same movements on the other.  
!start from (n,n) and then work up the nth column until the first row has been reached.  Then go in and reduce the column by 1 and reduce the row also.  
do i=n,1,-1
    currVal=UMatrix(i,i) !Assign currval to what the identity matrix has to be divided by.  
    iMatrix(i,:)=iMatrix(i,:)/currVal !divide identity matrix by value likewise
    UMatrix(i,:)=UMatrix(i,:)/currVal
    selectedRow=UMatrix(i,:)
    selectedRow2=iMatrix(i,:)
    do j=i-1,1,-1 !Moves up the rows in specific column
        nextVal=UMatrix(j,i)!nextval now grabs correct value
        UMatrix(j,:)=UMatrix(j,:)-selectedRow*nextVal
        iMatrix(j,:)=iMatrix(j,:)-selectedRow2*nextVal 
    end do
end do
UMatrix = iMatrix
!call printMatrix(LMatrix,n,m)
!call printMatrix(iMatrix,n,m)
!Now that both U and L are inversed as LMatrix and UMatrix, the actual matrix's inverse can be found by multiplying the two together.  
invMat = matmul(UMatrix,LMatrix)

end function inv
end program matrixOperations
