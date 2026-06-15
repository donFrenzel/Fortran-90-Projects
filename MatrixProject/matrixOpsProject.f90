program matrixOperations
implicit none
!! Takes file as input for matrix values, has output file as well. Creates an allocatable matrix determined by user input.  
!! Should determine Reduced Row Echelon Form using a radix sort of values
!! NOTE: MAKE SURE THAT THE FILE IS WRITTEN IN ORDER: FIRST ROW->NthRow n=noRows, m = noCols
integer::n,m,i,j,retty !Assigns values of row length and column height
real,dimension(:,:),allocatable::mymatrix,returnMatrix,inv,Q,R,reye,outerprod,retMat2, U, VT, S
real,dimension(:),allocatable::vec,eigs,vec2,vec3,S1
real::determinant,normal

write(*,*)'Please input the n and m of your matrix.  FORMAT: Num Rows [ENTER] Num Cols [ENTER]' !Takes input and allocates values of the matrix to them.  
read(*,*)n,m
allocate(mymatrix(n,m)) !!Should allocate the necessary memory for an array of a given size.  Need to review booleans.  
allocate(returnMatrix(n,m))
allocate(inv(n,m))
allocate(vec(m))
allocate(Q(m,m))
allocate(R(n,m))
allocate(reye(n,m))
allocate(eigs(m))
allocate(vec2(n))
allocate(vec3(n))
allocate(outerProd(m,n))
allocate(retMat2(n,m))

allocate(U(m,m))
allocate(S1(m))
allocate(S(m,n))
allocate(VT(n,n))

!remember that n = nRows, m = nCols
!! Insert file values into the matrix from input 10.
open(20,file='matrix4.txt', status='old')
!! Create file-read method for this.  Basically just from one input to another, or allow multi-input?  File works better I think.  
do i=1,n,1
    do j=1,m,1
        !!!loop through and fill matrix values in.  
        read(20,*)mymatrix(i,j)
    end do
end do
!! From there present the options of what can be done, such as finding the span, inverse, etc.  

!! Step 1: Get Gauss-Jordan Elim, RREF - Write subroutines or functions for these.  
!! Step 2: Span, Rank, Inverse[X], Determinant[X], Identity[X].   
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

!returnMatrix = GaussJordan(mymatrix,n,m)
!write(*,*)'Gaussian Elim Matrix:'
!call printMatrix(returnMatrix,n,m)
!write(*,*)

returnMatrix = RREF(mymatrix) !switch this back
write(*,*)'RREF of Matrix:'
call printMatrix(returnMatrix,n,m)
write(*,*)

if (n==m) then
determinant = det(mymatrix,n,m)
write(*,*)'The Determinant of the Matrix is:', determinant
write(*,*)

inv = inverse(mymatrix,n,m)
write(*,*)'The Inverse of the Matrix is:'
call printMatrix(inv,n,m)
!mymatrix = retMatrix
end if

retty = linearIndependence(mymatrix,n,m)
!call QR(mymatrix,n,m,Q,R)
!call printMatrix(Q,m,m)
!call printMatrix(R,n,m)


!!First vector of the matrix
vec = mymatrix(1,:)
normal = Norm(vec)
write(*,*)'Normal of the first vector:',normal

vec2 = mymatrix(:,1)
vec3 = myMatrix(:,2)
outerprod = outerProduct(vec2,vec3)

!!Test Hessenberg Reduction Function here: 
!eigs = eigenvals(mymatrix)
!write(*,*)'The Eigenvalues of the Matrix are:'
!write(*,*)eigs

call SVD2(myMatrix,U,S,VT)

!call printMatrix(mymatrix,n,m)
close(20)
deallocate(mymatrix)
stop
!primary program is done!
contains 

!!!BASIC FUNCTIONS INCLUDE printMatrix, norm, and custom sign function.  
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

!!Function Takes the Norm
real function Norm(vector) RESULT(vectNormal)
integer::n,i
real, intent(in)::vector(:)
real::sum
n = size(vector)
sum=0.0
!Take square root of the squared sum.  
do i=1,n,1
    sum=sum+(vector(i)**2)
end do 
vectNormal = sqrt(sum)
end function Norm

!!!sign function: Return 1 if positive or 0 val, return -1 if negative.
function sgn(val) RESULT(retVal)
implicit none
real::val,retVal
if (val.ge.0.0) then 
    retVal=1.0
else 
    retVal=-1.0
end if
end function sgn

!Create a givens function here, 'returns' modified values of c and s, they are necessary as inputs.  
subroutine givens(a,b,c,s)
implicit none
real, intent(in)::a,b
real, intent(out)::c,s 
real::t
if (b==0.0) then 
    c=1.0
    s=0.0
else
    if(abs(b)>abs(a)) then 
        t=-a/b
        s=1/sqrt(1+t**2)
        c=s*t
    else
        t=-b/a
        c=1/sqrt(1+t**2)
        s=c*t
    end if
end if
end subroutine givens

!Function for rounding values.  
function roundSmalls(inMatrix,precision) RESULT(outMatrix)
implicit none
real, intent(in),optional::precision
real, intent(in)::inMatrix(:,:)
integer::n,m,i,j
real::actPrecision
real, dimension(:,:),allocatable::outMatrix
!Check first if precison before allocating.  Autoset value if not.  
if(present(precision)) then
    actPrecision=precision
else
    actPrecision = 0.00001 !aka 1e-5
end if
!!Get size and allocate output as that value.  
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
allocate(outMatrix(n,m))
outMatrix = inMatrix !Make sure the two are equal to start.  
!Check firs
!!Now comb through and find the small values.  Take abs(val @ i,j) and then see if its greater than or less than the current one.  If less, round to zero.
do i=1,n,1
    do j=1,m,1
        !Check if absolute value (positive or negative is less than precision.  if so, round it)
        if (abs(outMatrix(i,j))<actPrecision) then 
            outMatrix(i,j) = 0.0
        else 
            continue
        end if
    end do 
end do
end function roundSmalls

!Eye function analog to numpy.eye
function eye(n,m,rightShift,downShift) RESULT(eyeMatrix)
implicit none
integer, intent(in)::n,m
integer, intent(in), optional::rightShift,downShift
integer::actRightShift,actDownShift,j
real,dimension(n,m)::eyeMatrix
eyeMatrix=0
if (present(rightShift)) then
    actRightShift=rightShift
else
    actRightShift=0
end if
if (present(downShift)) then 
    actDownShift=downShift
else
    actDownShift=0
end if
i=1+actDownShift
do j=1+actRightShift,m,1
    eyeMatrix(i,j)=1.0
    i=i+1
end do
end function

!Simple sorting algorithm (roughly O(n))
function simpleSort(array) RESULT(outArray)
implicit none
real, intent(in)::array(:)
real, dimension(:), allocatable::curArray,outArray, temp
integer::i,n,l
integer,dimension(1)::mloc
real::max
n = size(array)
allocate(outArray(n))
allocate(curArray(n))
curArray = array
do i=1,n,1
    mloc = maxloc(curArray) !max location
    max = curArray(mloc(1)) !max val
    l = size(curArray)
    !allocate temp if unallocated, else deallocate. (Will always be allocated.  )
    if (allocated(temp)) then 
        deallocate(temp)
    end if
    allocate(temp(l-1))
    !Use slicing here.  
    temp = (/curArray(1:mloc(1)-1), curArray(mloc(1)+1:l)/)
    !Reallocate current, modifiable array to be 1 smaller so that the size fits.  
    if(allocated(curArray)) then 
        deallocate(curArray)
    end if
    allocate(curArray(l-1))
    !Assign newly allocated array the values of temp (which has now excluded the old max)
    curArray = temp
    outArray(i)=max
end do
end function simpleSort

function rowswap(inMatrix,rowInd1,rowInd2) RESULT(outMatrix)
integer, intent(in)::rowInd1,rowInd2
real, intent(in)::inMatrix(:,:)
integer::n,m
real, dimension(:,:),allocatable::outMatrix
real, dimension(:),allocatable::row1, row2, temp
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2)
allocate(outMatrix(n,m))
allocate(row1(m))
allocate(row2(m))
allocate(temp(m))
outMatrix = inMatrix
row1 = outMatrix(rowInd1,:)
row2 = outMatrix(rowInd2,:)
outMatrix(rowInd2,:) = row1
outMatrix(rowInd1,:) = row2
end function rowswap

!END Simple operations.  
!!!MATRIX OPERATIONS BEGIN HERE:  

!!!Create GJ Subroutine here:
function GaussJordan(matrix,n,m) RESULT(retMatrix)
implicit none
integer, intent(in)::n,m
integer::i,k
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
end function GaussJordan
!Remember, a subroutine returns the value in place.  Return var specified at the top. 
function RREF(inMatrix) RESULT(outMatrix)
implicit none
real, intent(in)::inMatrix(:,:)
integer::n,m,i,j,pivotR,nextPivotR
real::pivotVal,curVal
real, dimension(:,:), allocatable::outMatrix
real, dimension(:), allocatable::curRow,nextRow
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
allocate(outMatrix(n,m))
allocate(curRow(m))
allocate(nextRow(m))
pivotR=1
nextPivotR=1
outMatrix = inMatrix
!write(*,*)'This is now running!'
!Step one, keep index of the pivot row.  Currently it is 1.
do j=1,m,1
    !sends all zero rows to the bottom for every step. 
    !first check if the matrix is empty, in which case it should return itself
    if (sum(outMatrix)==0.0) then
        exit
    end if
    !move all zero rows to the bottom per column iteration.  
    do i=1,n-1,1
        curRow = outMatrix(i,:)
        nextRow = outMatrix(i+1,:)
        if(sum(curRow)==0.0) then 
            !send this to the bottom
            if (sum(nextRow)/=0.0) then 
                outMatrix(i,:) = outMatrix(i+1,:)
                outMatrix(i+1,:) = curRow
            end if
        end if
    end do
    !now that all zero rows are automatically sent to the bottom, find the first column which is a pivot.  
    !current pivot position is 1, aka first position.  after finding the next pivot the two should swap. 
    nextPivotR = 0
    do i=pivotR,n,1
        if(outMatrix(i,j)/=0.0) then 
            !we know in this case that our current pivot is GOOD
            nextPivotR=i
            exit !exit the loop.  
        end if 
    end do
    !cycle to next value in the loop if there are no valid pivots found.  
    if (nextPivotR==0.0) then 
        cycle
    end if
    !so now we have two pivot row indices, one for the current pivot (in this case 1) and another for the next one, currently next.  
    !now, if nextPivotR.lt.pivotR then swap them.  
    if (pivotR.lt.nextPivotR) then
        curRow = outMatrix(nextPivotR,:)
        outMatrix(nextPivotR,:) = outMatrix(pivotR,:)
        outMatrix(pivotR,:) = curRow
    end if
    !now that we have the pivot row
    pivotVal = outMatrix(pivotR,j)
    !normalize the pivot row to the pivot value.  
    outMatrix(pivotR,:)=outMatrix(pivotR,:)/pivotVal !should normalize it. 
    !go to all rows below the pivot and reduce them by the pivot value.  
    do i=1,n,1
        if (i/=pivotR) then
            curVal=outMatrix(i,j) 
            outMatrix(i,:) = outMatrix(i,:)-(outMatrix(pivotR,:)*curVal)
        end if
    end do
    !so far so good.  
    pivotR = pivotR+1 !increase by 1 at the end, always.  
    if (pivotR>n) then
        exit
    end if
end do
!Now gets rref.  Correct -0 values 
do i=1,n,1
    do j=1,m,1
        if (outMatrix(i,j)==-0.0) then 
            outMatrix(i,j)=0.0
        end if
    end do
end do
deallocate(curRow)
deallocate(nextRow)
end function RREF
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
function inverse(matrix,n,m) RESULT(invMat)
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
!Now that both U and L are inversed as LMatrix and UMatrix, the actual matrix's inverse can be found by multiplying the two together. 
!A=LU so Ainv=U_inv*L_inv 
invMat = matmul(UMatrix,LMatrix)
end function inverse

!Linear Independence Checker 
integer function linearIndependence(matrix,n,m) RESULT(output)
integer::n,m,i,j,zeroFlag,countPivots
real,dimension(n,m)::matrix,retMatrix
real::currVal
countPivots=0
zeroFlag = 0
!!Use Gaussian Elim since its the quickest. Output is a 0 or a 1; 1 being linearly independent, 0 being linearly dependent. 
!call printMatrix(retMatrix,n,m)
!!Take gaussian elimination and then count rows which consist of only zeroes.  
retMatrix = GaussJordan(matrix,n,m)
!call printMatrix(retMatrix,n,m)
!we know that there are m columns in the matrix.  
!Now comb through and make sure that every row posesses a leading/pivot value on the diagonal.  If it does not, raise the zero flag and immediately return a 0.  
!Every column must have a pivot value.  So keep a count of the pivot values.  
!find first nonzero entry in the row.  Get the column it corresponds to.  
do i=1,n,1
    do j=1,m,1
        currVal=retMatrix(i,j)
        if (currVal/=0.0) then 
            countPivots=countPivots+1
            exit
        end if
    end do
end do
if (countPivots<m) then 
    output=0
    !write(*,*)'LINEARLY DEPENDENT'
else 
    output=1
    !write(*,*)'LINEARLY INDEPENDENT'
end if
end function linearIndependence

!Rank does Gaussian Elim, then checks to see how many nonzero rows exist.  Rank is the count of such vectors. 
subroutine QR(matrix,n,m,Q,R)
integer::n,m,isLinInd
real,dimension(n,m)::matrix
real,dimension(m,m),intent(out)::Q
real,dimension(n,m),intent(out)::R
real,dimension(n)::v !used for the QR decomposition portion of this
Q=0!Sets all values of the matrix to zero.  
R=0
!!Check first to see if the given matrix is linearly independent or not
isLinInd = linearIndependence(matrix,n,m)
!exit subroutine fail if the matrix is lin dep. 
if (isLinInd==0) then
    write(*,*)'FAILURE: Cannot perform Graham Schmidt QR Decomp if matrix is linearly dependent.'
    return
end if
!selectedCol=matrix(:,n) !should get last column of values
!Modified Graham Schmidt Method - more generally applicable. 
do j=1,m,1
    v = matrix(:,j)
    do i=1,j,1
        R(i,j)=dot_product(Q(:,i),matrix(:,j))
        v = v-R(i,j)*Q(:,i)
    end do
    R(j,j) = Norm(v)
    Q(:,j) = v/(R(j,j))
end do
end subroutine QR !Must be subroutine to return both Q and R.  

!Create function for the outer Product of two vectors; useful in Hesseberg stuff.  
function outerProduct(v1,v2) RESULT(outMatrix)
implicit none
integer::vecLen1,vecLen2
real, intent(in)::v1(:),v2(:)!NOTE: This is how you make the function automatically allocate variable inputs
real, dimension(:,:),allocatable::outMatrix,v1R,v2T
vecLen1 = size(v1)
vecLen2 = size(v2)
!Allocate memory to fit the given sizes.  !Both are represented as columns with form (:,1), with form (1,:).  Pass both as columns and then transpose. 
allocate(outMatrix(vecLen1,vecLen2))!Allocates for the outmatrix
allocate(v2T(1,vecLen2))
allocate(v1R(vecLen1,1))
v2T = reshape(v2, shape=[1,vecLen2])!Performs a transpose.  
v1R = reshape(v1, shape=[vecLen1,1])
!now perform the multiplication
!call printMatrix(v1,vecLen1,1)
!call printMatrix(v2T,1,vecLen2)
outMatrix = matmul(v1R,v2T)
!call printMatrix(outMatrix,vecLen1,vecLen2)
deallocate(v1R,v2T)
end function outerProduct

!Kronecker Product (Outer product for two matrices)
function kroneckerProduct(m1,m2) RESULT(outMatrix)
implicit none
integer::m1Cols,m1Rows,m2Cols,m2Rows !For matrix heights and widths; necessary in this.
integer::i,j,k,q
real, intent(in)::m1(:,:),m2(:,:)
real, dimension(:,:),allocatable::outMatrix
m1Cols = size(m1, dim=2) !#columns
m1Rows = size(m1, dim=1) !#rows
m2Cols = size(m2,dim=2)
m2Rows = size(m2,dim=1)
allocate(outMatrix(m1Rows*m2Rows,m1Cols*m2Cols)) !Allocates towards size of output.  
do i=1,m1Rows,1
    do j=1,m1Cols,1
        do k=1,m2Rows,1
            do q=1,m2Cols,1
                outMatrix((i*m2Rows+k),(j*m2Cols+q)) = m1(i,j)*m2(k,q)
            end do
        end do
    end do
end do
end function kroneckerProduct

!Hessenberg now works as intended.  
function hessenberg(inMatrix) RESULT(outMatrix)
real, intent(in)::inMatrix(:,:)
real, dimension(:,:),allocatable::outMatrix
real, dimension(:),allocatable::e1,v,x
integer::n,m,i
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
allocate(outMatrix(n,m))
allocate(x(m))
allocate(v(m))
outMatrix = inMatrix !assign in matrix to out matrix.  
!Now translate the python code to do Hessenberg reduction.  
do i=1,m-2,1
    !Check if already preallocated
    if (allocated(e1)) then 
        deallocate(e1)
    end if
    x = outMatrix(i+1:m,i) !x is supposed to be a column vector.  !WORKS
    allocate(e1(size(x)),source=0.0) !Allocate and fill e1 to mimic size of x. 
    e1(1)=1.0 !error here, calling 0 index when should be 1.
    v=sgn(x(0))*Norm(x)*e1+x
    v=v/Norm(v)
    !These two are signaling segmentation fault.  
    !first one has to start at second array position while the 
    !write(*,*)matmul(reshape(v,[1,size(v)]),outMatrix(i+1:m,i:m))!this works to flatten into 1D array as req'd by outer product
    outMatrix(i+1:m,i:m)=outMatrix(i+1:m,i:m) - 2.0*outerProduct(v,matmul(v,outMatrix(i+1:m,i:m)))
    outMatrix(1:m,i+1:m)=outMatrix(1:m,i+1:m) - 2.0*outerProduct(matmul(outMatrix(1:m,i+1:m),v),v)
end do
outMatrix = roundSmalls(outMatrix)
!call printMatrix(outMatrix,n,m)
end function hessenberg

!Updated eigenvalues function with Hessenberg Reduction and Schur stuff.  
function eigenvals(inMatrix, iterations) RESULT(eigens)
implicit none
real, intent(in)::inMatrix(:,:)
integer, intent(in), optional::iterations
integer::n,m,actIters,i
real::s, tol
real, dimension(:),allocatable::eigens
real, dimension(:,:),allocatable::eigenMatrix,shift, Q, R
!Check if iterations input is present, if not, default to 1000.  
tol = 1.00e-5
if (present(iterations)) then 
    actIters = iterations
else
    actIters = 100
end if
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
!Allocate the vars 
allocate(eigens(m))
allocate(eigenMatrix(n,m))
allocate(shift(n,m))
allocate(R(n,m))
allocate(Q(m,m))
!!!CREATE CONDITION CHECK FOR MATRIX SIZE, CHOOSE WISELY.
!if matrix size is less than 4 go with reg qr algorithm, else go with the shifted version.  Use Hessenberg for both regardless.  
!!convert to upper hessenberg form.
eigenMatrix = inMatrix  
eigenMatrix = Hessenberg(eigenMatrix)
if (n<5) then 
    !Add previous code here. 
    do i=1,actIters,1
    !these two work well enough with 1000 iters for 3x3, 4x4 as well.  Tolerance very good.  
    call QR(eigenMatrix,n,m,Q,R)
    eigenMatrix= matmul(R,Q)
    end do 
else
!!Now do the algorithm.  Should be more accurate and converge quicker.  
    do i=1,actIters,1
        !s is iterative shift
        s = eigenMatrix(n,n)
        shift = s*eye(n,n)
        !call printMatrix(shift,n,n)
        call QR(eigenMatrix-shift,n,m,Q,R) !assigns both Q and R

        eigenMatrix = matmul(R,Q)+shift !Check if there's one for matrix addition.  SHIFTS AREN'T ADDING UP.  
    end do
end if 
!Now isolate and sort the eigenvals.  
!Eigenvals are going to lie along the diagonal, so grab them like before.  
!use array eigens.  
do i=1,n,1
    if (abs(eigenMatrix(i,i))<tol) then 
        eigens(i)=0.0
    else
        eigens(i)=eigenMatrix(i,i)
    end if
end do
!then sort them.  
eigens = simpleSort(eigens)
end function eigenvals
!! Pseudo-Inverse() should return the pseudo-inverse of a given matrix.  Any size.  Function should work.  

!also make Span(), Rank()
!!Subroutine for SVD USE JACOBI.  
subroutine SVD(inMatrix,U,S,VT) 
implicit none
real, intent(in)::inMatrix(:,:)
integer::n,m, count, sweep, sweepMax, i, j, k, sorted, orthog, noisyA, noisyB
real::DBLeps, tol, p1, q1, a1, b1, sine, cosine, v, aerrorA, aerrorB, normC, prevNorm, aij, aik, qij, qik
real, intent(out), dimension(:,:)::U,VT
real, intent(out), dimension(:)::S
real,dimension(:,:), allocatable::A,Q
real,dimension(:), allocatable::t, cj, ck, col
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
allocate(A(n,m)) !copy of inMatrix
allocate(Q(m,m)) !copy V
allocate(t(m), source=0.0) !copy S !the zeroes thing.  
DBLeps = 1.0e-15
A = inMatrix
Q = eye(m,m)
!counters & setup
count=1
sweep=0
sweepMax = max(5*m,12)
tol = 10*n*DBLeps

!fill column vector values
do j=1,m,1
    t(j) = Norm(A(:,j)) !t should be a collection of the norms of all of the column vectors.  
end do

!continue from here. try to use more .gt. type ones as this looks cooler.  
do while((count.gt.0.0).AND.(sweep.le.sweepMax))
    count = (m*(m-1))/2 !rotation counter.  
    do j=1,(m-1),1
        do k=j+1,m,1
            !may have to write allocation procedure here so that they deallocate if already allocated, then reallocate every loop.  
            !Will need to be done incrementally.  Def them as dimension(:) types and only alloc/dealloc here.  a & b are reals.  
            if (allocated(cj)) then 
                deallocate(cj)
            end if
            if (allocated(ck)) then 
                deallocate(ck)
            end if
            allocate(cj(m))
            allocate(ck(m))
            cj=A(:,j)
            ck=A(:,k)
            p1=2*dot_product(cj,ck)
            a1=Norm(cj)
            b1=Norm(ck)
            !test for orthogonality or error in the cols.  
            aerrorA = t(j) !these two only need to be defined as reals since they're reading from the arrays.  
            aerrorB = t(k)
            q1=(a1*a1)-(b1*b1)

            v=(p1**2)+(q1**2)

            !If switches are here for the various values
            if (a1.ge.b1) then 
                sorted=1 !def sorted as an integer.  
            else
                sorted=0
            end if
            if (abs(p1).le.(tol*(a1*b1))) then
                orthog=1 !def Orthog auch Integer. 
            else
                orthog=0
            end if
            if (a1.lt.aerrorA) then
                noisyA=1
            else
                noisyA=0
            end if
            if (b1.lt.aerrorB) then 
                noisyB=1
            else
                noisyB=0
            end if 
            !end switches.  
            !check
            if (sorted==1.AND.(orthog==1.OR.noisyA==1.OR.noisyB==1)) then
                count = count-1
                continue
            end if
            !continue from here.  
            if (v==0.OR.sorted==0) then 
                cosine=0.0
                sine=1.0
            else
                cosine=sqrt(((v+q1)/(2.0*v)))
                sine=(p1/(2.0*v*cosine))
            end if
            !Apply rotation to A
            do i=1,n,1
                Aik=A(i,k)
                Aij=A(i,j)
                A(i,j)=Aij*cosine+Aik*sine
                A(i,k)=-Aij*sine +Aik*cosine
            end do
            !update the singular values
            t(j) = abs(cosine)*aerrorA + abs(sine)*aerrorB
            t(k) = abs(sine)*aerrorA + abs(cosine)*aerrorB

            !Apply rotation to Q now
            do i=1,m,1
                Qij=Q(i,j)
                Qik=Q(i,k)
                Q(i,j)=(Qij*cosine)+(Qik*sine)
                Q(i,k)=(-Qij*sine)+(Qik*cosine)
            end do
        end do
    end do
    sweep = sweep+1
end do
!now compute the singular values
prevNorm = -1.0
do j=1,m,1
    col = A(:,j) !by ref
    normC = Norm(col)
    !det if the sing val is zero
    if ((normC==0.0).OR.prevNorm==0.0.OR.((j.gt.0.0).AND.(normC.le.(tol*prevNorm)))) then
        t(j)=0.0
        do i=1,size(col),1
            col(i)=0.0 !updates A indirectly
        end do
        prevNorm=0.0
    else
        t(j)= normC
        do i=1,size(col),1
            col(i)=(col(i)*(1.0/normC))
        end do
        prevNorm = normC
    end if
end do
if (count.gt.0.0) then 
    write(*,*)'DOES NOT CONVERGE WITH JACOBI ITERS'
end if
U=A
S=t
Vt = transpose(Q) 

if(n.lt.m) then 
    U=U(:,0:n)
    S=t(0:n)
    Vt=Vt(0:n,:)
end if
!U, S, & Vh are automatically assigned/returned. 
write(*,*)
write(*,*)'U Matrix:'
call printMatrix(U,m,m)
write(*,*)'S Matrix:'
call printMatrix(S,m,1)
write(*,*)'V Transpose Matrix:'
call printMatrix(Vt,n,n)
end subroutine 
!Truncated SVD might be more interesting to look into.  

subroutine SVD2(inMatrix,U,S,VT) 
implicit none
real, intent(in)::inMatrix(:,:)
real, intent(out), dimension(:,:)::U,VT,S
integer::n,m
real,dimension(:,:), allocatable::A, AT, AAT, ATA, eigenVecs
real,dimension(:), allocatable::eigensAAT,eigensATA,singVals
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
allocate(A(n,m))
allocate(AT(m,n))
allocate(AAT(n,n))
allocate(ATA(m,m))
allocate(eigensAAT(n))
allocate(singVals(n))
allocate(eigenVecs(n,m))
!Manually calculate all of the SVD's since the fucking Jacobi won't work.  
A = inMatrix
AT = transpose(A)
!Need to first calculate AAT
AAT = matmul(A,AT)
ATA = matmul(AT,A)
!take the eigens of AAT
eigensAAT = eigenvals(AAT)
singVals = sqrt(eigensAAT)
!call printMatrix(eigensAAT,n,1) !for checking purposes.  
!call printMatrix(singVals,n,1)
!now that the eigens of AAT have been retrieved.  
!eigensATA = eigenvals(ATA)
!Try to round the values to the closest
!write(*,*)'ATA Eigenvalues',eigensATA
eigenvecs = eigenvectors(ATA,eigensAAT)
end subroutine 

!!Might also be good to do alternate *forward solving* version using other website.  
!create an eigenvectors function; returns the matrix of the eigenvectors.  
function eigenvectors(inMatrix,eigens) RESULT(eigenMatrix)
implicit none
real, intent(in)::inMatrix(:,:) !for the input
real, intent(in)::eigens(:) !for all eigenvals
integer::n,m,i,j,k,prevPivRow,ticker
real::curVal
real, dimension(:,:), allocatable::A,AR,eigenMatrix
real, dimension(:), allocatable::eigenVector !for the output.  
integer,dimension(:),allocatable::freeCols
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
k = size(eigens)
allocate(eigenMatrix(n,m))
allocate(AR(n,m))
allocate(A(n,m+1),source=0.0)
allocate(eigenVector(m)) !use for eigenvector column.  
allocate(freeCols(m))
!give us A with an extra row of zeroes to solve for.  
!make an eye matrix

!Now put this into a do loop where k=1,size(eigens),1
do k=1,size(eigens),1
    AR = inMatrix-eigens(k)*eye(n,m) !swap eigens in position i.
    !write(*,*)'This part of the code is running!'
    !Here is where it can be functionalized.  
    AR = RREF(AR)
    A=A+AR
    A(:,m+1) = 0.0 !set the last col to zero.
    !keep track of the pivot column indices.  
    prevPivRow=0
    !Should only go across the existing matrix, last column not to be measured AT ALL.  
    do j=1,m,1
        do i=1,n,1
            if (prevPivRow/=i.AND.A(i,j)==1.0) then 
                !check and allocate then reallocate
                !we know that i & j are the positions o
                prevPivRow=prevPivRow+1
            else
                curVal=A(i,j)
                A(i,m+1)=A(i,m+1)-curVal
                A(i,j)=A(i,j)-curVal
            end if
        !mark the column as 'free variable' and move to the far right (solution bracket) by subtracting that value. 
        !and if the previous Pivot row is the same as i and the value is nonzero, then label it as a free variable and subtract to the left.  
        !now assigns free variables.  
        !since j is now a column of free vals, we can just assign it to be equal to the first eigenvalue.  A character should be used for this.  
        end do
    end do
    !Now use eigenVector array
    !call printMatrix(A,n,m+1)
    ticker=0
    do j=1,m,1
        if(all(A(:,j)==0.0)) then 
            freeCols(ticker+1) = j !stores indices of free variables.  
            ticker=ticker+1
        else
            !should be all 1's at this point assuming RREF
            eigenVector(j) = A(j,m+1)
        end if
    end do
    freeCols(ticker+1:m)=0.0

    !now that all free variables are stored in memory, store their place in the resulting eigenvector as 1.  
    do j=1,m,1
        if(freeCols(j)/=0.0) then 
            eigenVector(freeCols(j)) = 1.0
        end if
    end do
    !write(*,*)eigenVector
    !write(*,*)
    eigenVector = eigenVector/Norm(eigenVector) !normalize eigenvectors.  
    eigenMatrix(k,:) = eigenVector
    A = 0.0 !reset A, otherwise residuals could cause all sorts of problems.  !oohhhh was only resetting the last one duh.  
    eigenVector=0.0
end do
!write(*,*)'Normalized:',eigenVector
!eigenvectors are now normalized. 
eigenMatrix(n,:)=0.0

eigenVector = eigenVector/Norm(eigenVector)
eigenMatrix(n,:) = eigenVector

!!!Basically need to figure this part out depending on how many eigenvectors there are present within the whole thing.  

!Flow is basically this: 
!Input Matrix, Eigenvalues=> Eigenvectors determined with regular system->If there are less eigenvalues than the num of variables, follow
!the same process to get them using previous solution->Output matrix of eigenvectors, which then are used for SVD.  

!It might be prudent to incorporate the solution mechanism as its own function/subroutine.  Like maybe have it return the solution column for a 
!matrix which might reduce otherwise.  maybe solve RREF?  BC Some RREF reduces to just a free one.  Maybe call the function 

!function solveSystem(inMatrix) RESULT(solutions)
!real,intent(in)::inMatrix(:,:)]
!integer::n,m
!real,dimension(:),allocatable::output
!n = size(inMatrix,dim=1) !#rows
!m = size(inMatrix,dim=2) !#cols
eigenVector = solveSystem(eigenMatrix(1:2,:))
eigenVector = eigenVector/Norm(eigenVector)
eigenMatrix(n,:) = eigenVector

!call printMatrix(A,n,m+1)
eigenMatrix = transpose(eigenMatrix) !put into form where the eigenvectors are vertical column vectors.  

!Now check if there are more eigenvectors than eigenvalues, if so, an extra step must be done to figure out the last one.  
write(*,*)
write(*,*)'Matrix of Eigenvectors (columns):'
call printMatrix(eigenMatrix,n,m) !check for final printing.  !technically already transposed.  

end function eigenVectors

!Function for solving the system of equations (used in eigenvectors function)
function solveSystem(inMatrix) RESULT(solutions)
implicit none
real,intent(in)::inMatrix(:,:)
integer::n,m,i,j,prevPivRow,ticker
real::curVal
integer,dimension(:),allocatable::freeCols
real,dimension(:),allocatable::solutions
real,dimension(:,:),allocatable::A
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
allocate(A(n,m+1))
allocate(freeCols(m))
allocate(solutions(m))
A=0.0!reset again
A=A+inMatrix
A(:,m+1)=0.0
!begin here.  
A = RREF(A)
prevPivRow=0
do j=1,m,1
    do i=1,n,1
        if (prevPivRow/=i.AND.A(i,j)==1.0) then 
            prevPivRow=prevPivRow+1
        else
            curVal=A(i,j)
            A(i,m+1)=A(i,m+1)-curVal
            A(i,j)=A(i,j)-curVal
        end if       
    end do
end do
freeCols=0
ticker=0
do j=1,m,1
    if(all(A(:,j)==0.0)) then 
        freeCols(ticker+1) = j !stores indices of free variables.  
        ticker=ticker+1
    else
        !should be all 1's at this point assuming RREF
        solutions(j) = A(j,m+1)
    end if
end do
freeCols(ticker+1:m)=0.0
do j=1,m,1
    if(freeCols(j)/=0.0) then 
        solutions(freeCols(j)) = 1.0
    end if
end do
deallocate(A)
deallocate(freeCols)
end function solveSystem



end program matrixOperations
