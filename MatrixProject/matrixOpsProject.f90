program matrixOperations
implicit none
!! Takes file as input for matrix values, has output file as well. Creates an allocatable matrix determined by user input.  
!! Should determine Reduced Row Echelon Form using a radix sort of values
!! NOTE: MAKE SURE THAT THE FILE IS WRITTEN IN ORDER: FIRST ROW->NthRow n=noRows, m = noCols
integer::n,m,i,j,retty !Assigns values of row length and column height
real,dimension(:,:),allocatable::mymatrix,returnMatrix,inv,Q,R,reye,outerprod
real,dimension(:),allocatable::vec,eigs,vec2,vec3
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

returnMatrix = RREF(mymatrix,n,m)
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

write(*,*)'The Eigenvalues of the Matrix are:'
eigs = eigenvals(mymatrix,n,m)
call printMatrix(eigs,1,m)
!reye = eye(n,m)

!!First vector of the matrix
vec = mymatrix(1,:)
normal = Norm(vec,m)
write(*,*)'Normal of the first vector:',normal

vec2 = mymatrix(:,1)
vec3 = myMatrix(:,2)
outerprod = outerProduct(vec2,vec3)

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

!!Function Takes the Norm
real function Norm(vector,n) RESULT(vectNormal)
integer::n,i
real, dimension(n)::vector
real::sum
sum=0.0
!Take square root of the squared sum.  
do i=1,n,1
    sum=sum+(vector(i)**2)
end do 
vectNormal = sqrt(sum)
end function Norm

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
function RREF(matrix, n, m) RESULT(retMatrix)
implicit none
integer::n,m,i,k,g,zeroFlag !g is the row index of the last value of the RREF
real::currVal,nextVal,normal
real,dimension(n,m)::matrix, retMatrix
real,dimension(m)::selectedRow
retMatrix = matrix !swaps returnMatrix for the original input
zeroFlag=0
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
!take the norm of ANY of the rows; if zero then set zero flag.  Assign g to all row ind
do i=1,n,1
    selectedRow = retMatrix(i,:)
    normal = Norm(selectedRow, m)
    if (normal==0) then
        g=i-1 !g is the index of the row where the first zero row is.  It sets g=i-1. 
        zeroFlag=1
    end if
end do

!Problem is here.  
!Check to see if a zero row has been detected before proceeding.  
if(zeroFlag==1) then
    do k=g-1, 1, -1
        do i=k-1, 1, -1
            !write(*,*)retMatrix(i,k)
            selectedRow=retMatrix(k,:)
            nextVal = retMatrix(i,k)
            retMatrix(i,:)=retMatrix(i,:)-(selectedRow*nextVal)
        end do
    end do

!add new small loop if the matrix is square to clean up final column values
    if (n==m) then
        retMatrix(g,:)=retMatrix(g,:)/retMatrix(g,g) !sets it to 1 if 1
        selectedRow = retMatrix(g,:)
        do i=g-1,1,-1
            nextVal = retMatrix(i,g)
            retMatrix(i,:)=retMatrix(i,:)-(selectedRow*nextVal)
        end do
    end if
else
    do k=m-1, 1, -1
        do i=k-1, 1, -1
        !write(*,*)retMatrix(i,k)
            selectedRow=retMatrix(k,:)
            nextVal = retMatrix(i,k)
            retMatrix(i,:)=retMatrix(i,:)-(selectedRow*nextVal)
        end do
    end do
    if (n==m) then
        retMatrix(n,:)=retMatrix(n,:)/retMatrix(n,n) !sets it to 1 if 1
        selectedRow = retMatrix(n,:)
        do i=n-1,1,-1
            nextVal = retMatrix(i,n)
            retMatrix(i,:)=retMatrix(i,:)-(selectedRow*nextVal)
        end do
    end if
end if
!Clean up negative zero values; could be a problem later on.  s
do i=1, n, 1
    do k=1, m ,1
        if (abs(retMatrix(i,k))==0) then
            retMatrix(i,k)=0
        end if
    end do
end do
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
    R(j,j) = Norm(v,n)
    Q(:,j) = v/(R(j,j))
end do
end subroutine QR !Must be subroutine to return both Q and R.  

!! Pseudo-Inverse() should return the pseudo-inverse of a given matrix.  Any size.  Function should work.  
!! Method: Graham Schmidt for QR Decomposition

!!!Eigenvalues and shittt.  
!!For storing the values of the eigenstuff
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
!call printMatrix(eyeMatrix,n,m)
end function

!Take the eigen QR of a matrix
function eigenvals(inMatrix,n,m,iterations) RESULT(eigens)
implicit none
integer,intent(in)::n,m
integer,intent(in),optional::iterations
real, dimension(n,m), intent(in)::inMatrix
integer::i,actIters
real, dimension(n,m)::eigenMatrix,R
real, dimension(m,m)::Q
real, dimension(m)::eigens
!!!Check input options to make sure that they're properly there.  
if(present(iterations)) then
    actIters=iterations
else
    actIters = 100000
end if
eigenMatrix=inMatrix !assign copy.   
!now the actual algorithm implementation
do i=1,actIters,1
    !these two work well enough with 1000 iters for 3x3, 4x4 as well.  Tolerance very good.  
    call QR(eigenMatrix,n,m,Q,R)
    eigenMatrix= matmul(R,Q)
end do
!Isolate eigenvalues from main diagonal.  CURRENTLY UNSORTED
do i=1,m,1
    eigens(i)=eigenMatrix(i,i)
end do
end function eigenvals

!Add a sorting algorithm of some variety for the eigenvalues
!also make Span(), Rank()
!!Subroutine for SVD
subroutine SVD(inMatrix, n, m)
implicit none
integer::n,m
real, dimension(n,m)::inMatrix
end subroutine 


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

function hessenberg(inMatrix) RESULT(outMatrix)
real, intent(in)::inMatrix(:,:)
real, dimension(:,:),allocatable::outMatrix
real, dimension(:),allocatable::x
integer::n,m,i,j
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols
allocate(outMatrix(n,m))
allocate(x(m))
!Now translate the python code to do Hessenberg reduction.  
do i=1,m-2,1
    x = inMatrix(i+1:m,i)
    !Continue from here.  
end do

end function hessenberg

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

end program matrixOperations
