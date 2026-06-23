module libdonmat
    implicit none
contains
    function readMatrix(filename,n,m) RESULT(matrix)
    implicit none
    integer,intent(in)::n,m
    character(len=*),intent(in)::filename
    real,dimension(:,:),allocatable::matrix
    integer::i,j
    allocate(matrix(n,m))
    open(20,file=filename, status='old') !process still works for matrix4.  
    !! Create file-read method for this.  Basically just from one input to another, or allow multi-input?  File works better I think.  
    do i=1,n,1
        do j=1,m,1
            !!!loop through and fill matrix values in.  
            read(20,*)matrix(i,j)
        end do
    end do
    close(20)
    end function readMatrix

    !Define print matrix subroutine; inputs are the matrix and number of rows.  Rets matrix as-is. 
    subroutine printMatrix(inMatrix)
    implicit none
    real,intent(in)::inMatrix(:,:)
    integer::n,k
    n = size(inMatrix,dim=1) !#rows
    do k = 1, n
        write(*,*)inMatrix(k,:)
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
            t=a/b
            s=1/sqrt(1+t**2)
            c=s*t
        else
            t=b/a
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
    integer::actRightShift,actDownShift,i,j
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
    j=1
    do i=1+actDownShift,n,1
        if (j.ge.1.AND.j.le.m) then
            eyeMatrix(i,j+actRightShift) = 1.0
            j=j+1
        end if
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
    !Gaussian Elimination goes here.  
    function Gaussian(inMatrix) RESULT(outMatrix)
    implicit none
    real, intent(in)::inMatrix(:,:)
    integer::n,m,i,j,pivotR,nextPivotR
    real::pivotVal,curVal
    real,dimension(:),allocatable::curRow,nextRow
    real,dimension(:,:),allocatable::outMatrix
    n = size(inMatrix,dim=1) !#rows
    m = size(inMatrix,dim=2) !#cols
    allocate(outMatrix(n,m))
    allocate(curRow(m))
    allocate(nextRow(m))
    pivotR=1
    nextPivotR=1
    outMatrix = inMatrix
    !Basically copy over the logic from RREF BUT only allow it to go downwards instead of both directions from the current pivot.
    !Next idea, to help eliminate in the solveSystems function, the pivot values should be maintained at their original values after
    !complete multiplication.  
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
        if(pivotR.lt.m) then
        !now that we have the pivot row
            pivotVal=outMatrix(pivotR,j)
            !normalize the pivot row to the pivot value.  
            outMatrix(pivotR,:)=outMatrix(pivotR,:)/pivotVal !should normalize it. 
        end if
        !go to all rows below the pivot and reduce them by the pivot value.  
        do i=pivotR+1,n,1
            curVal=outMatrix(i,j) 
            outMatrix(i,:) = outMatrix(i,:)-(outMatrix(pivotR,:)*curVal)
        end do
        !now multiply row back by its pivotVal ONLY if it's less than m, as the mth or the last pivot rather will be the final value.  
        if (pivotR.lt.m) then 
            outMatrix(pivotR,:)=outMatrix(pivotR,:)*pivotVal
        end if
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
    end function Gaussian

    !RREF rewritten to be more universal with pivots.  
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
    !!!Function for determinant !needs a rework. 
    real function det(inMatrix) RESULT(r)
    implicit none
    real, intent(in)::inMatrix(:,:)
    real,dimension(:),allocatable::selectedRow !allocate to m
    real,dimension(:,:),allocatable::LMatrix,workMat
    integer::i,j,n,m
    real::a,b,c,d,detL,detU, currVal, nextVal, inputVal  !return for determinant
    n = size(inMatrix, dim=1) !#rows
    m = size(inMatrix, dim=2) !#columns
    allocate(workMat(n,m))
    allocate(selectedRow(m))
    if(n/=m) then
        write(*,*)'Cannot take the determinant of a nonsquare matrix.'
        return
    end if
    workMat = inMatrix
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
    function inverse(inMatrix) RESULT(invMat)
    implicit none
    real,intent(in)::inMatrix(:,:)
    real,dimension(:,:),allocatable::LMatrix,UMatrix,invMat,iMatrix,iMatrix2 !invMat is return value
    real,dimension(:),allocatable::selectedRow,selectedRow2 !second specifically for U inverse.  
    integer::n,m,i,j
    real::currVal,nextVal,inputVal
    n = size(inMatrix, dim=1) !#rows
    m = size(inMatrix, dim=2) !#columns
    allocate(UMatrix(n,m)) !allocate memory for the input/used variables:
    allocate(invMat(n,m))
    allocate(iMatrix(n,m))
    allocate(iMatrix2(n,m))
    allocate(selectedRow(m))
    allocate(selectedRow2(m))
    allocate(LMatrix(n,m), source=0.0) !allocates the memory to it and fills all values with zeroes. 
    UMatrix = inMatrix !Assigns UMat
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

    !LU Decomposition
    subroutine LU(inMatrix,L,U)
    implicit none
    !variable declarations (Order is input/output, then declare allocatables, then declare reg types)
    real, intent(in)::inMatrix(:,:)
    real, dimension(:,:), allocatable, intent(out)::L,U
    real, dimension(:),allocatable::selectedRow
    integer::n,m,i,j
    real::currVal,nextVal,inputVal
    !Measurement for allocation
    n = size(inMatrix, dim=1) !#rows
    m = size(inMatrix, dim=2) !#columns
    !now that the size of the input variable has been defined, allocate the memory depending on the size of the matrix
    !allocation
    if (n/=m) then
        allocate(L(n,n))
        allocate(U(n,m))
    else
        allocate(L(n,m))
        allocate(U(n,m))
    end if
    allocate(selectedRow(m))
    !now that they have been properly allocated begin assignment
    U=inMatrix !assign inmatrix
    do i=1,n,1 !assign the diagonals to 1.0 to construct identity matrix
            L(i,i)=1.0
    end do
    !perform decomposition into upper and lower triangular matrices.  
    do i=1, m-1, 1
        selectedRow = U(i,:)
        currVal = U(i,i)
        do j=i, n-1, 1
            nextVal = U(j+1,i) !!Nextval will be input into its place in the identity matrix.  
            inputVal=(nextVal/currVal) !inputVal is to be entered into the precise place in the iMatrix
            L(j+1,i)=inputVal
            U(j+1,:)=U(j+1,:)-(selectedRow*inputVal)
        end do
    end do
    !L & U should be properly decomposed now.  
    end subroutine LU

    !Linear Independence Checker !all good now.  
    integer function linearIndependence(inMatrix) RESULT(output)
    real, intent(in)::inMatrix(:,:)
    real, dimension(:),allocatable::selectedRow
    real,dimension(:,:),allocatable::retMatrix !allocate to n,m
    integer::n,m,k,i,j,zeroFlag,countPivots, curPivRow
    real::currVal, nextVal
    n = size(inMatrix, dim=1) !#rows
    m = size(inMatrix, dim=2) !#columns
    allocate(retMatrix(n,m))
    countPivots=0
    curPivRow=0
    zeroFlag = 0
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
    !we know that there are m columns in the matrix.  
    !check all of the rows and the columns.  
    do i=1,n,1
        do j=1,m,1
            currVal=retMatrix(i,j)
            if (currVal/=0.0) then 
                countPivots=countPivots+1
                exit
            end if
        end do
    end do
    !write(*,*)'Count Pivots:',countPivots
    !call printMatrix(retMatrix)
    if (countPivots<m) then 
        output=0
        !write(*,*)'LINEARLY DEPENDENT'
    else 
        output=1
        !write(*,*)'LINEARLY INDEPENDENT'
    end if
    end function linearIndependence

    !Rank does Gaussian Elim, then checks to see how many nonzero rows exist.  Rank is the count of such vectors. 
    subroutine gramSchmidtQR(inMatrix,Q,R)
    real,intent(in)::inMatrix(:,:)
    real,dimension(:,:),allocatable,intent(out)::Q,R!n,n !n,m
    real,dimension(:),allocatable::v !used for the QR decomposition portion of this
    integer::n,m,i,j,isLinInd
    !take measurements
    n = size(inMatrix, dim=1) !#rows
    m = size(inMatrix, dim=2) !#columns
    allocate(Q(n,n))
    allocate(R(n,m))
    allocate(v(n))
    Q=0!Sets all values of the matrix to zero.  
    R=0
    !!Check first to see if the given matrix is linearly independent or not
    isLinInd = linearIndependence(inMatrix)
    !write(*,*)'Is lin Ind'
    !call printMatrix(inMatrix)
    !exit subroutine fail if the matrix is lin dep. 
    if (isLinInd==0) then
        write(*,*)'FAILURE: Cannot perform Graham Schmidt QR Decomp if matrix is linearly dependent.'
        return
    end if
    !selectedCol=matrix(:,n) !should get last column of values
    !Modified Graham Schmidt Method - more generally applicable. 
    do j=1,m,1
        v = inMatrix(:,j)
        do i=1,j,1
            R(i,j)=dot_product(Q(:,i),inMatrix(:,j))
            v = v-R(i,j)*Q(:,i)
        end do
        R(j,j) = Norm(v)
        Q(:,j) = v/(R(j,j))
    end do
    end subroutine gramSchmidtQR !Must be subroutine to return both Q and R.  

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
    outMatrix = matmul(v1R,v2T)
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
    deallocate(x)
    deallocate(v)
    end function hessenberg

    !! Pseudo-Inverse() should return the pseudo-inverse of a given matrix.  Any size.  Function should work.  
    !Attempt 2 at Householder QR
    subroutine QR(inMatrix,Q,R)
    implicit none
    real,intent(in)::inMatrix(:,:)
    real,dimension(:,:),allocatable,intent(out)::Q,R
    real,dimension(:),allocatable::w
    real,dimension(:,:),allocatable::tempQ
    integer::n,m,j
    real::normX, tau, u1,s
    n = size(inMatrix,dim=1) !#rows
    m = size(inMatrix,dim=2) !#cols
    !allocation step
    allocate(R(n,m))
    allocate(Q(n,n))
    R=inMatrix
    Q=eye(n,n)
    do j=1,m,1
        normX = Norm2(R(j:n,j))
        if(normX<1.00e-10) cycle
        s=sgn(R(j,j)) !s is just integer
        u1=R(j,j)+(s*normX)!u1 is a real
        if(abs(u1)<1.00e-10) cycle !this is what's causing the NaN's
        if(allocated(w)) then 
            deallocate(w)
        end if 
        allocate(w(n-j+1))
        w=R(j:n,j)/u1 !w is a column vector with length n-j
        w(1)=1.0 !w(1) is the first value
        tau=s*(u1/normX) !tau is just a real
        R(j:n,:)= R(j:n,:)-matmul(reshape(tau*w,[size(w),1]),(matmul(reshape(w,[1,size(w)]),R(j:n,:))))
        if(allocated(tempQ)) then 
            deallocate(tempQ)
        end if
        allocate(tempQ(size(Q,1),m-j+1))
        tempQ=matmul(Q(:,j:n),(reshape(w,[size(w),1]))) !def as temp
        Q(:,j:n)=Q(:,j:n)-matmul(tempQ,reshape(w*tau,[1,size(w)]))
    end do
    !now run roundsmalls on them
    Q=roundSmalls(Q)
    R=roundSmalls(R)
    end subroutine QR

    !Rewrite eigenValues 
    function eigenvalsHouseholder(inMatrix,iterations) RESULT(eigens)
    implicit none
    real, intent(in)::inMatrix(:,:)
    integer, intent(in), optional::iterations
    integer::n,m,actIters,i
    real::s, tol, trace,disc
    real, dimension(:),allocatable::eigens
    real, dimension(:,:),allocatable::eigenMatrix,shift, Q, R
    tol = 1.00e-5
    if (present(iterations)) then 
        actIters = iterations
    else
        actIters = 50
    end if
    n = size(inMatrix,dim=1) !#rows
    m = size(inMatrix,dim=2) !#cols
    !Allocate the vars 
    allocate(eigens(m)) !eigens being allocated to M.  This could be a problem as it doesn't check how many eigens there really are.  
    allocate(eigenMatrix(n,m))
    allocate(shift(n,m))
    allocate(R(n,m))
    allocate(Q(m,m))
    !begin
    eigenMatrix = inMatrix  
    !check size of the matrix input: !for smaller size, calculate the eigenvalues with the determinant.  Cannot have a 1x2 matrix, at that point
    !it's just a vector. 
    if (n==2) then
        !calculate the trace first (add all vals along the diagonal)
        trace=0.0
        do i=1,n,1
            trace=trace+eigenMatrix(i,i)
        end do
        !once trace is calculated, you can enter the equation for the eigenvals. 
        !calculate the discriminant
        disc=trace**2 - (4*det(eigenMatrix))
        eigens(1)=(trace+sqrt(disc))/2
        eigens(2)=(trace-sqrt(disc))/2
    end if


    if (n<50.AND.n>2) then 
        !Add previous code here. 
        do i=1,actIters,1
        !these two work well enough with 1000 iters for 3x3, 4x4 as well.  Tolerance very good.  
            call QR(eigenMatrix,Q,R)
            eigenMatrix= matmul(R,Q)
            Q=0 !reset values
            R=0 !reset values afterwards; junk values appear more frequently over time.  
        end do 
    else
    !!Now do the algorithm.  Should be more accurate and converge quicker.  
        do i=1,actIters,1
            !s is iterative shift
            s = eigenMatrix(n,n)
            shift = s*eye(n,n)
            call QR(eigenMatrix-shift,Q,R) !assigns both Q and R
            eigenMatrix = matmul(R,Q)+shift !Check if there's one for matrix addition.  SHIFTS AREN'T ADDING UP.  
        end do
    end if 
    !only fulfill this condition if larger than 2, otherwise skip right to sorting.  
    if (n.gt.2) then
        do i=1,n,1
            if (abs(eigenMatrix(i,i))<tol) then 
                eigens(i)=0.0
            else
                eigens(i)=eigenMatrix(i,i)
            end if
        end do
    end if
    !then sort them.  
    eigens = simpleSort(eigens)
    deallocate(eigenMatrix)
    deallocate(shift)
    deallocate(Q)
    deallocate(R)
    end function
    !also make Span(), Rank()
    !Truncated SVD might be more interesting to look into.  
    !SVD Function is here!!!
    subroutine SVD(inMatrix,U,S,VT) 
    implicit none
    real, intent(in)::inMatrix(:,:)
    real, intent(inout),dimension(:,:)::U,VT,S
    integer::n,m,i
    real,dimension(:,:), allocatable::A, AT, AAT, ATA, eigenVecs,identS,Uproto
    real,dimension(:), allocatable::eigensAAT,singVals
    n = size(inMatrix,dim=1) !#rows
    m = size(inMatrix,dim=2) !#cols
    allocate(A(n,m))
    allocate(AT(m,n))
    allocate(AAT(n,n))
    allocate(ATA(m,m))
    allocate(eigensAAT(m))
    allocate(singVals(m))
    allocate(eigenVecs(m,m))
    allocate(identS(n,m))
    allocate(Uproto(n,n))
    !Manually calculate all of the SVD's since the fucking Jacobi won't work.  
    A = inMatrix
    AT = transpose(A)
    !Need to first calculate AAT
    AAT = matmul(A,AT)
    ATA = matmul(AT,A)
    !take the eigens of AAT
    eigensAAT = eigenvalshouseholder(AAT)
    !write(*,*)eigensAAT
    singVals = sqrt(eigensAAT)
    !Assign the singVals to their own matrix.  
    identS = eye(n,m)
    do i=1,size(singVals),1 !tie this to the number of singular values instead of to the columns.  Previously had an indexing error.  
        identS(i,i)=identS(i,i)*singVals(i)
    end do
    !works right now
    eigenvecs = eigenvectors(ATA,eigensAAT) !gets the V matrix essentially.  
    VT=transpose(eigenvecs)!and this is officially VT.  
    !write(*,*)'V Transpose Matrix:'
    deallocate(AT)
    deallocate(AAT)
    deallocate(ATA)
    deallocate(eigensAAT)
    !Compute left singular vectors using VT, singVals, and the originally A (input Matrix)
    !THIS IS WHAT TO DO NEXT. 
    !leftsingular vectors 
    !u_i=1/(sv_i)*A*v_i
    !basically the ith vector of u is equal to 1/ith singVal multiplied by the matrix multiplication of the original matrix and the eigenvectors.  
    Uproto=0.0
    do i=1,n,1 
        Uproto(:,i) = matmul(A,eigenvecs(:,i))*(1/singVals(i))
    end do
    U=Uproto
    S = identS !Sigma Matrix is now in full effect.  Works amazingly!
    end subroutine SVD

    !!Might also be good to do alternate *forward solving* version using other website.  
    !create an eigenvectors function; returns the matrix of the eigenvectors.  
    function eigenvectors(inMatrix,eigens) RESULT(eigenMatrix)
    implicit none
    real, intent(in)::inMatrix(:,:) !for the input
    real, intent(in)::eigens(:) !for all eigenvals
    integer::n,m,i,k,numEigs,counter
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
    numEigs = size(eigens)
    !Now put this into a do loop where k=1,size(eigens),1
    do k=1,size(eigens),1
        AR = inMatrix-eigens(k)*eye(n,m) !swap eigens in position i.
        !write(*,*)'This part of the code is running!'
        !Here is where it can be functionalized.  
        eigenVector = solveSystem(AR)
        eigenVector = eigenVector/Norm(eigenVector) !normalize eigenvectors.  
        eigenMatrix(k,:) = eigenVector
        eigenVector=0.0
    end do
    !if the number of eigenvalues obtained is less than the number of columns in the matrix (also rows) then it knows to orthogonalize.  
    counter=0
    if (numEigs.lt.m) then
        do i=numEigs+1,n,1
            eigenVector = solveSystem(eigenMatrix(1:numEigs+counter,:))
            eigenVector = eigenVector/Norm(eigenVector)
            eigenMatrix(i,:) = eigenVector ! I believe this is where the problem is
            counter=counter+1
        end do
    end if
    eigenMatrix = transpose(eigenMatrix) !put into form where the eigenvectors are vertical column vectors.  
    deallocate(AR)
    deallocate(A)
    deallocate(eigenVector)
    deallocate(freeCols)
    end function eigenVectors

    !Function for solving the system of equations (used in eigenvectors function primarily)
    !Input is a matrix of any size and returns an array of solutions.  
    function solveSystem(inMatrix) RESULT(solutions)
    implicit none
    real,intent(in)::inMatrix(:,:)
    integer::n,m,i,j,prevPivRow,ticker,pivCount,pivCts
    real::curVal,pivAvg
    integer,dimension(:),allocatable::freeCols
    real,dimension(:),allocatable::solutions,pivArray
    real,dimension(:,:),allocatable::A
    n = size(inMatrix,dim=1) !#rows
    m = size(inMatrix,dim=2) !#cols
    allocate(A(n,m+1))
    allocate(freeCols(m))
    allocate(solutions(m))
    allocate(pivArray(m))
    A=0.0!reset again
    !call printMatrix(inMatrix,n,m) !note, try truncating the final value after performing Gauss-Jordan.  This should be a standard test in the 
    ! solution.  If the corner value is too small, it could affect calculations and cause NANs.  
    A=A+inMatrix
    A(:,m+1)=0.0
    !begin here.  
    A = Gaussian(A) !take the gauss-jordan elim of it, then, eliminate the values smaller than the tolerance value.  
    !check the pivots within the tolerance range
    !define tolerance as the rough average between the absolute values of the constituent pivots.
    !get the variance of each pivot from the average and then check to tell the magnitude of that difference?
    !find ones along the main diagonal? nth nonzero column entry?  so for like 
    prevPivRow=0
    pivCount=0
    pivArray=0.0
    do j=1,m,1
        do i=1,n,1
            if (prevPivRow.lt.i.AND.A(i,j)/=0.0) then 
                !what we would want to do here is disregard the elements which are nonpivotelements
                pivArray(prevPivRow+1)=A(i,j) !set this to be the main pivot element
                prevPivRow=prevPivRow+1
                pivCount=pivCount+1
                exit !effectively cycles current loop.  
            end if       
        end do
    end do
    !zero them out by taking the average, then subtracting the average from all of them, rounding to nearest int, and then subtracting from abs.
    !near-zero values will effectively become zero.  
    pivArray = abs(pivArray)
    pivAvg=sum(pivArray(1:pivCount))/pivCount 
    pivArray(1:pivCount)=pivArray(1:pivCount)-pivAvg
    !int ->round to zero.  nint-> round to nearest number.  
    pivArray=nint(pivArray)
    pivAvg=nint(pivAvg)
    pivArray=abs(pivArray)-pivAvg
    !now, the incredibly small values should be zeroed out. 
    !we know that the values in the matrix will shrink over time so that all of the near null values will be on the bottom. Also for array.  
    pivCts=0 
    do i=1,pivCount,1
        if (pivArray(i)==0.0) then 
            exit
        end if
        pivCts=pivCts+1
    end do
    !piv cts should be the last point before zero
    do i=pivCts+1,n,1
        A(i,:) = 0.0
    end do
    A=RREF(A)
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

    !Function that gives the determinant for nonsquare matrices (where m<n)
    function gramDeterminant(inMatrix) RESULT(gramDet)
    implicit none
    real,intent(in)::inMatrix(:,:)
    real,dimension(:,:),allocatable::AT,ATA
    integer::n,m
    real::gramDet
    n = size(inMatrix, dim=1) !#rows
    m = size(inMatrix, dim=2) !#columns
    allocate(AT(m,n))
    allocate(ATA(m,m))
    !check for condition
    if(m.gt.n) then 
        gramDet=0.0
        return
    end if
    AT=transpose(inMatrix)
    ATA=matmul(AT,inMatrix)
    gramDet=det(ATA)
    gramDet=sqrt(gramDet)!set it to its sqrt
    deallocate(AT)
    deallocate(ATA)
    end function gramDeterminant
    !Next step after SVD is to create Moore-Penrose Pseudoinverse, which should be the last NEW thing for a long minute.  
    function pseudoinverse(inMatrix) RESULT(invMatrix)
    implicit none
    real,intent(in)::inMatrix(:,:)
    real,dimension(:,:),allocatable::U,S,VT,invMatrix
    integer::n,m,i,j
    n = size(inMatrix, dim=1) !#rows
    m = size(inMatrix, dim=2) !#columns
    allocate(U(n,n))
    allocate(S(n,m))
    allocate(VT(m,m))
    allocate(invMatrix(m,n)) !specifically for output
    call SVD(inMatrix,U,S,VT) !fills and assigns values.  
    write(*,*)'Sigma Matrix'
    call printMatrix(S)
    do i=1,n,1
        do j=1,m,1
            if(S(i,j).gt.0.0) then 
                S(i,j)=1.0/S(i,j)
            end if
        end do
    end do
    !So now what I can do is calculate exactly how to multiply them together.
    invMatrix = matmul(transpose(VT),matmul(transpose(S),transpose(U)))
    deallocate(U)
    deallocate(S)
    deallocate(VT)
    end function pseudoinverse
end module libdonmat
