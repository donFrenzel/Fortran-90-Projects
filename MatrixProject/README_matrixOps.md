**Feature List**

*Matrix Functions:*
- RREF -> RREF(inputMatrix)
- Gauss-Jordan Elimination -> GaussJordan(inputMatrix)
- Determinant (only for Square for now) -> det(matrix,n,m) *n,m autodetection to come later
- Inverse -> inverse(matrix,n,m) 
- Outer Product -> outerProduct(vector1,vector2)
- Kronecker Product ->kroneckerProduct(matrix1,matrix2)
- Linear Independence Check -> linearIndependence(matrix,n,m) *n,m autodetection to come later.  
- Hessenberg Reduction -> hessenberg(inputMatrix)
- RoundSmalls (rounds to zero for small vals in a matrix) -> roundSmalls(inputMatrix, precision=(DEFAULT=0.00001))
- Eigenvalues (QR Algorithm & QR with shifts) -> eigenvals(inputMatrix, iterations=(DEFAULT=100))

*Matrix Subroutines:
- LU Decomposition (It's there but it's rudimentary - still needs implementation on its own as a subroutine)
- QR Decomposition -> QR(inputMatrix,n,m,Q,R) *Q & R must be predefined for subroutine output.  n,m autodetection to come later)
- SVD (limited functionality) -> SVD(inputMatrix,U,S,VT) *U,S,VT must be predefined in program for subroutine output.
- Givens -> givens(a,b,c,s) *c,s must be predefined in program for subroutine output.
- Print Matrix -> printMatrix(inputMatrix,n,m) *n,m autodetection to come later.  



*General Functions (misc.):*
- Norm -> Norm(vector)
- Sign -> sgn(val)
- Eye (equiv. to numpy python function, allows for custom identity matrix construction) -> eye(n,m,rightShift=(DEFAULT 0),downShift=(DEFAULT 0))
- simpleSort (sorting algorithm used for eigens but is general purpose) -> simpleSort(array)
