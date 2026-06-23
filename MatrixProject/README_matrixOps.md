**Feature List**

*Matrix Functions:*
- RREF -> RREF(inputMatrix)
- Gauss-Jordan Elimination -> GaussJordan(inputMatrix)
- Determinant (only for Square for now) -> det(matrix)
- Gram Determinant -> gramDet(matrix)
- Inverse -> inverse(matrix) 
- Outer Product -> outerProduct(vector1,vector2)
- Kronecker Product ->kroneckerProduct(matrix1,matrix2)
- Linear Independence Check -> linearIndependence(matrix)
- Hessenberg Reduction -> hessenberg(inputMatrix)
- RoundSmalls (rounds to zero for small vals in a matrix) -> roundSmalls(inputMatrix, precision=(DEFAULT=0.00001))
- Eigenvalues (QR Algorithm & QR with shifts) -> eigenvals(inputMatrix, iterations=(DEFAULT=100))
- Eigenvectors -> eigenvectors(inputMatrix) 
- Moore-Penrose Pseudo-Inverse -> pseudoinverse(inputMatrix)

*Matrix Subroutines:
- LU Decomposition (It's there but it's rudimentary - still needs implementation on its own as a subroutine)
- QR Decomposition -> QR(inputMatrix,Q,R) *Q & R must be predefined for subroutine output
- Singular Value Decomposition -> SVD(inputMatrix,U,S,VT) *U,S,VT must be predefined in program for subroutine output.
- Givens -> givens(a,b,c,s) *c,s must be predefined in program for subroutine output.
- Print Matrix -> printMatrix(inputMatrix)



*General Functions (misc.):*
- Norm -> Norm(vector)
- Sign -> sgn(val)
- Eye (equiv. to numpy python function, allows for custom identity matrix construction) -> eye(n,m,rightShift=(DEFAULT 0),downShift=(DEFAULT 0))
- simpleSort (sorting algorithm used for eigens but is general purpose) -> simpleSort(array)
