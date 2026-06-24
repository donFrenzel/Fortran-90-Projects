# An Official Guide to libdonmat

libdonmat is a matrix operations library built in and for Fortran 90, and though it is general use, I built it specifically for a few machine learning projects.  That being said, it is open and free to use.  I have included the module .f90 file as well, so that it can be better understood.  

## **How to use:**

Navigate to the libdonmat.a file in this repository folder, then click on it and find the download button that says 'Download Raw File' when you hover over it.  
Once downloaded, it is ready for use, in which you merely write 'use libdonmat' in your program file right after the program declaration (implicit none should still declared in the program body, just after the 'use' statement.  (Refer to sample.f90 for a sample of how to call/run a program using this library)

When you would normally compile your code in the terminal use the format: *gfortran -I. -o sample.exe sample.f90 .\libdonmat.a*

Once your file is fully compiled, you may run it with *./filename.exe* as you would any other file.  You should be able to enter and use the functions present in the library as you please.  

## **Functions, Subroutines, & Syntax:** 

For this portion of the guide, I will be going over the various functions and subroutines present in the library in their order of appearance/complexity.  

### Functions: 

**readmatrix(filename.txt,#rows,#columns)**
  - This function takes a file name (only .txt for now) and two integer inputs, those being the number of rows and the number of columns present in your input matrix. 
  - It is recommended that, per the sample file, you take them as read(*,*)m,n and then allocate your dimension datatype to them for proper program input.  Though I must note that you can also basically just hardcode the matrix out in the program body as well.
  - It returns a matrix with the exact specifications given as input from the text file, which needs to have each value on a different line.

**norm(vector)**
  - Takes a vector of any length as input and outputs its norm as a type real number.  

**sgn(value)**
  - Takes its input as a single value, whether it be integer, real, or double precision, and converts it to either a 1 or a -1 depending on its sign.  Simple sign function.

**givens(a,b,c,s)**
  - Computes the givens rotations values for two real numbers (a,b) and outputs their cosine (c) and sine (s) accordingly.  

**roundSmalls(matrix,precision=1.00e-5)**
  - Takes a matrix as input and rounds down all values smaller than the optional precision definition to zero.  Outputs the rounded matrix.  
  - optional precision value input, with the default value being 1.00e-5.  You may set it as low or as high as you would like.

**eye(#rows,#columns,rightShift=0,downShift=0)**
  - My equivalent of the python eye() function.  It basically creates an identity matrix of the given shape, with the rightshift and downshift values selecting where the diagonal starts.  Output is an identity matrix with 1's along the diagonal (or in a diagonal pattern if specified)

**simplesort(array)**
  - Simplesort is just my own simple sorting algorithm I made for one of the other functions.  It isn't very efficient, but it basically extracts the maximum value from one array and places it into the sorted array.  It sorts them in descending order and returns the sorted array.  

**rowswap(matrix,row1,row2)**
  - Rowswap takes an input matrix and two row indices, then swaps the rows at those indices in the input array, and then outputs a copy of the input array with those rows swapped.  

**Gaussian(matrix)**
  - Gaussian is simply the function name for Gaussian Elimination.  It converts the matrix into upper triangular form, though does not reduce to leading 1's for the output (this can be done easily enough later on).  The output for this function is the upper triangular matrix without leading 1's.

**RREF(matrix)**
  - RREF computes the Reduced Row Echelon Form for the input matrix, basically just the full reduction into leading 1's form for a square matrix.

**det(matrix)**
  - Det computes the determinant of a square matrix and outputs the Determinant found as a type real value.  

**gramdeterminant(matrix)**
  - This function computes the Gram Determinant of a given matrix, which is just taking the determinant of A<sup>T</sup>A.  It returns it as a type real value.  

**inverse(matrix)**
  - Inverse returns the inverse matrix of the input matrix, so long as that matrix is invertible.  It uses LU Decomposition for this.  It returns the inverse matrix of the input.  

**pseudoinverse(matrix)**
  - This function takes the Moore-Penrose Pseudo Inverse of a matrix which would not normally be invertible.  The output is the aforementioned matrix, whose size   will be the transpose of the initial matrix.  

**linearindependence(inMatrix)**
  - This function just outputs a 0 if linearly dependent and a 1 if linearly independent and takes a matrix as its only argument.

**outerproduct(vector1,vector2)**
  - This function outputs the matrix outer product of two vectors.  Input is two vectors, the output is their outer product in matrix form.

**kroneckerproduct(matrix1,matrix2)**
  - This function computes the kronecker product of two matrices.  It takes two matrices as input and has their product as output.  

**hessenberg(matrix)**
  - This function takes a matrix as input and converts it into Upper Hessenberg Form, which happens to be the output.  

**eigenvalshouseholder(matrix, iterations=50)**
  - This function computes the eigenvalues of any matrix using the QR Algorithm with Householder Reflections.  It outputs the eigenvalues in sorted order, from greatest to least, in a single array.  It has a default of 50 iterations, but may be set to whatever you may please. 

**eigenvectors(matrix,eigenvalues)**
  - This function takes a matrix and an array of eigenvalues as input arguements and outputs the matrix of eigenvectors.  

**solvesystem(matrix)**
  - This function takes a matrix as input and basically reduces it to gaussian form to then solve for the x1,x2,...,xn present in the matrix.  It does so in the far-right column, and outputs the values there as an array of type real 'solutions'.  Useful for solving systems of equations.  

### Subroutines:
*Note, all subroutines require you to write 'call (name of subroutine)'*
*Also, subroutine outputs must be entered as inputs, as subroutines can only use existing input values as outputs, hence why all of the below (except print) have some of their inputs listed as outputs*

**printmatrix(matrix)**
  - This subroutine simply prints the matrix to the terminal.  It has no output and it's input is just the matrix.

**LU(matrix,L,U)**
  - This subroutine computes the LU Decomposition of the input matrix.  It takes three matrices are arguments, though two are rewritten by the subroutine (L & U), and returns L & U in place back to the program in their current forms.  

**QR(matrix,Q,R)**
  - This subroutine takes a matrix as input and computes the QR Decomposition using Householder Reflections.  It takes Q & R as inputs and uses them as outputs.  Works better and is more generally applicable than gramschmidt QR.  Utilized in eigenvalshouseholder().  

**gramschmidtQR(matrix,Q,R)**
  - This subroutine computes the Gram-Schmidt Method for QR Decomposition.  It takes three matrices as arguments, those being the input matrix, Q, and R, which are reallocated to match the needs of the input matrix and are the outputs of the subroutine.  

**SVD(matrix,U,S,VT)**
  - This subroutine computes the Singular Value Decomposition of the input matrix and has input/output's U, S, and VT (U, Sigma, & V-transpose, respectively).  
