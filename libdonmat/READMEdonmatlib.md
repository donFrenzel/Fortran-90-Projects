# An Official Guide to libdonmat

libdonmat is an extensive matrix operations library built for Fortran 90, and though it is general use, I built it specifically for a few machine learning projects.  That being said, it is open and free to use.  I have included the module .f90 file as well, so that it can be better understood.  

## **How to use:**

Navigate to the libdonmat.a file in this repository folder, then click on it and find the download button that says 'Download Raw File' when you hover over it.  
Once downloaded, it is ready for use, in which you merely write 'use libdonmat' in your program file right after the program declaration (implicit none should still declared in the program body, just after the 'use' statement.  

When you would normally compile your code in the terminal use the format: *gfortran -I. -o filename.exe filename.f90 .\libdonmat.a*

Once your file is fully compiled, you may run it with ./filename.exe, as you would any other file.  You should be able to enter and use the functions present in the library as you please.  

## **Functions, Subroutines, & Syntax:** 

For this portion of the guide, I will be going over the various functions and subroutines present in the library in their order of appearance/complexity.  

### Functions: 

**readmatrix(filename.txt,#rows,#columns)**
  - This function takes a file name (only .txt for now) and two integer inputs, those being the number of rows and the number of columns present in your input matrix. 
  - It is recommended that, per the sample file, you take them as read(*,*)m,n and then allocate your dimension datatype to them for proper program input.  Though I must note that you can also basically just hardcode the matrix out in the program body as well.
  - It returns a matrix with the exact specifications given as input from the text file, which needs to have each value on a different line.

**Norm(vector)**
  - Takes a vector of any length as input and outputs its norm as a type real number.  

**sgn(value)**
  - Takes its input as a single value, whether it be integer, real, or double precision, and converts it to either a 1 or a -1 depending on its sign.  Simple sign function.  
### Subroutines:

**printmatrix(matrix)**
  - This function simply prints the matrix to the terminal.  It has no output and it's input is just the matrix.  
