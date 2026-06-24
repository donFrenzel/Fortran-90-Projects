# An Official Guide to libdonmat

libdonmat is an extensive matrix operations library built for Fortran 90, and though it is general use, I built it specifically for a few machine learning projects.  That being said, it is open and free to use.  

**How to use:**

Navigate to the libdonmat.a file in this repository folder, then click on it and find the download button that says 'Download Raw File' when you hover over it.  
Once downloaded, it is ready for use, in which you merely write 'use libdonmat' in your program file right after the program declaration (implicit none should still declared in the program body, just after the 'use' statement.  

When you would normally compile your code in the terminal use the format: gfortran -I. -o filename.exe filename.f90 .\libdonmat.a
Once your file is fully compiled, you may run it with ./filename.exe, as you would any other file.  You should be able to enter and use the functions present in the library as you please.  

**Functions, Subroutines, & Syntax:** 

For this portion of the guide, I will be 
