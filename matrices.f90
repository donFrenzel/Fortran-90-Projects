program matrices
!! When creating a matrix, the count for rows and for columns begins at ONE and not zero like other arrays
implicit none
integer,dimension(2,3)::mymatrix1
integer,dimension(2)::q
integer::b
q=maxloc(myMatrix1)  !!Prints in row,column form.  
myMatrix1(2,3)=5 !!Means you can set the value within a specific matrix to a specific value
b=myMatrix1(2,3)
write(*,*)q,b
stop
end program matrices