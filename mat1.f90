program matrixMult
implicit none
integer,dimension(2,2)::A,B,O !!Request two 2x2 integer matrices, A & B.  A & B Input, O output.
write(*,*)'Please enter the information for the first matrix [FOUR VALUES]'
read(*,*)A(1,1),A(1,2),A(2,1),A(2,2)
write(*,*)'Please enter the information for the second matrix. [FOUR VALUES]'
read(*,*)B(1,1),B(1,2),B(2,1),B(2,2)
O(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)!!Now we need to multiply them together to get the total
O(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)
O(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)
O(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)
write(*,*)O(1,1),O(1,2)
write(*,*)O(2,1),O(2,2)
stop
end program matrixMult