program convertBinary
implicit none
real::i,j,k
integer::r,c,q,l
integer,allocatable::binary(:)
write(*,*)'Please input the number you wish to be converted into binary.'
read(*,*)i
q=0
if(i.lt.0) then
    q=1 !Signifies if it's a negative. 
    i=-i
    i=int(i)
end if
r=0
k=i
do while(i>0.9) !In essence, similar to both python and c
    j=modulo(i,2.0) !Modulo takes the current integer and then a divisor, BOTH MUST BE REAL()
    i=i/2
    r=r+1
end do
l=4-modulo(r,4)! l is the effective displacement, line below puts it into effect.  
r=r+l
!!Allocate the array
allocate(binary(r))
c=r
do while(k>0.9)
    j=modulo(k,2.0) !Modulo takes the current integer and then a divisor, BOTH MUST BE REAL()
    k=k/2
    binary(c)=int(j)
    c=c-1
end do
!Handles filling in the bit-extended pieces.
do while(l>-1)
    binary(l)=0
    l=l-1
end do
!Check if it's positive or negative using sign variable
if(q.ne.1) then 
    write(*,*)binary
else
    call twosComplement(binary) !!Call statement very important for a subroutine
end if
deallocate(binary)
stop
!
contains 
subroutine twosComplement(arr)
implicit none
integer,dimension(:), intent(inout)::arr
integer::i
do i=0,size(arr),1 !This loop flips all of the bits
    if(arr(i).ne.1) then
        arr(i)=arr(i)+1
    else
        arr(i)=arr(i)-1
    end if
end do
do i=size(arr),0,-1 !This loop handles adding 1 to the bit-flipped binary number.
   if(arr(i).eq.0) then
        arr(i)=arr(i)+1
        exit
   else
        arr(i)=0
   end if
end do
write(*,*)arr
return
end subroutine twosComplement
end program convertBinary