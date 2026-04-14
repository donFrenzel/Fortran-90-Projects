program orderOfOps
implicit none
real::a,b,c,d,e,f,g,x,y,z
a=10.5
b=5.5
c=8.9 !!not sure why its doing this
d=9.2
e=11.8
f=14.6
g=20.26
x=(a*b)+(c*d)+(e/(f**g))
y=(a*(b+c)*d)+((e/f)**g)
z=(a*(b+c))*((d+e)/(f**g))
write(*,*)'x is ',x,'y is ',y,'z is ',z
stop
end program orderOfOps
