program findresults
implicit none
integer,parameter::numStudents=5 !Wrote the number of students as a parameter.
character*5,dimension(numStudents)::students
character*15,dimension(5)::passfail !!Allows for max 15 chars
integer,dimension(5)::grades
integer::numPasses,i
real::avg
passfail=''!!Mass assignment of blank variable to all array variables.
students(1)='Fred'!!Assignment of all names and grades.
students(2)='Susie'
students(3)='Tom'
students(4)='Anita'
students(5)='Peter'
grades(1)=64
grades(2)=57
grades(3)=49
grades(4)=71
grades(5)=37
numPasses=Count(grades>50)
avg=real(sum(grades))/real(size(grades))
write(*,11)students
11 format('Students:'1x5A7)
close(11)
write(*,12)grades
12 format('Grades:'5i7)
close(12)
where(grades>=50) passfail='P'
write(*,13)passfail
13 format('Pass?'7x5A7)
close(13)
write(*,9)numPasses
9 format('Number of Passes: 'i1)!! Integer of length 1
close(9)
write(*,10)avg
10 format('Average Grade: '1f6.1)
i=maxloc(grades, dim=1)
close(10)
write(*,14)students(i)
14 format('Best grade: 'A5)
close(14)
stop
end program findresults