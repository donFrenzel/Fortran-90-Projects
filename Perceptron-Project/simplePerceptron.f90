module csvtotxt
implicit none
contains
!fill here for module.
subroutine csv_to_txt(filename,filenameC,rows,cols) 
implicit none 
character(len=*), intent(in)::filename
integer::rows,cols
integer::eof,suffixPos,i,commaPosition,startPos
integer::numCols,numRows
character(len=4)::oldSuf,newSuf
character(len=255)::filenameC,field
character(len=2048)::firstLine,curLine !for first line of the file and current file line (for dividing up and entering as data into txt.  )
!define old and new suffixes for transfusion
suffixPos=0
oldSuf=".csv"
newSuf=".txt"
filenameC=filename
suffixPos=index(filenameC,oldSuf)!find position of the old .csv suffix AND replace with new suffix, new filename to open up. 
if (suffixPos.gt.0) then 
    filenameC=filenameC(:suffixPos-1)//newSuf//filenameC(suffixPos+len(oldSuf):)
end if
eof=0
numCols=1 !will always have at least 1 data column for the csv not to be empty.  If you feed it an empty CSV its your fault.  
numRows=0
open(15,file=filename,status='old')
open(20,file=filenameC,status='unknown',action="write")
read(15, '(A)', IOSTAT=eof) firstLine !read to the first line.  
!NOTE: len_trim() intrinsic is the trimmed (no nullspace) length of the character array being entered.  Loop thru and get comma delimit count
do i=1,len_trim(firstLine)
    if(firstLine(i:i)==',') then 
        numCols=numCols+1
    end if
end do
do 
    read(15,'(A)',iostat=eof)curLine !now we need to take the number of columns and separate the data into those x columns PER UNIT
    if(eof<0) then 
        exit
    end if
    startPos=1
    commaPosition=1
    !now that it is curLine holding the line data, we take the numCols and divide the data up into them by finding the index
    do
        !get the index of the comma position and check if there is no next comma, in which case just grab the data until the end. 
        commaPosition=index(curLine(startPos:),',')
        if(commaPosition==0) then 
            field=curLine(startPos:)
            write(20,*)trim(adjustl(field)) !write to the txt the last bit (since that's what we're at in this if statement)
            exit
        else
            !now basically we're gonna need to get the comma position and read everything from that until the next comma.  
            commaPosition= startPos+commaPosition-1 !basically here because you're slicing a small portion with INDEX() above, must be added to 
            !total position.  Caused many problems.  Not any more.  
            field=curLine(startPos:commaPosition-1) !grab between the position of the start and the comma (-1 necessary because it does index that value.  )
            write(20,*)trim(adjustl(field)) 
            startPos=commaPosition+1
        end if
    end do
    numRows=numRows+1
end do
rows=numRows
cols=numCols
!close out the connections.
close(15)
close(20)
end subroutine csv_to_txt
end module csvtotxt
!Perceptron module here.  
module perceptronModule
implicit none
contains
!write primary logic (for matrix input)
!input is a matrix nxm which is tailored so that one may be processed as the label (so long as it is either a 1 or 0)
!and the others are to be taken as the 'coordinates' or the classifying variables.  
!predict takes a single row along with input weights and produces it already.  
function predict(row,weights) RESULT(prediction)
implicit none
!inputs here
real, dimension(:),intent(in)::weights, row
integer::i,n
real::prediction,activation
n=size(row)
activation = weights(1)
!loop should skip the first value (for this instance, since it is the classifying binary one)
do i=2,n,1
    activation = activation + weights(i)*row(i)
end do
!now assign the activation
if(activation.ge.0.0) then 
    prediction=1.0
else
    prediction=0.0
end if
end function predict

!train the weights vector for the perceptron
subroutine train(inMatrix, weights, learningRate, numIters)
implicit none
!input/optionals
integer,intent(in),optional::numIters !create checks
real,intent(in),optional::learningRate !create checks
real, dimension(:,:),intent(in)::inMatrix
real,dimension(:),intent(inout)::weights
!used by function main body
integer::n,m, actNumIters, i,j,k
real::actLearningRate,sumErr,prediction, err
real,dimension(:),allocatable::curRow !allocate 

!Get input matrix size
n = size(inMatrix,dim=1) !#rows
m = size(inMatrix,dim=2) !#cols

!assign default values
if(present(learningRate)) then 
    actLearningRate=learningRate
else
    actLearningRate=0.001
end if
if(present(numIters)) then 
    actNumIters=numIters
else
    actNumIters=100
end if
!allocatables
allocate(curRow(m)) !take all row values

weights=0.0
!main body of the function
do i=1, actNumIters, 1
    sumErr = 0.0
    !for each row in the dataset
    do j=1,n,1
        !grab the row data
        curRow=inMatrix(j,:) !shouldn't overindex EVER.  IF SO BIG PROBLEM.
        prediction=predict(curRow,weights)
        err=curRow(1)-prediction !since val 1 will always be the label
        sumErr=sumErr+err**2
        weights(1)=weights(1)+actLearningRate*err
        do k=2,m,1
            weights(k) = weights(k)+actLearningRate*err*curRow(k) !nevermind this value.  
        end do
    end do
    if(modulo(i,5)==0) then 
        write(*,*)"Epoch: ",i,"learnRate: ",actLearningRate,"Error: ",sumErr
    end if
end do


!deallocate and wrap up
deallocate(curRow)
end subroutine train

!actual perceptron (takes a training dataset AND a testing dataset)
subroutine perceptron(trainSet,testSet,predictions,learningRate,numIters)
implicit none
!def inputs
integer,intent(in),optional::numIters
real,intent(in),optional::learningRate
real, dimension(:,:),intent(in)::trainSet,testSet
!def outputs
real,dimension(:),allocatable,intent(out)::predictions
!define program variables
integer::n1,m1,n2,m2,actNumIters,i !1 is for train, 2 is for test
real::actLearningRate,prediction
real, dimension(:),allocatable::weights,testRow
!get the sizes of the training and testing set
n1=size(trainSet,dim=1) !#rows
m1=size(trainSet,dim=2) !#cols
n2=size(testSet,dim=1)
m2=size(testSet,dim=2)
!Check to make sure column values match; if not could really mess up the program.  Idiot-proofing.
if(m1/=m2)then 
    write(*,*)"Train/Test feature count does not match.  Exiting program."
    stop
end if
!fill in default values
if(present(numIters)) then
    actNumIters=numIters
else
    actNumIters=100
end if
if(present(learningRate)) then
    actLearningRate=learningRate
else
    actLearningRate=0.001
end if
!allocatables
allocate(predictions(n2)) !allocate to match the number of rows present in the test set. 
allocate(weights(m1)) !allocate weights for m1-1 (both m's same)
allocate(testRow(m2)) !allocate to length of test set feature and label combination vector. 
!main body
call train(trainSet,weights,learningRate,numIters) !problem is with weights
do i=1,n2,1
    testRow=testSet(i,:) !get whole testing row
    prediction=predict(testRow,weights)
    predictions(i)=prediction
end do
deallocate(weights)
deallocate(testRow)
end subroutine perceptron
end module perceptronModule

!for calling all modules + the general predictive model.  
!to compile: gfortran -I. -o simplePerceptron.exe simplePerceptron.f90 .\libdonmat.a
!to run: ./simplePerceptron.exe
program simplePerceptron
!library/module declarations
use csvtotxt
use perceptronModule
use libdonmat !for matrix conversion & input
implicit none
!load in the data now
integer::numR_train,numC_train,numR_test,numC_test
character(len=255)::trainTxt,testTxt
real,dimension(:),allocatable::testPreds,testActual
real,dimension(:,:),allocatable::trainMatrix,testMatrix
!read the training data
call csv_to_txt('sampleperceptrontrain.csv',trainTxt,numR_train,numC_train)
!read the testing data
call csv_to_txt('sampleperceptrontest.csv',testTxt,numR_test,numC_test)
!allocate
allocate(trainMatrix(numR_train,numC_train))
allocate(testMatrix(numR_test,numC_test))
allocate(testPreds(numR_test))
allocate(testActual(numR_test))
!convert them both to matrices
trainMatrix=readmatrix(trainTxt,numR_train,numC_train)
testMatrix=readmatrix(testTxt,numR_test,numC_test)
!then call the perceptron test 
call perceptron(trainMatrix,testMatrix,testPreds,0.001,50)
!and now compare the testPreds with the testActual
testActual=testMatrix(:,1) !get only the labels
!length of testPreds should be exactly equal so now loop through and keep track of the amount of errors in comparison to the whole.  
write(*,*)testActual
write(*,*)testPreds
!comparison of test and predictions to come soon.  

end program simplePerceptron
