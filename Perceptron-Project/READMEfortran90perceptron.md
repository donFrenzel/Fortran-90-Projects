# Fortran-90-Perceptron Project

This is going to be my main repository for my Fortran 90 Perceptron Project.  The goal here is to create a simple perceptron which can take label and feature inputs and then predict unknowns based on the learned weights.  It is, in essence, a recreation of a Perceptron I wrote in Python in a Machine Learning Course I took.  

The components are: 
- A Fortran 90 implementation of the classic Perceptron Learning Algorithm 
- A CSV File Reader which takes a .csv file and converts it to a .txt file so that I can then use libdonmat's matrix feature for Perceptron input.  

To use the Perceptron, you must update the file names in the .f90 file to be your exact .csv's.  It takes .csv input and converts to .txt format for entry into libdonmat, where then the resulting matrices (which are the converted csv data) can be be used as input for the perceptron.  Current readouts do not include the confusion matrix, recall, accuracy, etc.
