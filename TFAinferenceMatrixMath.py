import sys
from gurobipy import *
import time

from TFAinferenceIO import *

"""
matrix multiplication
Someday I will import numpy and this won't be used anymore ...
"""
def matrixMultiply(C, A):
  if len(C[0]) == len(A):
    result = []
    for i in range(len(C)):
      result.append([])
      for j in range(len(A[0])):
        val = 0
        for k in range(len(A)):
          val += C[i][k]*A[k][j]
        result[i].append(val)
    return result
  else:
    print "matrices are the wrong dimensions for multiplication"
    print len(C[0]), len(A)
    return []

"""
Input:
  activity matrix, and a binary matrix of factor labels (factor x sample matrix)
  a list of lists of factor values
Output:
  activity matrix with factor values applied

This is for datasets such that certain samples can be grouped as having some factor
such that the factor is expected to have a consistent effect on TF activity
Currently not used anywhere
"""
def factorMultiplyA(A, F, fVals):
  newA = [x[:] for x in A]
  for i in range(len(fVals[0])): #TF
    for j in range(len(A[i])): #sample
      for fIndex in range(len(F)): #factor
        if F[fIndex][j]==1:
          newA[i][j] *= fVals[fIndex][i]
  return newA
          

"""
Input:
  two matrices, where the first matrices #col = #row of second
Output:
  two rescaled matrices, such that the non-zero values of each row of the second is 1

one possible way to normalize the learned CS and activity matrices
currently not used anywhere
"""
def rescale(C, A):
  xMatrix = []
  xMatrixInverse = []
  for i in range(len(A)):
    xMatrix.append([0]*len(A))
    xMatrixInverse.append([0]*len(A))
    xMatrix[i][i]=1.0/(sum(abs(x) for x in A[i])/len([x for x in A[i] if x != 0.0]))
    xMatrixInverse[i][i]=(sum(abs(x) for x in A[i])/len([x for x in A[i] if x != 0.0]))

  Crescaled = matrixMultiply(C, xMatrixInverse)
  Arescaled = matrixMultiply(xMatrix, A)
  return [Crescaled, Arescaled]

"""
Input:
  two data matrices
  the first is assumed to be true data
  the second is assumed to be predicted/learned/fitted data
Output:
  list of var explained and SSE

calculates total variance of data, SSE of dataLearned-data
and variance explained of dataLearned
"""
def calcError(data, dataLearned, printFlag):
  numerator = 0
  denominator = 0
  error = 0
  numValues = 0
  for i in range(len(dataLearned)):
    geneMean = sum(data[i])/len(data[i])
    for j in range(len(dataLearned[i])):
      numerator += (data[i][j]-dataLearned[i][j]) * (data[i][j]-dataLearned[i][j])
      denominator += (data[i][j] - geneMean) * (data[i][j] - geneMean)
      numValues += 1
  if printFlag:
    print "variance: ", denominator/numValues
    print "variance explained: ", 1-(numerator/denominator)
    print "Sum of Squared Errors: ", numerator
  return [1-(numerator/denominator), numerator]

"""
Input:
  a list of lists
  a list of column indices
Output:
  a list of lists

Returns a matrix of only the columns listed from the input matrix
If a column index in the list is beyond the input matrix, it's ignored
"""
def grabColumns(matrix, cols):
  newMatrix = []
  for row in matrix:
    newMatrixRow = []
    for j in range(len(row)):
      if j in cols:
  newMatrixRow.append(row[j])
    newMatrix.append(newMatrixRow)
  return newMatrix






  
