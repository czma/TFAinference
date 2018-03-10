#!/usr/bin/python

#learn TFA matrix values from 
#	symbolic TFA matrix, 
#	CS matrix, 
#	and measured log expression values

import sys
from gurobipy import *
import time
import random

from TFAinferenceIO import *
from TFAinferenceMatrixMath import *


"""
Input:
  binary activity matrix A
  latest control strength matrix C
  expression matrix data
  list of model parameters
  fileLabel for use w/logging gurobi output
Output:
  False if learning failed
  learned A matrix if learning succeeded

makes a least squares optimization problem for gurobi to optimize
"""
def learnTFA(A, C, data, dataLabels, modelParams, fileLabel):
  [tfList, geneList, sampleList] = dataLabels

  numGenes = len(geneList)
  numSamples = len(sampleList)
  numTFs = len(tfList)

  #if validating, don't use scaling constraint
  #if zeroFlag, don't use zero constraint
  validating, zeroFlag = modelParams

  print "learning activity values"

  # initialize gurobi model
  model = Model()
  model.setParam('LogToConsole', False)
  model.setParam('LogFile', 'logFiles/'+fileLabel+".log")

  # Add tfa variables to the model
  varsMatrix = []   # holds the activity matrix, with pointers to coeff where relevant
  for i in range(numTFs):
    constraintCounter = 0   # counts the number of coeff in a row
    varsMatrix.append([])   # start a new row in the activity matrix
    constraint = LinExpr()    # initialize the constraint that each row's avg coeff value is 1
    for j in range(numSamples):
      if A[i][j]==0:
        varsMatrix[i].append(0)
      else:
	if zeroFlag:
          v = model.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name='A['+str(i)+','+str(j)+']')
	else:
          # crude lower bound to avoid learning of 0 for activity values
          v = model.addVar(lb=0.0001, vtype=GRB.CONTINUOUS, name='A['+str(i)+','+str(j)+']')
        varsMatrix[i].append(v)
        constraint += v
        constraintCounter += 1
    # add the scaling constraint if not validating
    if not validating:
      model.addConstr(constraint/constraintCounter, GRB.EQUAL, 1.0, "c"+str(i))
    model.update()

  # Populate objective
  obj = QuadExpr()
  for i in range(numGenes):
    for j in range(numSamples):
      geneExpr = LinExpr()
      geneExpr += C[i][numTFs]
      for k in range(numTFs):
        if (type(varsMatrix[k][j])==Var) and not (geneList[i] == sampleList[j] == tfList[k]):
          geneExpr += C[i][k]*varsMatrix[k][j]
      geneError = data[i][j] - geneExpr
      obj += geneError * geneError

  model.setObjective(obj)
  model.update()

  # Solve
  try:
  	model.optimize()
  except:
  	return False

  # Write model to a file
  # model.write('modelFiles/learnTFA'+fileLabel+'.lp')

  # check that optimization succeeded
  if model.status != GRB.Status.OPTIMAL:
    return False

  #convert back to matrix
  Atemp = []
  for i in range(numTFs):
    Atemp.append([])
    for j in range(numSamples):
      if A[i][j] == 0:
        Atemp[i].append(0)
      else:
        Atemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])
  Atemp.append([1]*numSamples)

  return Atemp


"""
Input:
  latest activity matrix A
  binary cs matrix C, or latest cs matrix if validating
  expression matrix data
  list of model parameters
Output:
  False if learning failed
  learned CS matrix if learning succeeded

makes a least squares optimization problem for gurobi to optimize
"""
def learnCS(A, C, data, dataLabels, modelParams, fileLabel):
  [tfList, geneList, sampleList] = dataLabels

  numGenes = len(geneList)
  numSamples = len(sampleList)
  numTFs = len(tfList)

  #csFlag is bool, whether or not cs signs are constrained to binary, ignored if validating
  #lassoWall is numeric value of upper bound for lasso constraint, ignored if validating
  #maFlag is bool, whether or not microarray data with no sign constraint on baseline value
  #validating is book, whether or not to learn TF-target gene influences at all
  #zeroFlag is bool, whether or not to allow zero coeff
  csFlag, lassoWall, maFlag, validating, zeroFlag = modelParams
  if validating: #C matrix has pseudoTF column, need to remove from count
    numTFs-=1

  print "learning CS with", numGenes, "genes,", numSamples, "samples,", numTFs, "TFs"
  print "lasso constraint:", lassoWall

  # Initialize the model
  model = Model()
  model.setParam('LogToConsole', False)
  model.setParam('LogFile', "logFiles/"+fileLabel+".log")
  
  # Add cs variables to the model
  varsMatrix = []   # holds the cs matrix, with pointers to coeff where relevant
  lassoConstraint = LinExpr()   # Intialize the LASSO constraint
  for i in range(numGenes):
    varsMatrix.append([])   # add a row to the cs matrix
    for j in range(numTFs+1):
      if j==numTFs: #learning baseline expression
        if maFlag:
          v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, \
            name='C['+str(i)+','+str(j)+']')
        else:
          v = model.addVar(lb=0.0001, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, \
            name='C['+str(i)+','+str(j)+']')
        varsMatrix[i].append(v)
      else: #learning an influence between a TF and gene
	if validating:
	  varsMatrix[i].append(C[i][j])
	else:
          if C[i][j]==0:  #no influence
            varsMatrix[i].append(0)
          elif C[i][j] > 0: #an influence to be learned, activating
            if csFlag:
              v = model.addVar(lb=0.0001, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, \
                name='C['+str(i)+','+str(j)+']')
            else:
                v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, \
                  name='C['+str(i)+','+str(j)+']')
            varsMatrix[i].append(v)
            v2 = model.addVar(name='|C['+str(i)+','+str(j)+']|')
            model.addGenConstrAbs(v2, v, "absconstr"+str(i)+"-"+str(j))
	    if not csFlag and not zeroFlag:
	      model.addConstr(v2 >= 0.0001, name='zeroconstr'+str(i)+'-'+str(j))
            lassoConstraint += v2
          else: #an influence to be learned, repressing
            if csFlag:
              v = model.addVar(lb=-GRB.INFINITY, ub=-0.0001, vtype=GRB.CONTINUOUS, \
                name='C['+str(i)+','+str(j)+']')
            else:
              v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, \
                name='C['+str(i)+','+str(j)+']')
            varsMatrix[i].append(v)
            v2 = model.addVar(name='|C['+str(i)+','+str(j)+']|')
            model.addGenConstrAbs(v2, v, "absconstr"+str(i)+"-"+str(j))
	    if not csFlag and not zeroFlag:
	      model.addConstr(v2 >= 0.0001, name='zerocontr'+str(i)+str(j))
            lassoConstraint += v2
  if lassoWall and not validating:
    model.addConstr(lassoConstraint <= lassoWall, "lasso")
  model.update()

  # Populate objective
  obj = QuadExpr()
  for i in range(numGenes):
    for j in range(numSamples):
      geneExpr = LinExpr()
      geneExpr += varsMatrix[i][numTFs]
      for k in range(numTFs):
        if type(varsMatrix[i][k])==Var and not (geneList[i] == sampleList[j] == tfList[k]):
          geneExpr += varsMatrix[i][k]*A[k][j]
      geneError = data[i][j] - geneExpr
      obj += geneError * geneError

  model.setObjective(obj)
  model.update()

  # Solve
  try:
    model.optimize()
  except:
    return False

  # Write model to a file
  model.write('modelFiles/learnCS'+fileLabel+'.lp')

  # check that optimization succeeded
  if model.status != GRB.Status.OPTIMAL:
    return False

  #convert back to matrix
  if validating:
    Ctemp = []
    for i in range(numGenes):
      Ctemp.append([])
      for j in range(numTFs+1):
        if j<numTFs:
          Ctemp[i].append(C[i][j])
        else:
          Ctemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])
  else:
    Ctemp = []
    for i in range(numGenes):
      Ctemp.append([])
      for j in range(numTFs+1):
        if j<numTFs and C[i][j] == 0:
          Ctemp[i].append(0)
        else:
          Ctemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])
    
  return Ctemp

"""
reads all input files
checks for discrepencies
returns false if discrepency found
otherwise return data from files
"""
def processInputFiles(inputFiles):
  startFile, csFile, tfaFile, dataFile, tfFile, geneFile, sampleFile = inputFiles
  # Put model data into matrices and lists
  # this function is in TFAinferenceIO.py
  A = readMatrixFromFile(tfaFile)
  C = readMatrixFromFile(csFile)
  Ctemp = readMatrixFromFile(startFile)
  data = readMatrixFromFile(dataFile)
  tfList = [x.strip() for x in open(tfFile, 'r').readlines()]
  geneList = [x.strip() for x in open(geneFile, 'r').readlines()]
  sampleList = [x.strip() for x in open(sampleFile, 'r').readlines()]
 
  #error checking
  if not (numGenes == len(C) == len(data) == len(Ctemp)):
    print "inconsistent number of genes\n gene list:", len(geneList), \
    "C matrix length:", len(C), "starting CS length:", len(Ctemp), "data length:", len(data)
    return False
  if not (numSamples == len(A[0]) == len(data[0])):
    print "inconsistent number of samples\n sample list:", len(sampleList), \
    "A matrix width:", len(A[0]), "data width:", len(data[0])
    return False
  if not (numTFs == len(A) == len(C[0]) == len(Ctemp[0])-1):
    print "inconsistent number of TFs\n TF list:", len(tfList), \
    "C matrix width:" len(C[0]), "starting CS width:", len(Ctemp[0]), "A matrix length:", len(A)
    return False

  return [A, C, Ctemp, data, tfList, geneList, sampleList]



"""
Input:
  list of input file names
  string for output file labeling
  integer value numIterations for how many iterations of optimization
  list of boolean and numerical flags
Output:
  the variance explained of the final model learned

Executes whole process of TFA inference learning
"""
def tfaInference(inputFiles, fileLabel, numIterations, modelParams):

  #zeroRelease is when to allow learning of 0 coeff values
  csFlag, lassoFlag, maFlag, validating, zeroRelease = modelParams

  print "validating?", validating
  
  # Put model data into matrices and lists
  [A, C, Ctemp, data, tfList, geneList, sampleList] = processInputFiles(inputFiles)

  numGenes = len(geneList)
  numSamples = len(sampleList)
  numTFs = len(tfList)
 
  cProgression = []
  aProgression = []
  if lassoFlag:
    # sample splits for 10 fold cross validation
    trainColsList, testColsList = cvSampleSplits(numSamples, 10)
    numCoeffCS = sum([sum([abs(x) for x in y]) for y in C])
  
  start = time.time()

  for itr in range(numIterations):
  
    print "\niteration ", itr, "\n"
  
    Atemp = learnTFA(A, Ctemp, data, [validating, itr > zeroRelease], \
      [tfList, geneList, sampleList], fileLabel)
    if Atemp == False:
      print "Could not learn the activity matrix"
      return

    if validating:
      Ctemp = learnCS(Atemp, Ctemp, data, [tfList, geneList, sampleList], \
        [csFlag, 0, maFlag, validating, itr > zeroRelease], fileLabel)
    else:  
      Ctemp = learnCS(Atemp, C, data, [tfList, geneList, sampleList], \
        [csFlag, 0, maFlag, validating, itr > zeroRelease], fileLabel)
    
    if Ctemp == False:
      print "Could not learn the control strength"
      return

    if lassoFlag and not validating:
      lassoLog = open("logFiles/lassoLog"+fileLabel+".tsv", 'a')
      # calculate lasso constraint upper bound as sum of abs coeff, except baseline values
      coeffSum = 0
      for i in range(len(Ctemp)):
        coeffSum += sum([abs(x) for x in Ctemp[i][:-1]])
      lassoLog.write(str(round(coeffSum, 3))+"\t")
      print "coeffSum:", coeffSum
      trainA, testA, trainData, testData = \
        cvSamplingMatrices(Atemp, data, trainColsList, testColsList)

      bestParam = 10
      bestError = float('inf')
      for param in range(1,11):
          errorList = []
          trainErrorList = []
          testDataCompilation = []
          realDataCompilation = []
          for cross in range(10):
            print "param", param, "cross", cross
            Ctest = learnCS(trainA[cross], C, trainData[cross], [tfList, geneList, sampleList], \
              [csFlag, param*0.1*coeffSum, maFlag, validating, itr > zeroRelease], fileLabel)
            if Ctest != False:
              validationData = matrixMultiply(Ctest, testA[cross])
              var, l2 = calcError(testData[cross], validationData, False)
              errorList.append(l2)
              var, l2 = calcError(trainData[cross], matrixMultiply(Ctest, trainA[cross]), False)
              trainErrorList.append(l2)
              for row in map(list, zip(*validationData)):
                testDataCompilation.append(row)
              for row in map(list, zip(*testData[cross])):
                realDataCompilation.append(row)
            else:
              lassoLog.write("could not learn with param "+str(param)+" on fold "+str(cross)+"\t")
          if sum(errorList)/len(errorList) < bestError:
            bestParam = param
            bestError = sum(errorList)/len(errorList)
          # log the results of testing this param
          lassoLog.write(str(param)+"\t")   # param
          # list of validation errors
          lassoLog.write("{"+",".join([str(round(x)) for x in errorList])+"}\t")
          lassoLog.write(str(round(sum(errorList)/len(errorList), 3))+"\t") # avg validation error
          #list of training errors
          lassoLog.write("{"+",".join([str(round(x)) for x in trainErrorList])+"}\t")
          lassoLog.write(str(round(sum(trainErrorList)/len(trainErrorList), 3))+"\t") # avg training error
          testVar, testL2 = calcError(map(list, zip(*realDataCompilation)), \
            map(list, zip(*testDataCompilation)), False)
          lassoLog.write(str(round(testVar,3))+"\t")  # var explained over all test results

      lassoLog.write(str(bestParam)+"\t") #best param
      lassoLog.write(str(round(bestError, 3))+"\t") #best avg validation error
      # learn CS with lasso constraint
      Ctemp = learnCS(Atemp, C, data, [tfList, geneList, sampleList], \
        [csFlag, bestParam*0.1*coeffSum, maFlag, validating, itr > zeroRelease], fileLabel)
      if Ctemp == False:
        print "Could not learn the control strength"
        return
      var, l2 = calcError(data, matrixMultiply(Ctemp, Atemp), False)
      lassoLog.write(str(round(l2,3))+"\n")   # error over all data
      lassoLog.close()

    aProgression.append(Atemp)
    cProgression.append(Ctemp) 

    currentVarExplained, currentError = calcError(data, matrixMultiply(Ctemp, Atemp), True)
    # log the results every 10 iterations
    if itr%10 == 0:
      saveResults(Ctemp, Atemp, currentVarExplained, \
        "logFiles/csLog"+fileLabel+".csv", "logFiles/tfaLog"+fileLabel+".csv", \
        "logFiles/varExplainedLog"+fileLabel+".csv")
  
  end = time.time()
  
  print "\n\n\n"
  print "done learning, now review:"
  
  for iteration in range(numIterations):
    C = cProgression[iteration]
    A = aProgression[iteration]
    print "iteration ", iteration
    dataTemp = matrixMultiply(C,A)
    var, l2 = calcError(data, dataTemp, True)
  
  print "total run time (secs): ", end-start
  
  saveResults(Ctemp, Atemp, var, \
    "results/learnedCS"+fileLabel+".csv", "results/learnedTFA"+fileLabel+".csv", \
    "results/learnedVarExplained"+fileLabel+".csv")

  return var

"""
Input:
  integer values for number of samples and number of folds
Output:
  two lists of column indices

generates randomized partitionings of samples across needed folds
"""
def cvSampleSplits(numSamples, folds):
  train = []
  test = []
  if numSamples < folds:
    print "not enough samples to do", folds, "fold cross validation"
    for i in range(folds):
      testCol = random.choice(range(numSamples))
      trainCols = [x for x in range(numSamples) if x != testCol]
      test.append([testCol])
      train.append(trainCols)
  else:
    unsampled = range(numSamples)
    for i in range(folds):
      testCols = random.sample(unsampled, len(unsampled)/(folds-i))
      #print "test", testCols
      trainCols = [x for x in range(numSamples) if x not in testCols]
      #print "train", trainCols
      unsampled = [x for x in unsampled if x not in testCols]
      #print "unsampled", unsampled
      train.append(trainCols)
      test.append(testCols)
  return [train, test]

"""
Input:
  two matrices of the same number of columns
  two lists of lists of column indices
Output:
  four lists of matrices

uses the column indices (these should be output from the cvSampleSplits function)
to make training and testing matrices
"""
def cvSamplingMatrices(Amatrix, dataMatrix, trainColsList, testColsList):
  trainA = []
  testA = []
  trainData = []
  testData = []
  for i in range(len(trainColsList)):
      testCols = testColsList[i]
      #print "test", testCols
      trainCols = trainColsList[i]
      #print "train", trainCols
      trainA.append(grabColumns(Amatrix, trainCols))
      testA.append(grabColumns(Amatrix, testCols))
      trainData.append(grabColumns(dataMatrix, trainCols))
      testData.append(grabColumns(dataMatrix, testCols))
  return [trainA, testA, trainData, testData]



















