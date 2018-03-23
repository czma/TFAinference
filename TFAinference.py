#!/usr/bin/python

# learn TFA matrix values from
#	symbolic TFA matrix, 
#	CS matrix, 
#	and measured log expression values

from gurobipy import *
import numpy as np
import time
from TFAinferenceMatrixMath import *
from TFAinferenceIO import *

"""
Input:
  binary activity matrix A
  latest control strength matrix C
  expression matrix data
  list of model parameters
Output:
  False if learning failed
  learned A matrix if learning succeeded

makes a least squares optimization problem for gurobi to optimize
"""


def learnTFA(A, C, data, modelParams):
    numGenes = len(data)
    numSamples = len(data[0])
    numTFs = len(C[0]) - 1

    # currently none of these params relevant to learning TFA
    # csFlag, lassoFlag, maFlag = modelParams

    print("learning activity values")

    # initialize gurobi model
    model = Model()
    model.setParam('OutputFlag', False)

    # Add tfa variables to the model
    varsMatrix = []  # holds the activity matrix, with pointers to coeff where relevant
    for i in range(numTFs):
        constraintCounter = 0  # counts the number of coeff in a row
        varsMatrix.append([])  # start a new row in the activity matrix
        constraint = LinExpr()  # initialize the constraint that each row's avg coeff value is 1
        for j in range(numSamples):
            if A[i][j] == 0:
                varsMatrix[i].append(0)
            else:
                # currently does not allow learning of 0 for activity values
                v = model.addVar(lb=0.0001, vtype=GRB.CONTINUOUS, name='A[' + str(i) + ',' + str(j) + ']')
                varsMatrix[i].append(v)
                constraint += v
                constraintCounter += 1
        # add the scaling constraint
        model.addConstr(constraint / constraintCounter, GRB.EQUAL, 1.0, "c" + str(i))
        model.update()

    # Populate objective
    obj = QuadExpr()
    for i in range(numGenes):
        for j in range(numSamples):
            geneExpr = LinExpr()
            geneExpr += C[i][numTFs]
            for k in range(numTFs):
                if type(varsMatrix[k][j]) == Var:
                    geneExpr += C[i][k] * varsMatrix[k][j]
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
    # model.write('learnTFA.lp')

    # check that optimization succeeded
    if model.status != GRB.Status.OPTIMAL:
        return False

    # convert back to matrix
    Atemp = []
    for i in range(numTFs):
        Atemp.append([])
        for j in range(numSamples):
            if A[i][j] == 0:
                Atemp[i].append(0)
            else:
                Atemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])
    Atemp.append([1] * numSamples)

    return Atemp


"""
Input:
  latest activity matrix A
  binary cs matrix C
  expression matrix data
  list of model parameters
Output:
  False if learning failed
  learned CS matrix if learning succeeded

makes a least squares optimization problem for gurobi to optimize
"""


def learnCS(A, C, data, modelParams):
    numGenes = len(data)
    numSamples = len(data[0])
    numTFs = len(C[0])

    csFlag, lassoWall, maFlag = modelParams

    print("learning CS with", numGenes, "genes,", numSamples, "samples,", numTFs, "TFs")

    # Initialize the model
    model = Model()
    model.setParam('OutputFlag', False)

    # Add cs variables to the model
    varsMatrix = []  # holds the cs matrix, with pointers to coeff where relevant
    lassoConstraint = LinExpr()  # Intialize the LASSO constraint
    for i in range(numGenes):
        varsMatrix.append([])  # add a row to the cs matrix
        for j in range(numTFs + 1):
            if j == numTFs:  # learning baseline expression
                if maFlag:
                    v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                     name='C[' + str(i) + ',' + str(j) + ']')
                else:
                    v = model.addVar(lb=0.0001, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                     name='C[' + str(i) + ',' + str(j) + ']')
                varsMatrix[i].append(v)
            else:  # learning an influence between a TF and gene
                if C[i][j] == 0:  # no influence
                    varsMatrix[i].append(0)
                elif C[i][j] > 0:  # an influence to be learned, activating
                    if csFlag:
                        v = model.addVar(lb=0.0001, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                         name='C[' + str(i) + ',' + str(j) + ']')
                    else:
                        v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                         name='C[' + str(i) + ',' + str(j) + ']')
                    varsMatrix[i].append(v)
                    v2 = model.addVar()
                    model.addGenConstrAbs(v2, v, "absconstr" + str(i) + "-" + str(j))
                    lassoConstraint += v2
                else:  # an influence to be learned, repressing
                    if csFlag:
                        v = model.addVar(lb=-GRB.INFINITY, ub=-0.0001, vtype=GRB.CONTINUOUS,
                                         name='C[' + str(i) + ',' + str(j) + ']')
                    else:
                        v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                         name='C[' + str(i) + ',' + str(j) + ']')
                    varsMatrix[i].append(v)
                    v2 = model.addVar()
                    model.addGenConstrAbs(v2, v, "absconstr" + str(i) + "-" + str(j))
                    lassoConstraint += v2
    if lassoWall:
        model.addConstr(lassoConstraint <= lassoWall, "lasso")
    model.update()

    # Populate objective
    obj = QuadExpr()
    for i in range(numGenes):
        for j in range(numSamples):
            geneExpr = LinExpr()
            geneExpr += varsMatrix[i][numTFs]
            for k in range(numTFs):
                if type(varsMatrix[i][k]) == Var:
                    geneExpr += varsMatrix[i][k] * A[k][j]
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
    # model.write('learnCS.lp')

    # check that optimization succeeded
    if model.status != GRB.Status.OPTIMAL:
        return False

    # convert back to matrix
    Ctemp = []
    for i in range(numGenes):
        Ctemp.append([])
        for j in range(numTFs + 1):
            if j < numTFs and C[i][j] == 0:
                Ctemp[i].append(0)
            else:
                Ctemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])

    return Ctemp


"""
Input:
  latest activity matrix A
  validation cs matrix C
  expression matrix data
  list of model parameters
Output:
  False if learning failed
  learned CS matrix if learning succeeded

makes a least squares optimization problem for gurobi to optimize
"""


def learnCSValidation(A, C, data, modelParams):
    numGenes = len(data)
    numSamples = len(data[0])
    numTFs = len(C[0]) - 1

    csFlag, lassoWall, maFlag = modelParams

    print("learning CS with", numGenes, "genes,", numSamples, "samples,", numTFs, "TFs")

    # Initialize the model
    model = Model()
    model.setParam('OutputFlag', False)

    # Add cs variables to the model
    varsMatrix = []
    for i in range(numGenes):
        varsMatrix.append([])
        for j in range(numTFs + 1):
            if j == numTFs:  # learning baseline expression
                if maFlag:
                    v = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                     name='C[' + str(i) + ',' + str(j) + ']')
                else:
                    v = model.addVar(lb=0.0001, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS,
                                     name='C[' + str(i) + ',' + str(j) + ']')
                varsMatrix[i].append(v)
            else:  # learning an influence between a TF and gene
                varsMatrix[i].append(C[i][j])
    model.update()

    # Populate objective
    obj = QuadExpr()
    for i in range(numGenes):
        for j in range(numSamples):
            geneExpr = LinExpr()
            geneExpr += varsMatrix[i][numTFs]
            for k in range(numTFs):
                if type(varsMatrix[i][k]) == Var:
                    geneExpr += varsMatrix[i][k] * A[k][j]
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
    # model.write('learnCS.lp')

    # check that optimization succeeded
    if model.status != GRB.Status.OPTIMAL:
        return False

    # convert back to matrix
    Ctemp = []
    for i in range(numGenes):
        Ctemp.append([])
        for j in range(numTFs + 1):
            if j < numTFs:
                Ctemp[i].append(C[i][j])
            else:
                Ctemp[i].append(model.getAttr('x', [varsMatrix[i][j]])[0])

    return Ctemp


"""
Input:
  list of input file names
  string for output file labeling
  integer value numIterations for how many iterations of optimization
  list of boolean flags for CS constraints, LASSO contraints, and microarray/RNAseq data
Output:
  the variance explained of the final model learned

Executes whole process of TFA inference learning
"""


def tfaInference(inputFiles, fileLabel, numIterations, modelParams, identicalIterations = 3, foldChangeLimit = 0.1, errorChangeLimit=1.0):
    startFile, csFile, tfaFile, dataFile = inputFiles
    csFlag, lassoFlag, maFlag = modelParams

    # Put model data into matrices and lists
    # this function is in TFAinferenceIO.py
    A = readMatrixFromFile(tfaFile)
    Atemp = A
    C = readMatrixFromFile(csFile)
    Ctemp = readMatrixFromFile(startFile)
    data = readMatrixFromFile(dataFile)
    var = 0

    numGenes = len(data)
    numSamples = len(data[0])
    numTFs = len(A)

    cProgression = []
    aProgression = []
    if lassoFlag:
        # sample splits for 10 fold cross validation
        trainColsList, testColsList = cvSampleSplits(numSamples, 10)

    start = time.time()

    identicalCounter = 0
    prevError = 0


    for iter in range(numIterations):

        print("\niteration ", iter, "\n")

        Atemp = learnTFA(A, Ctemp, data, [])
        if Atemp == False:
            print("Could not learn the activity matrix")
            return

        Ctemp = learnCS(Atemp, C, data, [csFlag, 0, maFlag])

        if Ctemp == False:
            print("Could not learn the control strength")
            return

        if lassoFlag:
            lassoLog = open("logFiles/lassoLog" + fileLabel + ".tsv", 'a')
            # calculate lasso constraint upper bound as sum of abs coeff, except baseline values
            coeffSum = 0
            for i in range(len(Ctemp)):
                coeffSum += sum([abs(x) for x in Ctemp[i]])
                coeffSum -= abs(Ctemp[i][-1])
            lassoLog.write(str(round(coeffSum, 3)) + "\t")
            trainA, testA, trainData, testData = cvSamplingMatrices(Atemp, data, trainColsList, testColsList)

            bestParam = 10
            bestError = float('inf')
            for param in range(1, 11):
                errorList = []
                trainErrorList = []
                testDataCompilation = []
                realDataCompilation = []
                for cross in range(10):
                    print("param", param, "cross", cross)
                    Ctest = learnCS(trainA[cross], C, trainData[cross], [csFlag, param * 0.1 * coeffSum, maFlag])
                    if Ctest != False:
                        validationData = np.dot(Ctest, testA[cross])
                        var, l2 = calcError(testData[cross], validationData, False)
                        errorList.append(l2)
                        var, l2 = calcError(trainData[cross], np.dot(Ctest, trainA[cross]), False)
                        trainErrorList.append(l2)
                        for row in map(list, zip(*validationData)):
                            testDataCompilation.append(row)
                        for row in map(list, zip(*testData[cross])):
                            realDataCompilation.append(row)
                    else:
                        lassoLog.write("could not learn with param " + str(param) + " on fold " + str(cross) + "\t")
                #print(errorList)
                if sum(errorList) / len(errorList) < bestError:
                    bestParam = param
                    bestError = sum(errorList) / len(errorList)
                # log the results of testing this param
                lassoLog.write(str(param) + "\t")  # param
                # list of validation errors
                lassoLog.write("{")
                for x in errorList:
                    lassoLog.write(str(round(x, 3)) + ", ")
                lassoLog.write("}\t")
                lassoLog.write(str(round(sum(errorList) / len(errorList), 3)) + "\t")  # avg validation error
                # list of training errors
                lassoLog.write("{")
                for x in trainErrorList:
                    lassoLog.write(str(round(x, 3)) + ", ")
                lassoLog.write("}\t")
                lassoLog.write(str(round(sum(trainErrorList) / len(trainErrorList), 3)) + "\t")  # avg training error
                testVar, testL2 = calcError(list(map(list, zip(*realDataCompilation))), list(map(list, zip(*testDataCompilation))),
                                            False)
                lassoLog.write(str(round(testVar, 3)) + "\t")  # var explained over all test results

            lassoLog.write(str(bestParam) + "\t")  # best param
            lassoLog.write(str(round(bestError, 3)) + "\t")  # best avg validation error
            # learn CS with lasso constraint
            Ctemp = learnCS(Atemp, C, data, [csFlag, bestParam * 0.1 * coeffSum, maFlag])
            if Ctemp == False:
                print("Could not learn the control strength")
                return
            var, l2 = calcError(data, np.dot(Ctemp, Atemp), False)
            lassoLog.write(str(round(l2, 3)) + "\n")  # error over all data
            lassoLog.close()

        currentVarExplained, currentError = calcError(data, np.dot(Ctemp, Atemp), True)

        # calculate change in error percent
        identical = abs(currentError - prevError) <= errorChangeLimit
        # don't waste time and don't calculate if we don't already have a previous Atemp and Ctemp
        if identical and aProgression and cProgression:
            # make sure that the fold-change for each parameter in the matrix is below the limit
            identical = identical and foldChange(Atemp, aProgression[-1]) <= foldChangeLimit and \
                        foldChange(Ctemp, cProgression[-1]) <= foldChangeLimit
            if identical:
                identicalCounter += 1
                # check if we have reached the proper number of identical iterations
                if identicalCounter == identicalIterations:
                    # adjust the number of iterations and then break out of the loop
                    numIterations = iter + 1
                    break
            else:
                # otherwise, reset the counter
                identicalCounter = 0
        else:
            identicalCounter = 0

        if aProgression and cProgression:
            print('Fold change in A', foldChange(Atemp, aProgression[-1]))
            print('Fold change in C', foldChange(Ctemp, cProgression[-1]))
            print('Fold change in C*A', foldChange(np.dot(Ctemp, Atemp), np.dot(cProgression[-1], aProgression[-1])))
            print('Change in error:', abs(currentError - prevError))

        aProgression.append(Atemp)
        cProgression.append(Ctemp)
        prevError = currentError


        # log the results every 10 iterations
        if iter % 10 == 0:
            saveResults(Ctemp, Atemp, currentVarExplained, "logFiles/csLog" + fileLabel + ".csv",
                        "logFiles/tfaLog" + fileLabel + ".csv", "logFiles/varExplainedLog" + fileLabel + ".csv")

    end = time.time()

    print("\n\n\n")
    print("done learning, now review:")

    for iteration in range(numIterations):
        C = cProgression[iteration]
        A = aProgression[iteration]
        print("iteration ", iteration)
        dataTemp = np.dot(C, A)
        var, l2 = calcError(data, dataTemp, True)

    print("total run time (secs): ", end - start)

    saveResults(Ctemp, Atemp, var, "results/learnedCS" + fileLabel + ".csv", "results/learnedTFA" + fileLabel + ".csv",
                "results/learnedVarExplained" + fileLabel + ".csv")

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
        print("not enough samples to do", folds, "fold cross validation")
        for i in range(folds):
            testCol = random.choice(range(numSamples))
            trainCols = [x for x in range(numSamples) if x != testCol]
            test.append([testCol])
            train.append(trainCols)
    else:
        unsampled = range(numSamples)
        for i in range(folds):
            testCols = random.sample(unsampled, int(len(unsampled) / (folds - i)))
            # print("test", testCols)
            trainCols = [x for x in range(numSamples) if x not in testCols]
            # print("train", trainCols)
            unsampled = [x for x in unsampled if x not in testCols]
            # print("unsampled", unsampled)
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
        # print("test", testCols)
        trainCols = trainColsList[i]
        # print("train", trainCols)
        trainA.append(grabColumns(Amatrix, trainCols))
        testA.append(grabColumns(Amatrix, testCols))
        trainData.append(grabColumns(dataMatrix, trainCols))
        testData.append(grabColumns(dataMatrix, testCols))
    return [trainA, testA, trainData, testData]


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


"""
Input:
  list of input file names
  string for output file labeling
  integer value numIterations for how many iterations of optimization
  list of boolean flags for CS constraints, LASSO contraints, and microarray/RNAseq data
Output:
  the variance explained of the final model learned

Executes learning of activity matrix and baseline values in cs matrix
"""


def tfaInferenceValidation(inputFiles, fileLabel, numIterations, modelParams):
    startFile, csFile, tfaFile, dataFile = inputFiles
    csFlag, lassoFlag, maFlag = modelParams

    # Put model data into matrices and lists
    # this function is in TFAinferenceIO.py
    A = readMatrixFromFile(tfaFile)
    C = readMatrixFromFile(csFile)
    Ctemp = readMatrixFromFile(startFile)
    data = readMatrixFromFile(dataFile)

    # numGenes = len(data)
    # numSamples = len(data[0])
    # numTFs = len(A)

    cProgression = []
    aProgression = []

    start = time.time()

    for iter in range(numIterations):

        print("\niteration ", iter, "\n")

        Atemp = learnTFA(A, Ctemp, data, [])
        if Atemp == False:
            print("Could not learn the activity matrix")
            return

        Ctemp = learnCSValidation(Atemp, C, data, [csFlag, 0, maFlag])

        if Ctemp == False:
            print("Could not learn the control strength")
            return

        aProgression.append(Atemp)
        cProgression.append(Ctemp)

        currentVarExplained, currentError = calcError(data, np.dot(Ctemp, Atemp), True)
        # log the results every 10 iterations
        if iter % 10 == 0:
            saveResults(Ctemp, Atemp, currentVarExplained, "logFiles/csLog" + fileLabel + ".csv",
                        "logFiles/tfaLog" + fileLabel + ".csv", "logFiles/varExplainedLog" + fileLabel + ".csv")

    end = time.time()

    print("\n\n\n")
    print("done learning, now review:")

    for iteration in range(numIterations):
        C = cProgression[iteration]
        A = aProgression[iteration]
        print("iteration ", iteration)
        dataTemp = np.dot(C, A)
        var, l2 = calcError(data, dataTemp, True)

    print("total run time (secs): ", end - start)

    saveResults(Ctemp, Atemp, var, "results/learnedCS" + fileLabel + ".csv", "results/learnedTFA" + fileLabel + ".csv",
                "results/learnedVarExplained" + fileLabel + ".csv")

    return var
