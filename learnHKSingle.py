"""
In the directory that you run this, 
you need three subdirectories called 
	randStarts
	results
	logfiles
All input/output data files should be in csv format
The only exception is the lasso parameter log files, which are in tsv format
"""

from TFAinference import *

"""
input file names
	binaryCSFilename should be a file with 1, -1, and 0 values 
		organized like a genes x TFs matrix
		each row is a new line, with columns separated by commas
	binaryTFAFilename should be a file with 1 and 0 values 
		organized like a TFs x samples matrix
		each row is a new line, with columns separated by commas
	matrixEFilename should be a file with the logFC values from the microarray data
		organized like a genes x samples matrix
		each row is a new line, with columns separated by commas
"""
binaryCSFilename = "signedBinaryCS.csv"
binaryTFAFilename = "binaryTFASmall.csv"
matrixEFilename = "matrixESmall.csv"

#how long to run
iterations = 3

#an added tag to output file names to identify and group them
fileLabel = "HKwoLASSO"

"""
parameters for things that depends on input data, or model changes being tested
currently:
	the first boolean is whether or not to learn with CS constraints
	the second boolean is whether or not to learn with LASSO constraint
	the third boolean is whether or not the data is microarray (as opposed to RNAseq)
	the fourth boolean is whether or not we are validating over a known CS matrix
	the fifth is the number of iterations after which we can learn 0 coeff values
"""
modelParams = [True, False, True, False, iterations]

def learnSingleRandomStart():
		#which number random start this run is
		i = int(sys.argv[1])
		#checking for random start file
		randStartFile = "randStarts/randCS"+str(i)+".csv"
		try:
			open(randStartFile)
		except:
			#this function is in TFAinferenceIO.py
			makeRandStart(binaryCSFilename, randStartFile, modelParams[0])

		inputFiles = [randStartFile, binaryCSFilename, binaryTFAFilename, matrixEFilename]
		#this function is in TFAinference.py
		#returns the variance learned
		var = tfaInference(inputFiles, fileLabel+str(i), iterations, modelParams)

learnSingleRandomStart()
