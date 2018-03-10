from generateJobs import *

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
iterations = 100

#numStarts, parallel, mem
sbatchParams = [100, 30, 16]

#an added tag to output file names to identify and group them
fileLabels = ["learnNNNN", "learnLNNN", "learnNSNN", "learnNNNZ5", "learnLSNN", "learnLNNZ5", "learnNSNZ5", "learnLSNZ5"]

"""
parameters for things that depends on input data, or model changes being tested
currently:
	the first boolean is whether or not to learn with CS constraints
	the second boolean is whether or not to learn with LASSO constraint
	the third boolean is whether or not the data is microarray (as opposed to RNAseq)
	the fourth boolean is whether or not we are validating over a known CS matrix
	the fifth is the number of iterations after which we can learn 0 coeff values
"""
modelParams = [ [ True, False, True, False, iterations+1], \
		[ True,  True, True, False, iterations+1], \
		[False, False, True, False, iterations+1], \
		[ True, False, True, False, 	       5], \
		[False,  True, True, False, iterations+1], \
		[ True,  True, True, False, 	       5], \
		[False, False, True, False, 	       5], \
		[False,  True, True, False, 	       5]]


generateJobFiles([binaryCSFilename, binaryTFAFilename, matrixEFilename], fileLabels, iterations, modelParams, sbatchParams)
