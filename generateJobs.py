"""
In the directory that you run this, 
you need three subdirectories called 
	randStarts
	results
	logfiles
All input/output data files should be in csv format
The only exception is the lasso parameter log files, which are in tsv format
"""



def generateJobFiles(inputFiles, fileLabels, iterations, modelParams, sbatchParams):
	numStarts, parallel, mem = sbatchParams
	for i in range(len(fileLabels)):
		jobFile = open(fileLabels[i]+'.py', 'w')
		jobFile.write('from TFAinference import *\n')
		jobFile.write('i = int(sys.argv[1])\n')
		jobFile.write('randStartFile = "randStarts/randCS" + str(i) + ".csv"\n')
		jobFile.write('inputFiles = [randStartFile] + ' + str(inputFiles) + '\n')
		jobFile.write('try:\n open(randStartFile)\n')
		jobFile.write('except:\n print "cannot find file", randStartFile\n')
		jobFile.write('var = tfaInference(inputFiles, "' + fileLabels[i] 
				+ '" + str(i), ' + str(iterations) + ", " + str(modelParams[i]) + ")\n")
		jobFile.close()
		
		sbatchFile = open('../'+fileLabels[i]+'.sbatch', 'w')
		sbatchFile.write('#!/bin/bash\n')
		sbatchFile.write('#SBATCH --array=1-'+str(numStarts)+'%'+str(parallel)+'\n')
		sbatchFile.write('#SBATCH --mem='+str(mem)+'G\n')
		sbatchFile.write('#SBATCH -o slurmOut/'+fileLabels[i]+'_%A_%a.out\n')
		sbatchFile.write('#SBATCH -e slurmOut/'+fileLabels[i]+'_%A_%a.err\n')
		sbatchFile.write('module load gurobi\n')
		sbatchFile.write('ID=${SLURM_ARRAY_TASK_ID}\n')
		sbatchFile.write('gurobi.sh TFAinference/'+fileLabels[i]+'.py ${ID}\n')

