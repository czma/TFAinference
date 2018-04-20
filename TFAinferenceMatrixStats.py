from TFAinferenceMatrixMath import calcError, foldChange, computeCSPseudocount
from TFAinferenceIO import readMatrixFromFile
import argparse
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt


"""
These functions compute various data for successive iterations of a SINGLE RANDOM START
"""


'''
Computes fold changes with various pseudocounts

Param:
    matricies: list of matrixes
    lower_limit: lower limit of pseudocount
    upper_limit: upper limit of pseudocount
    step = step size for pseudocount
'''
def computeFoldChangesWithPseudocountRange(matricies, lower_limit=0, upper_limit=1, step=0.1):
    results = [] # list to hold results, format (pseudocount, [foldchange across subsequent matricies])
    for ps in np.arange(lower_limit, upper_limit, step):
        r = [] # intermediate results
        for i in range(0, len(matricies) - 1):
            r.append(foldChange(matricies[i], matricies[i+1], l=ps))
        results.append((ps, r))
    return results

'''
Computes fold changes with various pseudocounts for ONE matrix pair

Param:
    A: matrix
    B: matrix
    lower_limit: lower limit of pseudocount
    upper_limit: upper limit of pseudocount
    step = step size for pseudocount
'''
def computeFoldChangesWithPseudocountRange_SINGLE(A, B, lower_limit=0, upper_limit=1, step=0.1, method='max'):
    r = []
    for ps in np.arange(lower_limit, upper_limit, step):
        r.append((ps, foldChange(A, B, l=ps, method=method)))
    return r

'''
Computes changes in variance explained for successive iterations

Param:
    varsExplained: list of variances explained
'''
def computeChangesInVarianceExplained(varsExplained):
    results = [] # list to hold results, format (pseudocount, [foldchange across subsequent matricies])
    for i in range(0, len(varsExplained) - 1):
        results.append(varsExplained[i+1] - varsExplained[i])
    return results

'''
Computes normalized change in error for successive iterations

Param:
    errors: list of errors
'''
def computeNormalizedChangesInError(errors):
    results = [] # list to hold results, format (pseudocount, [foldchange across subsequent matricies])
    for i in range(0, len(errors) - 1):
        results.append((errors[i+1] - errors[i]) / ((errors[i+1] + errors[i]) / 2.0))
    return results

"""
Computes various data across results of ALL RANDOM RESTARTS
"""

'''
Finds all results within the bound of the optimal solution

Param:
    varsExplained: list of variances explained
    varianceExplainedBound: numerical bound for which to return solutions
'''
def computeResultsWithinVarianceBoundOfOptimalSolution(varsExplained, varianceExplainedBound = 0.01):
    v = list(enumerate(varsExplained))
    # v.sort(key=lambda x: x[1])
    opt = max(v, key= lambda x: x[1])
    return [i for i,x in filter(lambda x: abs(x[1] - opt[1]) <= varianceExplainedBound, v)], opt[0]


'''
Finds all results within the bound of the optimal solution

Param:
    errors: list of errors
    errorChangeBound: numerical bound for which to return solutions
'''
def computeResultsWithinNormalizedErrorChangeBoundOfOptimalSolution(errors, errorChangeBound = 0.01):
    v = list(enumerate(computeNormalizedChangesInError(errors)))
    # v.sort(key=lambda x: x[1])
    opt = min(errors)
    return [i for i,x in filter(lambda x: abs(x - opt) <= errorChangeBound, v)]



def computeNumberOfResultsWithinVarianceBoundRangeOfOptimalSolution(varsExplained, lower_limit=0.0, upper_limit=0.5, step=0.05):
    return [(bound, len(computeResultsWithinVarianceBoundOfOptimalSolution(varsExplained, varianceExplainedBound=bound))) for bound in np.arange(lower_limit, upper_limit, step)]


def zip_fill(*argv, fill=None):
    argv_s = [[i, list(l)] for i, l in sorted(enumerate(argv), key=lambda x: len(x[1]))]
    for l in argv_s:
        l[1] += [fill] * (len(argv_s[-1][1]) - len(l[1]))
    argv_c = [l for i, l in sorted(argv_s, key=lambda x: x[0])]
    for i in range(len(argv_s[0][1])):
        yield tuple(x[i] for x in argv_c)  # preserve original order


def read_all_matricies(filename, iterations):
    data = readMatrixFromFile(filename)
    if not (len(data) / iterations) % 1 == 0.0:
        raise ValueError('Data length does not match iterations!')
    mat_length = len(data) / iterations
    return [data[mat_length*i:mat_length*(i+1)] for i in range(iterations)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process matricies to generate plots')
    parser.add_argument('--cs', '-c', type=str, nargs='+', help='List of control strength matrix files', required=True)
    parser.add_argument('--tfa', '-t', type=str, nargs='+', help='List of transcription factor activity matricies', required=True)
    parser.add_argument('--variances', '-v', type=str, nargs='+', help='List of variance explained files', required=True)
    parser.add_argument('--limit', '-l', type=float, default=1000, help='Upper limit of graph')
    parser.add_argument('--output', '-o', type=str, default='matrix_stats', help='Output filename base')
    # parser.add_argument('--error', '-e', action='store_true', help='Output error plots')
    parser.add_argument('--iterations', '-i', type=int, action='store', default=100, help='Number of iterations in log files')
    parser .add_argument('--results', '-r', action='store_true', help='Using final results from random starts')
    parser .add_argument('--computePseudocounts', action='store_true', help='Use computed pseudocounts for CS matrix')

    args = parser.parse_args()

    y_limit = [0, args.limit]

    cs_matricies = []
    tfa_matricies = []
    varsExplained = []

    # final results, with one matrix per file
    if args.results:
        for file in args.cs:
            cs_matricies.append(readMatrixFromFile(file))
        for file in args.tfa:
            tfa_matricies.append(readMatrixFromFile(file))
        for file in args.variances:
            with open(file) as f:
                line = f.readline().strip()
                varsExplained.append(float(line))


        ####
        # Plot the number of results within the range of the variance explained of the optimal solution
        ####
        v = computeNumberOfResultsWithinVarianceBoundRangeOfOptimalSolution(varsExplained)
        x,y = [bound for bound, num in v], [num for bound, num in v]
        plt.plot(x, y, 'b-')
        plt.title('Number of Results Within Variance Explained Bound of Optimal Solution')
        plt.xlabel('Variance Explained Bound')
        plt.plot()
        plt.savefig(args.output + '_numResultsInVarRange.png')
        plt.clf()


        ####
        # Plot the fold change of the results within the variance explained bound
        ####
        v, opt = computeResultsWithinVarianceBoundOfOptimalSolution(varsExplained)
        data = [foldChange(cs_matricies[x], cs_matricies[opt], l=0.001) for x in v if not x == opt]
        plt.hist(data)
        plt.title('Fold Change of Results Within Variance Explained Bound of Optimal Solution in CS matricies')
        plt.xlabel('Fold Change')
        plt.ylabel('Frequency')
        plt.plot()
        plt.savefig(args.output + '_foldChangeResultsInVarRange_CS.png')
        plt.clf()
        data = [foldChange(tfa_matricies[x], tfa_matricies[opt], l=0.001) for x in v if not x == opt]
        plt.hist(data)
        plt.title('Fold Change of Results Within Variance Explained Bound of Optimal Solution in TFA matricies')
        plt.xlabel('Fold Change')
        plt.ylabel('Frequency')
        plt.plot()
        plt.savefig(args.output + '_foldChangeResultsInVarRange_TFA.png')
        plt.clf()



    # intermediate results, with multiple matricies per file
    else:
        if args.computePseudocounts:
            for file in args.cs:
                cs_matricies.append(read_all_matricies(file, args.iterations))
                r = []  # intermediate results
                for i in range(0, len(cs_matricies[-1]) - 1):
                    r.append(foldChange(cs_matricies[-1][i], cs_matricies[-1][i+1], l=computeCSPseudocount(cs_matricies[-1][i])))
                plt.plot([x + 1 for x in range(1, len(r) + 1)], r, '-')
                plt.title('Fold-change in CS for computed pseudocounts (Limited to {})'.format(str(y_limit)))
                plt.xlabel('Iteration')
                plt.ylim(y_limit)
                plt.savefig(file + '_foldChangesComputedPseudocounts.png')
                plt.clf()
        else:
            for file in args.cs:
                cs_matricies.append(read_all_matricies(file, args.iterations))
                ps_foldchanges = computeFoldChangesWithPseudocountRange(cs_matricies[-1])
                for ps, foldchanges in ps_foldchanges:
                    plt.plot([x+1 for x in range(1, len(foldchanges) + 1)], foldchanges, '-', label=str(ps))
                plt.title('Fold-change in CS for differing pseudocounts (Limited to {})'.format(str(y_limit)))
                plt.xlabel('Iteration')
                plt.legend()
                plt.ylim(y_limit)
                plt.savefig(file + '_foldChangesPseudocounts.png')
                plt.clf()
        for file in args.tfa:
            tfa_matricies.append(read_all_matricies(file, args.iterations))
            ps_foldchanges = computeFoldChangesWithPseudocountRange(tfa_matricies[-1])
            for ps, foldchanges in ps_foldchanges:
                plt.plot([x + 1 for x in range(1, len(foldchanges) + 1)], foldchanges, '-', label=str(ps))
            plt.title('Fold-change in TFA for differing pseudocounts (Limited to {})'.format(str(y_limit)))
            plt.xlabel('Iteration')
            plt.legend()
            plt.ylim(y_limit)
            plt.savefig(file + '_foldChangesPseudocounts.png')
            plt.clf()
        for file in args.variances:
            vars_tmp = []
            with open(file) as f:
                for line in f:
                    line = line.strip()
                    vars_tmp.append(float(line))
            varsExplained.append(vars_tmp)
            vars_diffs = computeChangesInVarianceExplained(varsExplained[-1])
            plt.plot([x+1 for x in range(1, len(vars_diffs) + 1)], vars_diffs, 'r-', label='Variance Explained Differences')
            plt.title('Change in Variance Explained from Previous Iteration')
            plt.xlabel('Iteration')
            plt.savefig(file + '_changesInVarianceExplained.png')
            plt.clf()

        ####
        # Plot the number of matricies that are already in their final positions in terms of variance explained
        ####
        numbered_variances = list(enumerate(varsExplained))
        final = sorted(numbered_variances, key=lambda x:x[1][-1])
        results = []

        for i in range(args.iterations):
            r = 0
            tmp = sorted(numbered_variances, key=lambda x:x[1][i])
            for (final_pos, final_vars),(pos, vars) in zip(final,tmp):
                if final_pos == pos:
                    r += 1
            results += r

        plt.plot([x + 1 for x in range(1, len(results) + 1)], results, 'r-')
        plt.title('Rank-Order Correlations between Variance Explained Rank after n Iterations')
        plt.xlabel('Iteration')
        plt.savefig(args.output + '_varExplainedRankCorrelation.png')
        plt.clf()
