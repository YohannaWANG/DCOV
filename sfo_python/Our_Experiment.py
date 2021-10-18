from sfo_python import Min_Norm_Point as mnp
from sfo_python import Miscellaneous_Functions 
#from Miscellaneous_Functions import logdet_by_cholesky
import numpy as np

#Sigma is a square matrix is the form of a normal list or a numpy array.
def fn_logdet(matrix,exception):
    if type(matrix) != type(np.array([])):
        matrix = np.array(matrix)
    def cholesky(sset):
        sset = [x-1 for x in sset]
        sset.append(exception-1)
        selected = matrix[np.ix_(sset,sset)]
        if len(sset) != 1:
            selected += 1e-10 * np.identity(len(sset))
        return Miscellaneous_Functions.logdet_by_cholesky(selected)
    return cholesky

def our_experiment(Sigma):
    n = len(Sigma)
    V = np.array([[x] for x in range(1,n+1)])
    Abest = []
    ibest = 0
    bestval = 0
    for i in range(1,n+1):
        Vi = V[V!=i]
        Vi = np.array([[x] for x in Vi])
        F = fn_logdet(Sigma,i)
        A = mnp.min_norm_point(F,Vi)[0]
        if i == 1 or F(A) < bestval:
            Abest = A
            ibest = [i]
            bestval = F(A)
    return np.array([Abest,ibest]).flatten() 

print(our_experiment([[1,0.9], [0.9,1]]))
