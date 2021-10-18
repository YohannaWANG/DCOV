#By Wu Zhenan (bonhomme@nus.edu.sg)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

#opt: option struct of parameters, which may contain the desired tolerance
#for this function.

#F: submodular function
#G: submodular function
#V: index set,not a list or array.
#Input will be cast into one of this form.

#Returns a locally optimal solution to the problem A = argmin_A F(A) - G(A)

import random
import numpy
import Min_Norm_Point as mnp

def ssp(F,G,V,opt={}):
    TOL = opt["ssp_tolerance"] if "ssp_tolerance" in opt else 1e-6
    
    N = len(V)
    pi = random.sample(V,len(V))
    bestVal = float("inf")
    A = set()
    Varray = numpy.array([[x] for x in V])
    while True:
        Hw = ssp_modular_approx(G,pi)
        
        FM = lambda x:F(x) - sum([Hw[y-1] for y in x])
        A = mnp.min_norm_point(FM,Varray,{"minnorm_init":A},0,set)[0]
        curVal = FM(A)
        D = V - A
        D = random.sample(D,len(D))        
        pi = list(A) + D
        if curVal < bestVal - TOL:
            bestVal = curVal
        else:
            break
    return A

def ssp_modular_approx(G,pi):
    H = [0] * max(pi)
    W = []
    oldVal = G(W)
    for i in range(len(pi)):
        W += [pi[i]]
        newVal = G(W)
        H[i] = newVal - oldVal
        oldVal = newVal
    return H
