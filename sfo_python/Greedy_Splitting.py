#By Wu Zhenan (bonhomme@nus.edu.sg)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

#Implements Greedy splitting by Zhao et al
#Repeatedly uses Queyranne's algorithm to find the optimal two-partition
#minimizing the energy E(A) + E(V\A)

#E is the energy function per cluster, here assumed to be defined on list only,
#not set.
#V is the list of indices.
#k is the number of clusters.

#Returns partition P. Factor 2 approximation to minimum energy k-partition

import Queyranne as q
def greedy_splitting(E,V,k):
    P = [V]
    for i in range(k-1):
        cand = []
        scores = [float("inf")] * len(P)
        minimum = float("inf")
        argmin = 0
        for j in range(len(P)):
            print("iteration " + str(i+1) + ", cluster "+str(j+1))
            V1 = P[j]
            if len(V1) == 1:
                cand.append([])
                continue
            Etotal = E(list(V1))
            #Now symmetrize the energy function. The line below can be modified
            #if the energy function is valid on sets.
            Fsym = lambda A:E(A) + E(list(V1 - set(A))) - Etotal
            #Find best split using Queyranne's algorithm 
            A1 = set(q.queyranne(Fsym,V1))
            cand.append([A1,V1 - A1])
            scores[j] = Fsym(A1)
            if scores[j] < minimum:
                minimum = scores[j]
                argmin = j
        #Now greedily pick up the best candidate partition. We can make copies  
        #for elements in P instead of implementing the line below for safety.
        P = P[:argmin] + cand[argmin] + P[argmin+1:]
    return P
