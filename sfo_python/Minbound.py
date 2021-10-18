#By Wu Zhenan (bonhomme@nus.edu.sg)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

#Get bound on suboptimality/certificate of optimality [Edmonds '71]

#F: Submodular function
#V: index list in the form [[number],[number],...] 
#A: set to be tested
#Returns bound <= min_A F(A)

import Min_Norm_Point as mnp
def minbound(F,V,A):
    w = mnp.charvector(V,A)
    xw = mnp.polyhedrongreedy(F,V,w)
    return sum(xw[xw<0])
