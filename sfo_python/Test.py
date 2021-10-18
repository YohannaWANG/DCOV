#By Wu Zhenan (bonhomme@nus.edu.sg)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

import numpy
import Miscellaneous_Functions as m
def test_min_norm_point():
    import Min_Norm_Point as mnp
    Iwata = m.Iwata
    A = mnp.min_norm_point(Iwata(100),numpy.array([[x] for x in range(1,101)]))
    print(A)
#Optimal Solution Selects 34 to 100
    
def test_queyranne():
    import Queyranne as q
    G_un = [[0,1,1,0,0,0],[1,0,1,0,0,0],[1,1,0,1,0,0],[0,0,1,0,1,1]
            ,[0,0,0,1,0,1],[0,0,0,1,1,0]]
    V_G=[x for x in range(1,len(G_un)+1)]
    cutfun = m.cutfun
    print(q.queyranne(cutfun(G_un),V_G))
#Best cut should be [4,5,6]

#The two tests below are on SSP.
#The result of the test below is very off.
def test_SSP_by_fn_mi():
    import csv
    import SSP
    with open("Merced_Data_Sigma.csv") as data:
        r = csv.reader(data)
        next(r)
        s = []
        for row in r:
            s.append(row)
    for x in s:
        for y in range(len(x)):
            x[y] = float(x[y])
    s = numpy.array(s)
    F = lambda x:len(x)
    G = lambda x:m.fn_mi(s,set(x),x)
    print(SSP.ssp(F,G,{x for x in range(1,87)}))
#Matlab code chooses 2, 3, 8, 24, 36, 42, 85. The test above selects many more
#and often does not include all the correct answers.

#This is the test of SSP with Iwata function. 
def test_SSP_by_Iwata():
    import SSP
    F = lambda x:len(x)
    G = m.Iwata(100)
    print(SSP.ssp(F,G,set(range(1,101))))
#Python seems to select smaller sets than Matlab.

def test_S_T_Mincut():
    G_dir=[[0,1,1.2,0,0,0],[1.3,0,1.4,0,0,0],[1.5,1.6,0,1.7,0,0],
           [0,0,1.8,0,1.9,2],[0,0,0,2.1,0,2.2],[0,0,0,2.3,2.4,0]]
    cutfun = m.cutfun
    import S_T_Mincut as stmc
    print(stmc.s_t_mincut(cutfun(G_dir),set(range(1,7)),1,6))
    print(stmc.s_t_mincut(cutfun(G_dir),set(range(1,7)),4,6))
#The output for the first line is [2,3,1] and for the second is [1,2,3,4]

#There is no test for Minbound package as it is simple.


#Below are two tests for Greedy_Splitting. For the test with entropy, Python
#codes often fail when doing Cholesky decomposition.

#Further down below is a test using Iwata function. This test result coincides
#with that by matlab.
def test_Greedy_Splitting_by_entropy():
    import Greedy_Splitting as gs
    entropy = m.entropy
    C1 = numpy.random.normal(5,1,16)
    C2 = numpy.random.normal(0,1,16)
    C3 = numpy.random.normal(-5,1,16)
    C = [C1,C2,C3]
    D = [[0] * 24] * 24
    #Calculation of distance
    for x in range(24):
        c = 2 * (x%8)
        x_coordinates = C[x//8][c : c+2]
        for y in range(x):
            c = 2 * (y%8)
            y_coordinates = C[y//8][c : c+2]
            D[x][y] = D[y][x] = ((y_coordinates[0] - x_coordinates[0]) **2 + (
                y_coordinates[1] - x_coordinates[1]) **2) / (-4)
    sigma_cl = numpy.e ** numpy.array(D)
    sigma_cl = sigma_cl + 0.01 * numpy.identity(24)
    print(gs.greedy_splitting(lambda x:entropy(sigma_cl,x),
                              set(range(1,25)),2))

def test_Greedy_Splitting_by_Iwata():
    import Greedy_Splitting as gs
    print(gs.greedy_splitting(m.Iwata(24),set(range(1,25)),2))
    print(gs.greedy_splitting(m.Iwata(24),set(range(1,25)),3))
#This test result coincides with that by matlab, making 24 and other indices as
#two clusters and 24, 23 and others as three clusters.

def test_Our_Experiment():
    import Our_Experiment as oe
    print(oe.our_experiment([[1,0.9], [0.9,1]]))
#The result is [[2],1].
    
