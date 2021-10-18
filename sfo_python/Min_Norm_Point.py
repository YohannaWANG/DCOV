#By Wu Zhenan (bonhomme@u.nus.edu)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

import numpy
#F: submodular function that also takes a list as a valid input
#V: index list which is a numpy array in the form [[number],[number],...]
#w: weight vector which is a numpy array in the same form as V
def polyhedrongreedy(F,V,w):
    n = len(V)
    v = [[w[x][0],V[x][0],x] for x in range(n)]
    v.sort(reverse=True)
    r = numpy.zeros((n,1))
    f = set()
    s=F(f)
    for x in v:
        f.add(x[1])
        g = F(f)
        r[x[2]] = [g-s]
        s = g
    return r

#charvector indicates whether an element in list V is in set A, and True is
#converted to 1, False to 0
def charvector(V,A):
    r = []
    for x in V:
        r.append([1 if x[0] in A else 0])
    return numpy.array(r)

#opt: option struct of parameters. referencing:

#minnorm_init: starting guess for optimal solution
#minnorm_stopping_thresh: stopping threshold for search
#minnorm_tolerance: numerical tolerance
#minnorm_callback: callback routine for visualizing intermediate solutions

#Returns an optimal solution A, bound on suboptimality
#The variable below stores the class type function which will be used later
f = type(lambda x:x)

def min_norm_point(F,V,opt={},display=0,return_type=list):
    n = len(V)
    Ainit = opt["minnorm_init"] if "minnorm_init" in opt else set()
    eps = opt["minnorm_stopping_thresh"] if "minnorm_stopping_thresh" in opt else 1e-10
    TOL = opt["minnorm_tolerance"] if "minnorm_tolerance" in opt else 1e-10
    
    #Step 1: Initialize by picking a point in the polytope.
    wA = charvector(V,Ainit)
    xw = polyhedrongreedy(F,V,wA)
    xhat = xw
    
    #For now I think inheriting xhat is OK as changing xhat completely to
    #another array should not affect S.
    S = xhat
    
    Abest = -1
    Fbest = numpy.inf
    c=1
    while True:
        #Step 2: Find a vertex in the base polytope that minimizes its
        #inner-product with xhat.
        squarednorm = sum(numpy.multiply(xhat,xhat))
        norm = numpy.sqrt(squarednorm)
        if norm < TOL: #Snap to zero
            xhat = numpy.zeros((n,1))
        #Get phat by going from xhat towards the origin until we hit the
        #boundary of the base polytope
        phat = polyhedrongreedy(F,V,-xhat)
        S = numpy.append(S,phat,axis=1)
        #Check current function value
        lessthanzero = xhat < 0
        s = V[lessthanzero]
        Fcur = F(s)
        if Fcur < Fbest:
            Fbest = Fcur
            Abest = s
            
        #Get suboptimality bound
        subopt = Fbest - sum(xhat[lessthanzero])
        if display:
            print("suboptimality bound: " + str(float(Fbest - subopt)) +
                  "<=min_A F(A)<=F(A_best)=" + str(float(Fbest)) +
                  "; delta<=" + str(float(subopt)))
            
        absolute = abs(sum(numpy.multiply(xhat,xhat-phat))[0])
        if absolute < TOL or subopt<eps:
            #We are done: xhat is already closest norm point
            if absolute < TOL:
                subopt = 0
            A = Abest
            break
        
        #Here's some code just for outputting the current state
        if "minnorm_callback" in opt and type(opt["minnorm"])==f:#Do something
            #with current state
            opt["minnorm_callback"](Abest)

        [xhat,S] = min_norm_point_update(xhat,S,TOL)
    if display:
        print("suboptimality bound: " + str(float(Fbest - subopt)) +
              "<=min_A F(A)<=F(A_best)=" + str(float(Fbest)) +
              "; delta<=" + str(float(subopt)))
    #Casts return value into the type we want
    if return_type:
        A=return_type(A)
    return [A,subopt]

#Helper function for updating xhat
def min_norm_point_update(xhat,S,TOL):
    while True:
        #Step 3: Find minimum norm point in affine hull spanned by S

        S0 = S[:,1:] - S[:,[0]*(len(S[0]) - 1)]
        firstcolumn = S[:,[0]]
        y = firstcolumn - numpy.matmul(numpy.matmul(S0,numpy.linalg.pinv(S0)),
                                       firstcolumn)#Now y in min norm

        #Get representation of y in terms of S. Enforce affine combination
        #(i.e. sum(mu)==1)
        pseudoinverse = numpy.linalg.pinv(numpy.append(S,
                                                    [[1] * len(S[0])],axis=0))
        mu = numpy.matmul(pseudoinverse,numpy.append(y,[[1]],axis=0))
        
        #y is written as positive convex combination of S <==> y in conv(S)
        if not numpy.size(mu[mu < -TOL]) and abs(sum(mu) - 1)<TOL:
            return [y,S]

        #Step 4: Project y back into polytope.

        #Get representation of xhat in terms of S;enforce that we get affine
        #combination (i.e. sum(L)==1)
        L = numpy.matmul(pseudoinverse,numpy.append(xhat,numpy.array([[1]]),
                                                    axis=0))

        #Now find z in conv(S) that is closest to y
        bounds = numpy.divide(L,L-mu)
        bounds = bounds[bounds>TOL]
        beta = min(bounds)
        z = (1 - beta) * xhat + beta * y
        gamma=(1 - beta) * L + beta * mu
        S=S[:,numpy.where(gamma > TOL)[0]]
        xhat = z
