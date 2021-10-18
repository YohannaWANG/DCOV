#By Wu Zhenan (bonhomme@nus.edu.sg)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

import numpy
#Iwata function
def Iwata(n):
    return lambda x:len(x) * (n - len(x))-sum(x) * 5 + 2 * n * len(x)

#Cut function that returns the sum of weights of paths between two components in
#a graph represented by a matrix.
def cutfun(graph):
    def f(component):
        if type(component) !=type(set()):
            component = set(component)            
        outside_component = set(range(1,len(graph)+1)) - component
        return sum([sum([graph[x - 1][y
                    - 1] for y in outside_component]) for x in component])
    return f

#lincomb is defined below but not used in SSP, unlike in the matlab codes by
#Andreas Krause.

#Returns a linear combination of functions.
def lincomb(list_of_functions,list_of_weights):
    return lambda x:sum([list_of_functions[y](x)
                * list_of_weights for y in range(len(list_of_functions))])

#The function below computes log of determinant of a matrix using base 2.
def logdet_by_cholesky(matrix):
    c = numpy.linalg.cholesky(matrix)
    return 2 * numpy.sum(numpy.log2(numpy.diagonal(c)))

#The function below computes the Gaussian mutual information between a set and
#its complement.

#sigma: covariance matrix in the form of a numpy array
#sset: set of rows to be separated with others
def fn_mi(sigma,sset,sset_list=None):
    if len(sset) == 0 or len(sset) == len(sigma):
        return 0
    if not sset_list:
        sset_list = list(sset)
    Ac = {x for x in range(1,len(sigma)+1)} - sset
    Ac_list = list(Ac)
    for x in range(len(sset_list)):
        sset_list[x] -= 1
    for x in range(len(Ac_list)):
        Ac_list[x] -= 1
    sigmaA = sigma[numpy.ix_(sset_list,sset_list)]
    #Somehow += does not work here.
    sigmaA = sigmaA + 1e-10 * numpy.identity(len(sigmaA))
    
    sigmaAcomp = sigma[numpy.ix_(Ac_list,Ac_list)]
    sigmaAcomp = sigmaAcomp + 1e-10 * numpy.identity(len(sigmaAcomp))
    sigmaAAcomp = sigma[numpy.ix_(sset_list,Ac_list)]

    sigmaAcond = sigmaA - numpy.matmul(sigmaAAcomp,numpy.matmul(
        numpy.linalg.inv(sigmaAcomp),numpy.transpose(sigmaAAcomp)))
    return 0.5 * (logdet_by_cholesky(sigmaA) - logdet_by_cholesky(sigmaAcond))

#The functions below are an entropy function and an infogain function.

#sigma is the covariance matrix in the form of a numpy array.
#sset contains no duplicate.
#import numpy
def entropy(sigma,sset):
    if not sset:
        return 0
    sset = [x - 1 for x in sset]
    cholA = numpy.linalg.cholesky(sigma[numpy.ix_(sset,sset)])
    # Somehow += does not work here.
    cholA = cholA + (1e-10) * numpy.identity(len(cholA))
    s = sigma[numpy.ix_(sset,sset)]
    s = s + (1e-10) * numpy.identity(len(s))
    return 0.5 * numpy.log2((2 * numpy.pi * numpy.e) ** len(
        cholA)) + sum(numpy.log2(numpy.diag(cholA)))

def infogain(sigma,sset,noise):
    if not sset:
        return 0
    cholA = numpy.linalg.cholesky(sigma[numpy.ix_(sset,sset)])
    #Somehow += does not work here.
    cholA = cholA + (1e-10) * numpy.identity(len(cholA))
    return 0.5 * numpy.log2((2 * numpy.pi * numpy.e) ** len(cholA)) + sum(
        numpy.log2(numpy.diag(cholA))) - numpy.log2(numpy.sqrt(2 * numpy.pi * 
                                                               numpy.e * noise)) * len(sset)


#Below is an energy function for ising model for image denoising

#img: n x m binary array in the form of a list of list (image)
#A: subset of the pixels set to 1 (ranging in 1:(n*m), ordered from top to bottom then left to
#right)
#coeffPix/H/V/Diag: negative log potentials for differing pixels, horizontal, vertical or 
#diagonal intensity mismatch

def ising(img,coeffPix,coeffH,coeffV,coeffD):
	val0 = evalIsing(img,{},coeffPix,coeffH,coeffV,coeffD)
	return lambda A:evalIsing(img,A,coeffPix,coeffH,coeffV,coeffD) - val0

def evalIsing(img,A,coeffPix,coeffH,coeffV,coeffD):
	if type(A) != type(set()):
		A = set(A)
	r,c = len(img),len(img[0])
	A = [[(x-1)%r for x in A],[(x-1)//r for x in A]]
	mask = numpy.zeros((r,c))
	mask[A[0],A[1]] = 1
	delta = abs(mask - img)
	Epix = numpy.sum(delta)#Energy through pixelwise disagreement
	Ehor = numpy.sum(abs(mask[:,1:] - mask[:,:c-1]))#Energy through "horizontal" pixel 
	#differences
	Evert = numpy.sum(abs(mask[1:,:] - mask[:r-1,:]))#Energy through "vertical" pixel
	#differences
	Ediag = numpy.sum(abs(mask[:r-1,:c-1] - mask[1:,1:]) +
                          abs(mask[1:,:c-1] - mask[:r-1,1:]))
	return coeffPix * Epix + coeffH * Ehor + coeffV * Evert + coeffD * Ediag


#Given a submodular function G, the function below represents the "inverse".
#F(A) = G(V\A)
#V: index set

def invert(G,V):
    return lambda A:G(V-A)

    
