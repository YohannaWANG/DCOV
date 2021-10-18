#By Wu Zhenan (bonhomme@nus.edu.sg)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

#F: submodular function.
#V: index list or a set.

#Returns an optimal solution to min F(A) s.t. 0<|A|<n

from collections import OrderedDict
def queyranne(F,V):
    def Fnew(a):
        r=[]
        for x in a:
            r += S[x - 1]
        return F(r)
    n = len(V)
    S = [[x] for x in V]
    s = []
    A = []
    inew = OrderedDict()
    for x in range(1,n+1):
        inew[x] = x
    minimum=float("inf")
    position_of_min=0
    for h in range(n-1):
        #Find a pendant pair
        [t,u] = pendentpair(Fnew,inew)

        #This gives a candidate solution
        A.append(S[u - 1].copy())
        s.append(Fnew({u}))
        if s[-1] < minimum:
            minimum = s[-1]
            position_of_min = len(s) - 1
        S[t - 1] += S[u - 1]
        del inew[u]
        for x in range(len(S[u - 1])):
            S[u - 1][x] *= -1
    return A[position_of_min]

#Implements the pendant pair finding subroutine of Queyranne's algorithm
#(Queyranne '95)
#F is the submodular function
#inds is an array of indices; (typically, 1:n)

def pendentpair(F,V):
    vstart = V.popitem(last=False)[0]
    vnew = vstart
    n = len(V)
    Wi = []
    used = [0] * n
    for i in range(n):
        vold = vnew
        Wi += [vold]
        #Now update the keys
        keys = [1e99] * n
        minimum = float("inf")
        counter = -1
        for j in V:
            counter += 1
            if used[counter]:
                continue
            Wi += [V[j]]
            keys[counter] = F(Wi) - F({V[j]})                
            del Wi[-1]
            if keys[counter] < minimum:
                minimum = keys[counter]
                argmin_key = j
                argmin_position = counter
            vnew = argmin_key
            used[argmin_position] = 1
    V[vstart] = vstart
    V.move_to_end(vstart,last=False)
    return [vold,vnew]
