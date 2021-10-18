#By Wu Zhenan (bonhomme@nus.edu.sg)
#Modified from the matlab codes by Andreas Krause (krausea@gmail.com)

#Finds the marginal gain in function value when list B integrates list x.
#[number],...] 
def residual(F,B):
    def f(x):
        if type(x) != type(set):
            x = set(x)
        return F(x.union(B)) - F(B)
    return f

#Finding the minimum A of a submodular function such that s in A and t not in A

#Returns A in mincut(F,V,s,t)
#F: submodular function
#V: index set
#s: element in V to include
#t: element in V to exclude
#opt(optional): options, depending on
#minnorm_stopping_thresh: threshold for stopping minimization(1e-5)

import Min_Norm_Point as mnp
def s_t_mincut(F,V,s,t,opt={}):
    if "minnorm_stopping_thresh" not in opt:
        opt["minnorm_stopping_thresh"] = 1e-5
    F2 = residual(F,{s})
    V2 = mnp.numpy.array([[x] for x in V - {s,t}])
    return  mnp.min_norm_point(F2,V2,opt,return_type=list)[0]+[s]
