"""
    Title:  Identifiability of chain graph models;
    Author: 
    Email: 
    Date:   2021-02-02
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import STAP

import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri

import data as data
from utils import args
from evaluate import count_accuracy
#from baseline import baseline_cg, baseline_cg_lcdlike


pandas2ri.activate()
rpy2.robjects.numpy2ri.activate()


packageNames = ('graph', 'RBGL', 'ggm', 'mgcv', 'pcalg')
utils = rpackages.importr('utils')
#utils.chooseCRANmirror(ind=1)

#packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]

"""
    Notes: Uncomment this part if you don't have the following packages installed on ur comp.
"""
# Running R in Python example installing packages:
#if len(packnames_to_install) > 0:
#    utils.install_packages(StrVector(packnames_to_install))

importr('mgcv')
importr('pcalg')
#importr('lcd')
robjects.r.source("utils.R")


def load_graph_data():
    """
    Function: get parameters for chain graph simulation data generation
    Param:
        n = 5  # number of nodes
        c = 2  # number of chain components
        a = 2  # average node degree
        d = 10000  # number of dimensions
    """    
    n, d = args.n, args.d
    
    if args.graph == "chain_graph":

        det, c, a = args.det, args.c, args.a
        G, G_undirect, D, U, sorted_D, sorted_U = data.gen_cg_adj(n, c, a, d)
        cc = data.connected_comp(n, G_undirect)

        membership, c_size, num = data.is_chain_graph(cc)

        simu, G, top_sort  = test_chain_graph(G, membership)
        if args.plot:    
            data.plot_chain_graph(cc, D, U)    
            
        if data.is_acyclic(G) == 1:
            print("Contains cycles in this graph.")
            
    if args.graph == "dag_graph":     
        """
        Function: generate DAG data (as a special case of chain graph models)
        """
        simu, G = data.load_dag_data()
        
        if args.plot:
            U = np.zeros((n, n))
            nx.draw_circular(G,arrowsize=15, with_labels=True)
            plt.show()    
            
        G = nx.to_numpy_array(G)   
        cc = data.connected_comp(n, G)

        membership, c_size, num = data.is_chain_graph(cc)

    
    if args.graph == "npvar_dag":
        """
        Function: Load NPVAR data for algorithm performance test.
        """

        with open('utils.R', 'r') as f:
            string = f.read()
        chain = STAP(string, "chain")
        data_list = chain.npvar_dag_data(n, d)
        """Return simulation DAG data and graph""" 
        simu = data_list[0]
        G = nx.DiGraph(np.array(data_list[1]))  
        npvar_result = data_list[2] -1

        U = np.random.uniform(low=0.5, high=1, size=[n, n])
        
        cc = data.connected_comp(n, nx.to_numpy_array(G)*U)
        membership, c_size, num = data.is_chain_graph(cc)
        
        nx.draw_circular(G, with_labels=True, arrows=True)
        plt.show()   
 
                
    return simu, G, membership, c_size, num, cc, top_sort

def test_chain_graph(G, membership):
    
    n, d = args.n, args.d
    det, c, a = args.det, args.c, args.a
    
    def get_adj_cg(G):

        U = np.random.uniform(low=0.5, high=1.5, size=[n, n])
        U[np.random.rand(n, n) < 0.5] *= -1
        
        G2 = (G != 0).astype(float) * U
        
        return G, G2 
    
    """
    G/G_perm: binary adjacency matrix
    G2/G1: Weighted adjacency matrix
    X: Initialize data 
    """
    
    def equal_det(n):
        import scipy.linalg as la
        """
        Notes: Generate size n eigen vector. The product of n elements equals to 'det';
            To achieve this, we first generate 'n-1' random variables (runif) from a uniform distribution;
            The last element equals to 'det/(prod(runif))';
            Eigenvector is formed by the above n elements;
    
        """
        determinate = args.det  # Determinant of the covariant matrix of a chain component
        if n > 1:
            eigen_val = np.random.uniform(low=float(0.5 / n), high=float((1.5) / n), size=n-1)
            last_val = determinate / np.prod(eigen_val)
            eigen_val = np.append(eigen_val, last_val)
        if n == 1:
            eigen_val = np.expand_dims(determinate, axis=0)
        eigen_vec = np.diag(eigen_val)
        """
        Notes: eigen_vec = diag(λ1, ..., λn) is a diagonal matrix. 
               Q is a nxn orthogonal matrix (Compute QR decomposition of a matrix.)
        """
        s = np.random.choice([-1, 1], size=(n, n), replace=True)
        q, _ = la.qr(np.random.rand(n, n) * s)
        semidef = q.T.dot(eigen_vec).dot(q)
    
        return semidef


    def topological_sort(membership, G2):
        G3 = G2.copy()
        top_sort = []
        
        for i in membership:
            if len(i) == 1:
                if G3[:,i].sum() == 0:
                    top_sort.append(i)
                    G3[i,:] = 0
            else:
                if len(i) > 1:
                    j = np.setdiff1d( range(n), i)
                    if G3[j , :][:,i].sum() == 0:
                        top_sort.append(i)
                        G3[i,:] = 0
            
        return top_sort
    
    def predecessors(G, cc):
        parents = []
        for nodes in cc:
            parent = np.setdiff1d((np.nonzero(G[:,nodes])), cc)
            parents.append(parent)
        return parents

    G, G2 = get_adj_cg(G)   
    X = np.zeros((d, n))
    
    top_sort = topological_sort(membership, G2)
    
    for cc in top_sort:
        """
        Z: undirected chain components (fully connected) generated from inverse precision matrix 
        """
        mu = [0] * len(cc)
        sigma = equal_det(len(cc))
        inv_sigma = np.linalg.inv(sigma) 

        Z = np.random.multivariate_normal(mean=mu, cov=inv_sigma, size=d)

        for nodes in cc:
            parent = np.setdiff1d((np.nonzero(G[:,nodes])), cc)
            if len(parent) == 0:
                X[:, cc] = Z        
            else: 
                M_pa = np.array((X[:, parent].dot(G2[parent, nodes])))
                X[:, cc] = np.expand_dims(M_pa, axis=1) + Z

          
    return X, G, top_sort

                   
"""
    Function: Our Non-parametric equal determiant of covariance matrix algorithm
"""

def NPCOV(simu, dag, cc):

    def est_cond_cov(data, descendants, ancestors):

        """
        Function: Estimate conditional covariance between ancestors and descendants

        Parameters
        ----------
        data : Input simulation data
        descendants : Generate a list of node index contains all descendants
        node_idx : Index of all nodes (eg. n=3, node_idx = [1, 2, 3])
        ancestors : Index of all ancestor nodes
    
        Returns
        -------
        cond_cov : Conditional covariance. COV(descendants | ancestors)
    
        """
        
        descendants = np.expand_dims(descendants, axis=0)
    
        
        for i in range(len(descendants)):
            current_node = descendants[i]
           
            with open('utils.R', 'r') as f:
                string = f.read()
            chain = STAP(string, "chain")    
    
            if args.regress_method == "mgcv":
                est = chain.calculate_mgcv(data, ancestors, current_node)
            if args.regress_method =="np":
                est = chain.calculate_np(data, ancestors, current_node)
    
        return np.array(est)
    
    def find_element_in_list(element, list_element):
        try:
            index_element = list_element.index(element)
            return index_element
        except ValueError:
            return None
        
    def Convert(lst):
        t, dic_lst = 0,  {}
        for it in lst:
            dic_lst[t] = it
            t +=1
        return dic_lst
    
    n = simu.shape[0]
    p = simu.shape[1]
    node_index = list(range(p))
    
    """
    Step1: Find the first source by minimum marginal covariance
    """
    membership, c_size, num = data.is_chain_graph(cc)
    dic_membership = Convert(membership)   
    
    # Notes: Get all chain components (Chain components drops into cc_section)
    cc_section = np.cumsum(num)
    cc_section = np.append(0, cc_section)
    c = len(cc_section) - 1
    
    """Get the determinant of all chain components """
    det_cg = {}  # determinant of chain graph
    
    for idx, item in enumerate(membership):
        
        data_temp = []

        for col in (item):
            data_temp.append(simu[:, col])

        data_temp = np.array(data_temp)
        data_temp = np.expand_dims(data_temp, axis=0)
        if data_temp.shape[0] == 1:
            det = np.var(data_temp)
        else:
            det = np.linalg.det(np.cov(data_temp))
        det_cg[idx] = det
        
    """print('Determinant all', det_cg)"""

    
    """ Store ancestors (ancesors[]) and layer-wise ancestors (layer_ancestors[])"""
    ancestors = []         # Starts with the root node
    layer_ancestors = {}   # Starts from layer 0
    
    """ Find the source with the minimum covariance"""
    min_cov  = min(det_cg.values())
    
    for key, val in det_cg.items():
        val -= min_cov
        if val < args.eta:
            det_cg[key] = 0   

    source_idx = list(det_cg.keys())[list(det_cg.values()).index(0)]
    source_idx = membership[source_idx]
    
    """ Store the initial results """
    layer_ancestors[0] = source_idx
    ancestors.append(source_idx)

    
    """ Then remove the source ancestor node from the node list"""
    
    for key, val in dic_membership.items():  
        if val == source_idx:
            del det_cg[key]
    
    """ Step 2: Main loop of the algorithm 
            a. Estimate conditional variance for each descendant 
            b. Choose node with minimum conditional variance and set it as source"""
                 
    """ Estimate conditional covariance"""

    l = 1   
    while len(ancestors) < (len(num) - 1):
        Cond_cov = []
        descendants_new = []    
        
        """ Update descendants """
        for key, val in det_cg.items():
            descendants_new.append(membership[key])
    
        for i in range(len(descendants_new)):
            res_matrix = []
            
            for descent_temp in descendants_new[i]:
                """ Calculate conditional covariance """
                cond_covariance = est_cond_cov(simu, descent_temp, ancestors)
                res_matrix.append(cond_covariance)

            if np.array(res_matrix).shape[0] == 1 and (simu[:, descendants_new[i]].T.shape[0]) == 1:
                cond_cov = np.var(simu[:, descendants_new[i]].T) - np.var(np.array(res_matrix))
                
            if np.array(res_matrix).shape[0] == 1 and (simu[:, descendants_new[i]].T.shape[0]) > 1:
                cond_cov = np.linalg.det( np.cov(simu[:, descendants_new[i]].T)) - np.var(np.array(res_matrix))
                
            if np.array(res_matrix).shape[0] > 1 and (simu[:, descendants_new[i]].T.shape[0]) > 1:
                cond_cov = np.linalg.det( np.cov(simu[:, descendants_new[i]].T) -  np.cov( np.array(res_matrix)) )  

            """ Calculate conditional covariance for all descendants """
            Cond_cov.append(cond_cov)
                
        """ Find the minimum conditional covariance"""
        
        min_cond_cov  = min(Cond_cov)
        min_cond_cov_idx = descendants_new[Cond_cov.index(min(Cond_cov))]
        
        ancestors.append(min_cond_cov_idx )
        layer_ancestors[l] = min_cond_cov_idx 
        
        for key, val in dic_membership.items():  
            if val == min_cond_cov_idx:
                del det_cg[key]
        l += 1  # Layer plus one
         

    import itertools
    from itertools import chain
    
    ancestors_unzip = list(chain.from_iterable(ancestors))    
    descendant_final = [x for x in node_index if x not in ancestors_unzip]

    ancestors = list(itertools.chain(ancestors, descendant_final))
    layer_ancestors[l+1] = descendant_final

    """
    Lastly, infer adjacency matrix by significance given by gam
    """
    return layer_ancestors, ancestors

"""
    Function: Convert topological order into adjacency matrix
"""
def prune(simu, layer_ancestors):
    
    top_sort = []
    for i in list(layer_ancestors.values()):
        for j in i:
            top_sort.append(j )
    print('topological sort', top_sort )
    
    with open('utils.R', 'r') as f:
        string = f.read()
    chain = STAP(string, "chain")    
    
    adj = chain.prune(simu , top_sort)

    return adj


def top_shd(layer_ancestors, G):
    """
    Function: calculate SHD from topological order
    """    
    # Pre-processing G
    n = G.shape[0]
    diag = np.diag(np.zeros(n) + 0.1)
    G_pre = G + diag
    
    from numpy.linalg import matrix_power
    path = matrix_power(G, n)

    # RE-order ancestors    
    top_sort = []
    for i in list(layer_ancestors.values()):
        top_sort.append(i)

    order=np.array(list(range(len(top_sort)-1, -1, -1)))
    top_reorder = [top_sort[i] for i in order]
    
    top_temp = top_reorder.copy()
    num = len(top_sort) - 1 
 
    SHD = 0
    for cc in top_reorder:
        if cc != top_reorder[num]:
            i = cc[0]
            del top_temp[0]
            for cc_next in top_temp:
                j = cc_next[0]
                if path[i][j] > 0:
                    SHD += 1

    return SHD

  
if __name__ == "__main__":
    
    SHD_cov = []
    
    for i in range(20):
        simu, G, membership, c_size, num, cc, top_sort = load_graph_data()    
           
        layer_ancestors, _ = NPCOV(simu, G, cc)
        shd_cov = top_shd(layer_ancestors, G)
        SHD_cov.append(shd_cov)

    print('SHD_cov = ',   SHD_cov)

