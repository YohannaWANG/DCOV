"""
@author: yohanna
@email: yohanna.wang0924@gmail.com
"""
import numpy as np
from utils import args

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages

from rpy2.robjects.packages import importr
from rpy2.robjects.packages import STAP

import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri

pandas2ri.activate()
rpy2.robjects.numpy2ri.activate()


packageNames = ('graph', 'ggm', 'mgcv', 'pcalg')
utils = rpackages.importr('utils')

#utils.chooseCRANmirror(ind=1)
#packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]

importr('mgcv')
importr('pcalg')
importr('lcd')
robjects.r.source("utils.R")
    
         
with open('utils.R', 'r') as f:
    string = f.read()
chain_graph = STAP(string, "chain_graph")     


import matlab
import matlab.engine
eng = matlab.engine.start_matlab()

from sfo_python import Our_Experiment
def DCOV_unknown(simu, dag, cc):
    '''
    This is the implementation of algorithm 2
            1: We use submodular function optimization to extract chain components from nodes;
            2: Then we recover the topological order use the method in algorithm 1;
            3. Next, we use 'prune' function to recover adjacency matrix from topological order
    
    Import Matlab function: Submodular function optimization;        
    
    Arguments:
        simu    : chain graph data (known_cc)
        dag     : chain graph structure
        cc      : chain components
        
    Returns:
        adj     : estimated adjacency matrix
    '''

    def _est_cond_cov(data, descendants, ancestors):
        '''
        Estimate conditional covariance between ancestors and descendants
    
        Arguments:
            data : Input simulation data
            descendants : Generate a list of node index contains all descendants
            ancestors : Index of all ancestor nodes
    
        Returns:
            cond_cov : Conditional covariance. Eg, COV(descendants | ancestors)
        '''
        
        descendants = np.expand_dims(descendants, axis=0)
        
        for i in range(len(descendants)):
            current_node = descendants[i]
            '''
            mgcv : Nonparametric regression using GAM / mgcv package
            np   : Nonparametric regression using Kernel Methods / np package
            '''
            if args.regress_method == "mgcv":
                est = chain_graph.calculate_mgcv(data, ancestors, current_node)
            if args.regress_method =="np":
                est = chain_graph.calculate_np(data, ancestors, current_node)
    
        return np.array(est)

    def _top_to_adj(simu, top_sort):
        '''  
        Get undirected connection
            Notes: note that undirected connections of chain components can be easily 
            achieved by using GLASSO;
                eg. X_tau = A*X_pa_tau  + Z
                where Z = X_tau - a*X_pa_tau;
                The undirected connection is inv(Z)
        '''
        n = simu.shape[1]
        adj = np.zeros((n, n))
        
        for val in top_sort.values():
            val = np.squeeze(val)
            adj[int(val[0]), int(val[1])] = adj[int(val[1]), int(val[0])] = 1

        '''
        Get directed connection:
            The directed egdes is discovered by performing prune between chain components;
        '''
        top_order = []
        for i in top_sort.values():
            top_order.append(i)
        top_order = np.array(top_order).flatten()   

        test = chain_graph.prune_unknown_cc(simu, top_order, len(top_order))
        adj += test
        adj = np.where(adj > 0, 1.0, 0.0)

        return adj


    node_index = list(range(simu.shape[1]))

    ''''Find the first root chain components by minimizing the SFO function'''
    Sigma = np.cov(simu.T)
    
    if args.sfo == 'sfo-python':
        cc = Our_Experiment.our_experiment(Sigma)
        cc = np.expand_dims(cc, axis=0)
        
    elif args.sfo == 'sfo-matlab':
        Sigma = matlab.double(Sigma.tolist())
        cc = np.array(eng.sfo_min_cg(Sigma))
    
    '''Re-arrange the data format: Matleb index starts from 1 while python starts from 0'''
    cc_new = []
    for i in cc: cc_new.append(i-1)
    cc_new = np.sort(cc_new)
    cc_new = np.squeeze(cc_new, axis=0)

    '''Store the initial ancestor. Set all other nodes to be descendant'''
    top_order = {}   
    top_order[0] = cc_new
    ancestors = np.squeeze(cc_new).astype(int)
    descendant = np.array(np.setdiff1d(node_index , np.array(cc_new)))    

    '''Calculate the conditional covariance'''

    l = 1   # Initialize layers (layer-wise recovery)
    while len(ancestors) < (len(node_index) - 1):
        descendant_new = []    
        res_matrix = []

        for descent_temp in descendant:
            cond_covariance = _est_cond_cov(simu, descent_temp, ancestors)
            res_matrix.append(cond_covariance)    
        
        if np.array(res_matrix).shape[0] == 1 and (simu[:, descendant].T.shape[0]) == 1:
            cond_cov = np.var(simu[:, descendant].T) - np.var(np.array(res_matrix))

        if np.array(res_matrix).shape[0] == 1 and (simu[:, descendant].T.shape[0]) > 1:
            cond_cov = np.linalg.det( np.cov(simu[:, descendant].T)) - np.var(np.array(res_matrix))
            
        if np.array(res_matrix).shape[0] > 1 and (simu[:, descendant].T.shape[0]) > 1:
            cond_cov = np.cov(simu[:, descendant].T) -  np.cov( np.array(res_matrix))  

        cond_cov = np.array(cond_cov)
        Sigma = matlab.double(cond_cov.tolist())
        
        '''
        Note: This will retur new node index from Matlab;
              need to rearrange the index for the descendant nodes
        '''
        if args.sfo == 'sfo-matlab':
            cc_temp = np.array(eng.sfo_min_cg(Sigma))
        
        elif args.sfo == 'sfo-python':
            cc_temp = Our_Experiment.our_experiment(Sigma)
            cc_temp = np.expand_dims(cc_temp, axis=0)
            
        cc_index = []
        for i in cc_temp: cc_index.append(i-1)
        cc_index = np.squeeze(cc_index, axis = 0)

        for i in cc_index:
            i = np.int(i)
            descendant_new = np.append(descendant_new, descendant[i])
            descendant_new = np.sort(descendant_new)
        
        top_order[l] = descendant_new
        
        ancestors = np.append(ancestors, descendant_new)
        descendant = np.array(np.setdiff1d(node_index , np.array(ancestors)))
        l += 1

    top_sort = []
    for i in list(top_order.values()):
        top_sort.append(i)
    
    print('top_sort', top_sort)
    adj = _top_to_adj(simu, top_order)
    
    return adj    
