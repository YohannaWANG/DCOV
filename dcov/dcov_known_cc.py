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

def DCOV_known(simu, dag, cc):
    '''
    Implementation of Algorithm 1:  learning the topological order 
    of a chain graph with known chain component decomposition. 
    
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
            data        : Input simulation data
            descendants : Generate a list of node index contains all descendants
            ancestors   : Index of all ancestor nodes
    
        Returns
            cond_cov    : Conditional covariance. COV(descendants | ancestors)
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
    
    def _find_element_in_list(element, list_element):
        
        try:
            index_element = list_element.index(element)
            return index_element
        except ValueError:
            return None
        
    def _convert(lst):
        '''
        Convert list into dictionary
        '''
        t, dic_lst = 0,  {}
        for it in lst:
            dic_lst[t] = it
            t +=1
        return dic_lst

    '''Find the first source chain component by minimumize marginal covariance.'''
    node_index = list(range(simu.shape[1]))

    dic_cc = _convert(cc)   
    
    '''Get the determinant of all chain components'''
    det_cg = {} 
    
    for idx, item in enumerate(cc):
        
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
    
    '''Store ancestors (ancesors[]) and layer-wise ancestors (layer_ancestors[])'''
    ancestors = []         # Starts with the root node
    layer_ancestors = {}   # Starts from layer 0
    
    '''Find the source cc with the minimum covariance'''
    min_cov  = min(det_cg.values())
    
    for key, val in det_cg.items():
        val -= min_cov
        if val < args.eta:
            det_cg[key] = 0   

    source_idx = list(det_cg.keys())[list(det_cg.values()).index(0)]
    source_idx = cc[source_idx]
    
    '''Store the initial result'''
    layer_ancestors[0] = source_idx
    ancestors.append(source_idx)
    
    '''Then remove the source ancestor node from the node list'''
    for key, val in dic_cc.items():  
        if val == source_idx:
            del det_cg[key]
    
    '''
    Main loop of the algorithm 
        a. estimate conditional variance for each descendant 
        b. select node with minimum conditional variance and set it as the source cc
    '''
                 
    '''Estimate conditional covariance'''

    l = 1   # Initialize layers (layer-wise recovery)
    num = []
    for i in range(len(cc)):
        num.append(len(cc[i]))
    while len(ancestors) < (len(num) - 1):
        Cond_cov = []
        descendants_new = []    
        
        ''' Update descendants in each iteration'''
        for key, val in det_cg.items():
            descendants_new.append(cc[key])
    
        for i in range(len(descendants_new)):
            res_matrix = []
            
            for descent_temp in descendants_new[i]:
                '''Calculate conditional covariance'''
                cond_covariance = _est_cond_cov(simu, descent_temp, ancestors)
                res_matrix.append(cond_covariance)

            if np.array(res_matrix).shape[0] == 1 and (simu[:, descendants_new[i]].T.shape[0]) == 1:
                cond_cov = np.var(simu[:, descendants_new[i]].T) - np.var(np.array(res_matrix))

            if np.array(res_matrix).shape[0] == 1 and (simu[:, descendants_new[i]].T.shape[0]) > 1:
                cond_cov = np.linalg.det( np.cov(simu[:, descendants_new[i]].T)) - np.var(np.array(res_matrix))
                
            if np.array(res_matrix).shape[0] > 1 and (simu[:, descendants_new[i]].T.shape[0]) > 1:
                cond_cov = np.linalg.det( np.cov(simu[:, descendants_new[i]].T) -  np.cov( np.array(res_matrix)) ) 

            '''Calculate conditional covariaNPnce for all descendants '''
            Cond_cov.append(cond_cov)
                
        ''' Find the minimum conditional covariance '''
        min_cond_cov_idx = descendants_new[Cond_cov.index(min(Cond_cov))]
        
        ancestors.append(min_cond_cov_idx )
        layer_ancestors[l] = min_cond_cov_idx 
        
        for key, val in dic_cc.items():  
            if val == min_cond_cov_idx:
                del det_cg[key]
        l += 1  # Layer plus one
         

    import itertools
    from itertools import chain
    
    ancestors_unzip = list(chain.from_iterable(ancestors))    
    descendant_final = [x for x in node_index if x not in ancestors_unzip]

    ancestors = list(itertools.chain(ancestors, descendant_final))
    layer_ancestors[l+1] = descendant_final
    
    print('layer ancestors', layer_ancestors)

    ''' Infer adjacency matrix by significance '''
    top_sort = []
    for i in list(layer_ancestors.values()):
        top_sort.append(i)

    adj = chain_graph.prune(simu, top_sort)
    
    return adj 

    
    
    
    