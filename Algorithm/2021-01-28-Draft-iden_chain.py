"""
    Title:  Identifiability of chain graph models;
    Author: Yuhao Wang
    Email:  yohanna.wang0924@gmail.com
    Date:   2021-02-02
"""

import argparse
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import STAP

import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri

#import os
#os.environ['R_HOME'] = '/home/yohanna/R/x86_64-pc-linux-gnu-library/4.0'


pandas2ri.activate()
rpy2.robjects.numpy2ri.activate()


packageNames = ('graph', 'RBGL', 'ggm', 'mgcv', 'pcalg')
utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

#packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]

"""
    Notes: Uncomment this part if you don't have the following packages installed on ur comp.
"""
# Running R in Python example installing packages:
#if len(packnames_to_install) > 0:
#    utils.install_packages(StrVector(packnames_to_install))

importr('mgcv')
importr('pcalg')
importr('lcd')
robjects.r.source("utils.R")

parser = argparse.ArgumentParser()
# ----------- data parameters --------


parser.add_argument('--det', type=int, default=1,
                    help='equal determinant for chain components')
parser.add_argument('--n', type=int, default=10, help='number of nodes')
parser.add_argument('--c', type=int, default=3, help='number of chain components')
parser.add_argument('--a', type=int, default=1, help='average node degree')
parser.add_argument('--d', type=int, default=1000, help="number of dimensions")

parser.add_argument('--operator', type=str, default='det', 
help="super-additive function on positive semidefinite matrices. 1.det, 2.trace")
parser.add_argument('--eta', type=float, default=0.01, help="residual error variance")

args = parser.parse_args()

"""
Function: generated chain graph adjacency matrix.
"""


def gen_cg_adj(n, c, a, d):
    """
    Args:
        n: number of nodes
        c: number of chain components
        a: average node degree
        d: number of data dimensions

    Returns: A, A_undirect, D, U, sorted_D, sorted_U
    """
    """Generate chain graph"""
    V = list(range(n))
    A = np.zeros((n, n))

    """Fill lower triangle"""
    p = float(d) / (n - 1)
    for j in range(n):
        for i in range(j + 1, n):
            r = np.random.uniform(0, 1)
            if r <= p: A[i][j] = 1  # problem is we have no prior knowledge of undirected edges
    A = A + A.T

    """Assign chain components"""
    # (4) Select an integer c randomly as from {1,...,n} as the number of chain components;
    # Seperate nindex np.array(len(0, n)) into R sections, forms R chain components;
    R = np.random.choice(n - 1, c, replace=False) + 1
    R = np.append(np.append(0, R), n)
    R.sort()

    C = {}
    V = list(range(n))
    for i in range(len(R) - 1):
        for v in V[R[i]: R[i + 1]]:
            C[v] = i

    # (6) Set A_{ij}=0 for any (i,j) pair such that i\in C(l), j\in C(m) (l>m);
    D = []  # Directed edges
    U = []  # Undirected edges
    for i in range(n):
        for j in range(n):
            if C[i] > C[j]:
                if A[i][j] == 1:
                    D.append((j, i))
                A[i][j] = 0
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1 and A[j][i] == 1:
                U.append((i, j))

    # Find nodes with directed connections "sorted_D[]" and undirected connections "sorted_U[]"
    sorted_D = []  # Nodes involved in directed connections
    sorted_U = []  # Nodes involved in undirected connections
    for i in range(len(D)):
        sort_D = np.sort(np.array(D[i]))
        sorted_D.append(sort_D)
    sorted_D = np.unique(sorted_D)

    for i in range(len(U)):
        sort_U = np.sort(np.array(U[i]))
        sorted_U.append(sort_U)
    sorted_U = np.unique(sorted_U, axis=0)

    """
    Function: Determine chain components (size, number)
    """
    A_undirect = A.copy()  # A_undirect: undirected edges from chain graph A
    # Remove directed connections, keep undirected connections and determine clusters
    for i in range(len(D)):
        j, k = D[i]
        A_undirect[(j, k)] = 0

    temp_u_s = []  # Undirected clusters with two nodes only;
    temp_u_m = []  # Undirected clusters with more than two nodes.
    for i in range(A_undirect.shape[0]):
        for j in range(i + 1, n):
            if np.sum(A_undirect[i]) == 1 and A_undirect[i, j] > 0:
                temp_u_s.append((i, j))
            else:
                if np.sum(A_undirect[i]) > 1:
                    if A_undirect[i, j] - 1 == 0:
                        temp_u_m.append((i, j))
    temp_u_m = np.unique(temp_u_m)

    # cc: connected chain components within each undirected connections.
    if len(temp_u_m) > 0 and len(temp_u_s) > 0:
        cc = temp_u_s.append(temp_u_m)
    else:
        if len(temp_u_m) > 0 and len(temp_u_s) == 0:
            cc = temp_u_m
        if len(temp_u_m) == 0 and len(temp_u_s) > 0:
            cc = temp_u_s
    return A, A_undirect, D, U, sorted_D, sorted_U


"""
Function: Determine Connected components in an undirected graph
          Determine if graph contains cycles or not
"""


class Graph:
    # init function to declare class variables
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]

    def DFSUtil(self, temp, v, visited):

        # Mark the current vertex as visited
        visited[v] = True

        # Store the vertex to list
        temp.append(v)

        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v]:
            if visited[i] == False:
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

    # method to add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        # self.adj[w].append(v)

    # Method to retrieve connected components
    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc

    # Check is the input graph contains cycles
    def isCyclicUtil(self, v, visited, recStack):

        # Mark current node as visited and
        # adds to recursion stack
        visited[v] = True
        recStack[v] = True

        # Recur for all neighbours
        # if any neighbour is visited and in
        # recStack then graph is cyclic
        for neighbour in self.adj[v]:
            if visited[neighbour] == False:
                if self.isCyclicUtil(neighbour, visited, recStack) == True:
                    return True
            elif recStack[neighbour] == True:
                return True

        # The node needs to be poped from
        # recursion stack before function ends
        recStack[v] = False
        return False

    # Returns true if graph is cyclic else false
    def isCyclic(self):
        visited = [False] * self.V
        recStack = [False] * self.V
        for node in range(self.V):
            if visited[node] == False:
                if self.isCyclicUtil(node, visited, recStack) == True:
                    return True
        return False


"""
Function: use undiretced graph as input, return connected components.
"""


def connected_comp(n, A_undirect):
    """
    Args:
        n: node number
        A_undirect: undirected adjancency matrix from chain graph A.
    Returns:
        cc: connected chain components
    """
    g = Graph(n)
    for v in range(A_undirect.shape[0]):
        for w in range(v + 1, A_undirect.shape[1]):
            if A_undirect[v, w] == 1:
                g.addEdge(v, w)
    cc = g.connectedComponents()
    return cc


def is_chain_graph(cc):
    """
    Args:
        cc: connected chain components

    Returns:
        membership: nodes in each chain components;
        csize     : size of each chain components;
        num       : number of chain components.
    """
    membership = cc
    c_size = len(cc)
    num = []
    for i in range(len(cc)):
        num.append(len(cc[i]))
    return membership, c_size, num


def is_acyclic(A):
    """
    Args:
        A: chain graph adjacency matrix
    Returns:
        Bool: True/ False- Whether graph contains cycles

    """
    g = Graph(n)
    for v in range(A.shape[0]):
        for w in range(v + 1, A.shape[1]):
            if A[v, w] == 1:
                g.addEdge(v, w)
    if g.isCyclic() == 1:
        print("Graph has a cycle")
        return True
    else:
        print("Graph has no cycle")
        return False


"""
Function: plot chain graph ground truth
"""


def plot_chain_graph(cc, D, U):
    """
    Args:
        cc: clustered chain components
        D:  directed edges in chain graph
        U:  undirected edges in chain graph

    Returns:

    """
    C = {}
    with open("colours.txt") as f:
        lines = f.readlines()

    for i in range(len(cc)):
        r, g, b = lines[i].strip().split(",")
        r, g, b = int(r), int(g), int(b)
        C[i] = (r, g, b)

    # Convert to hex
    def clamp(x):
        return max(0, min(x, 255))

    for i in C.keys():
        r, g, b = C[i]
        h = "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))
        C[i] = h

    # Add network vertices
    V = list(range(n))  # Nodes indices

    # Assign colors to each nodes in different chain components
    dic_cc = {}
    for i in range(len(cc)):
        for j in cc[i]:
            dic_cc[j] = i
    nc = [C[dic_cc[v]] for v in V]

    # Add networkx edges
    G = nx.DiGraph(directed=True)  # Initialize a graph
    for v in V: G.add_node(v)
    for (v, w) in U: G.add_edge(v, w, w=6, c='#bbbbbb')
    for (v, w) in D: G.add_edge(v, w, w=2, c='#000000')

    width = [G[u][v]['w'] for u, v in G.edges()]
    colours = [G[u][v]['c'] for u, v in G.edges()]

    pos = nx.nx_pydot.graphviz_layout(G)
    nx.draw(G, pos=pos, with_labels=True, font_size=10)

    nx.draw_networkx_nodes(G, pos=pos, node_color=nc, node_size=500)
    nx.draw_networkx_edges(G, pos=pos, width=width, edge_color=colours)
    plt.savefig('chain_graph.png', dpi=1000)
    plt.show()

    
"""
    Function: Generate chain graph adjacency matrix 
"""    

def gen_adj_cg(A, num, d):
    
    def symmetric(A, tol=1e-8):
        return np.all(np.abs(A - A.T) < tol)

    def PSD(A):
        M = np.matrix(A)
        return symmetric(M) and np.all(np.linalg.eigvals(M) >= 0)
    # Generate adjacency matrix B
    cc_section = np.cumsum(num)
    cc_section = np.append(0, cc_section)
    n, c = cc_section[-1], len(cc_section) - 1
    
    f, i = False, 0
    Ac, B = A.copy(), np.zeros((n, n))

    while not f:
        f = True
        i = i + 1
        for j in range(c):
            aa = n - cc_section[j]
            bb = cc_section[j + 1] - cc_section[j]

            runif = np.random.uniform(low=float(0.5 / bb), high=float((1.5) / bb), size=aa * aa)
            s = np.random.choice([-1, 1], size=aa * aa, replace=True)
            randdist = runif * s
            
            randdist = np.reshape(randdist, (randdist.size))
            randdist.shape = (randdist.size // aa, aa)
            randdist = np.maximum(randdist, randdist.transpose())
            np.fill_diagonal(randdist, 1)

            B[cc_section[j]:n, cc_section[j]:n] = randdist
            if j != c - 1:
                B[cc_section[j + 1]:n, cc_section[j]:cc_section[j + 1]] = 0

        Ac[range(n), range(n)] = 1
        B = B * (Ac.T)
        for j in range(c):
            I = range(cc_section[j], cc_section[j + 1])
            if len(I) > 1:
                M = B[cc_section[j]:cc_section[j + 1], cc_section[j]:cc_section[j + 1]]
                if not PSD(M):
                    f = False
                    break
        return B


"""
Function: generate simulation data (rnorm.cg) without equal determinant
"""


def gen_rnorm_data(A, num, d):
    """
    Parameters
    ----------
    A : Chain graph adjacency matrix
    num : Number of chain components
    d : number of observational data (X or Y)

    Returns
    -------
    X: Simulation data

    """
    def symmetric(A, tol=1e-8):
        return np.all(np.abs(A - A.T) < tol)

    def PSD(A):
        M = np.matrix(A)
        return symmetric(M) and np.all(np.linalg.eigvals(M) >= 0)

    # (7) Given X via the {rnorm.cg} function (LCD R package) from A.
    cc_section = np.cumsum(num)
    cc_section = np.append(0, cc_section)
    n, c = cc_section[-1], len(cc_section) - 1
    
    #f, i = False, 0
    #Ac, B = A.copy(), np.zeros((n, n))
    B = gen_adj_cg(A, num, d)
    '''get normal dist: https://cran.microsoft.com/snapshot/2015-04-22/web/packages/lcd/lcd.pdf
    while not f:
        f = True
        i = i + 1
        for j in range(c):
            aa = n - cc_section[j]
            bb = cc_section[j + 1] - cc_section[j]
            """
            TODO: generate data with equal determinant
            """
            runif = np.random.uniform(low=float(0.5 / bb), high=float((1.5) / bb), size=aa * aa)
            s = np.random.choice([-1, 1], size=aa * aa, replace=True)
            randdist = runif * s
            
            randdist = np.reshape(randdist, (randdist.size))
            randdist.shape = (randdist.size // aa, aa)
            randdist = np.maximum(randdist, randdist.transpose())
            np.fill_diagonal(randdist, 1)

            B[cc_section[j]:n, cc_section[j]:n] = randdist
            if j != c - 1:
                B[cc_section[j + 1]:n, cc_section[j]:cc_section[j + 1]] = 0

        Ac[range(n), range(n)] = 1
        B = B * (Ac.T)
        for j in range(c):
            I = range(cc_section[j], cc_section[j + 1])
            if len(I) > 1:
                M = B[cc_section[j]:cc_section[j + 1], cc_section[j]:cc_section[j + 1]]
                if not PSD(M):
                    f = False
                    break 
        '''
    X = np.zeros((n, d))
    for j in range(c):
        #I = range(cc_section[j], cc_section[j + 1])
        mu, sigma = [0] * (cc_section[j + 1] - cc_section[j]), B[cc_section[j]:cc_section[j + 1], cc_section[j]:cc_section[j + 1]]
        #print('determinant_rnorm = ', np.linalg.det(sigma))
        X[cc_section[j]:cc_section[j + 1]] = np.random.multivariate_normal(mean=mu, cov=sigma, size=d).T

    return np.linalg.solve(B, X).T
       
"""
    Function: Generate simulation data with equal determinant
"""    

def gen_data(A, num, d):
    
    '''
    Args:
        A: A represent the chain graph adjacency matrix;
        cc: list of end indices of chain components;
        d: number of observations
    Returns:
        X = AX+ N
        X: Represent the observational samples;
        N: Independent noise
    '''
    
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

    cc_section = np.cumsum(num)
    cc_section = np.append(0, cc_section)
    n, c = cc_section[-1], len(cc_section) - 1
    
    B = gen_adj_cg(A, num, d)
        
    X = np.zeros((n, d))
    for j in range(c):
        mu = [0] * (cc_section[j + 1] - cc_section[j])
        sigma = equal_det(cc_section[j + 1] - cc_section[j])
        #print('determinant_eql = ', np.linalg.det(sigma))
        X[cc_section[j]:cc_section[j + 1]] = np.random.multivariate_normal(mean=mu, cov=sigma, size=d).T

        return np.linalg.solve(B, X).T

    
"""
    Function: Performance evaluation matrics.
"""
 
def count_accuracy(G_true: nx.DiGraph,
                   G: nx.DiGraph,
                   G_und: nx.DiGraph = None) -> tuple:
    
    """Compute FDR, TPR, and FPR for B, or optionally for CPDAG B + B_und.
    Args:
        G_true: ground truth graph
        G: predicted graph
        G_und: predicted undirected edges in CPDAG, asymmetric
    Returns:
        fdr: (reverse + false positive) / prediction positive
        tpr: (true positive) / condition positive
        fpr: (reverse + false positive) / condition negative
        shd: undirected extra + undirected missing + reverse
        nnz: prediction positive
    """
    B_true = nx.to_numpy_array(G_true) != 0
    B = nx.to_numpy_array(G) != 0
    B_und = None if G_und is None else nx.to_numpy_array(G_und)
    d = B.shape[0]
    # linear index of nonzeros
    if B_und is not None:
        pred_und = np.flatnonzero(B_und)
    pred = np.flatnonzero(B)
    cond = np.flatnonzero(B_true)
    cond_reversed = np.flatnonzero(B_true.T)
    cond_skeleton = np.concatenate([cond, cond_reversed])
    # true pos
    true_pos = np.intersect1d(pred, cond, assume_unique=True)
    
    if B_und is not None:
        # treat undirected edge favorably
        true_pos_und = np.intersect1d(pred_und, cond_skeleton, assume_unique=True)
        true_pos = np.concatenate([true_pos, true_pos_und])
    # false pos
    false_pos = np.setdiff1d(pred, cond_skeleton, assume_unique=True)
    if B_und is not None:
        false_pos_und = np.setdiff1d(pred_und, cond_skeleton, assume_unique=True)
        false_pos = np.concatenate([false_pos, false_pos_und])
    # reverse
    extra = np.setdiff1d(pred, cond, assume_unique=True)
    reverse = np.intersect1d(extra, cond_reversed, assume_unique=True)
    # compute ratio
    pred_size = len(pred)
    if B_und is not None:
        pred_size += len(pred_und)
    cond_neg_size = 0.5 * d * (d - 1) - len(cond)
    fdr = float(len(reverse) + len(false_pos)) / max(pred_size, 1)
    tpr = float(len(true_pos)) / max(len(cond), 1)
    fpr = float(len(reverse) + len(false_pos)) / max(cond_neg_size, 1)
    # structural hamming distance
    B_lower = np.tril(B + B.T)
    if B_und is not None:
        B_lower += np.tril(B_und + B_und.T)
    pred_lower = np.flatnonzero(B_lower)
    cond_lower = np.flatnonzero(np.tril(B_true + B_true.T))
    extra_lower = np.setdiff1d(pred_lower, cond_lower, assume_unique=True)
    missing_lower = np.setdiff1d(cond_lower, pred_lower, assume_unique=True)
    shd = len(extra_lower) + len(missing_lower) + len(reverse)
    
    return fdr, tpr, fpr, shd, pred_size


"""
    Function: Chain graph baseline algorithms 
"""

def baseline_cg(data, graph):

    
    """ Function: use rpy for chain graph structure analysis"""

    with open('utils.R', 'r') as f:
        string = f.read()
    chain = STAP(string, "chain")
    
    """ Baseline 1: LCD-like algorithm """ 
    adj_lcd = chain.baseline_lcdlike(data, graph)
    
    """ Baseline 2: LDCG algorithm"""
    adj_ldcg = chain.baseline_ldcg(data, graph)
    
    """ Baseline 3: PC-like algorithm """
    adj_pclike = chain.baseline_pclike(data, graph)
    
    """
    Notes: Below are baseline algorithms for DAG structure learning
    """
    """ Baseline 4: pc algorithm for chain graph structure learning"""
    adj_pc = chain.baseline_pc(data, graph)

    """ Baseline 5: GES algorithm """
    adj_ges = chain.baseline_ges(data, graph)
    
    """ Baseline 6: RFCI algorithm """
    adj_rfci = chain.baseline_rfci(data, graph)
    
    
    fdr_lcd, tpr_lcd, fpr_lcd, shd_lcd, pred_size_lcd = count_accuracy(nx.DiGraph(graph), nx.DiGraph(adj_lcd))
    print('LCD_like Graph Accuracy: fdr', fdr_lcd, ' tpr ', tpr_lcd, ' fpr ', fpr_lcd, 'shd', shd_lcd, 'pred_size', pred_size_lcd)

    fdr_ldcg, tpr_ldcg, fpr_ldcg, shd_ldcg, pred_size_ldcg = count_accuracy(nx.DiGraph(graph), nx.DiGraph(adj_ldcg))
    print('LDCG Graph Accuracy: fdr', fdr_ldcg, ' tpr ', tpr_ldcg, ' fpr ', fpr_ldcg, 'shd', shd_ldcg, 'pred_size', pred_size_ldcg)

    fdr_pclike, tpr_pclike, fpr_pclike, shd_pclike, pred_size_pclike = count_accuracy(nx.DiGraph(graph), nx.DiGraph(adj_pclike))
    print('PC-LIKE Graph Accuracy: fdr', fdr_pclike, ' tpr ', tpr_pclike, ' fpr ', fpr_pclike, 'shd', shd_pclike, 'pred_size', pred_size_pclike)

    fdr_pc, tpr_pc, fpr_pc, shd_pc, pred_size_pc = count_accuracy(nx.DiGraph(graph), nx.DiGraph(adj_pc))
    print('PC Graph Accuracy: fdr', fdr_pc, ' tpr ', tpr_pc, ' fpr ', fpr_pc, 'shd', shd_pc, 'pred_size', pred_size_pc)
    
    fdr_ges, tpr_ges, fpr_ges, shd_ges, pred_size_ges = count_accuracy(nx.DiGraph(graph), nx.DiGraph(adj_ges))
    print('GES Graph Accuracy: fdr', fdr_ges, ' tpr ', tpr_ges, ' fpr ', fpr_ges, 'shd', shd_ges, 'pred_size', pred_size_ges)

    fdr_rfci, tpr_rfci, fpr_rfci, shd_rfci, pred_size_rfci = count_accuracy(nx.DiGraph(graph), nx.DiGraph(adj_rfci))
    print('RFCI Graph Accuracy: fdr', fdr_rfci, ' tpr ', tpr_rfci, ' fpr ', fpr_rfci, 'shd', shd_rfci, 'pred_size', pred_size_rfci)
    

def baseline_cg_lcdlike(data, graph):
    """ Function: use rpy for chain graph structure analysis
                 Two ways for performance evaluation, R version follows AMPCG paper
        TODO: To be decide later.
    """

    with open('utils.R', 'r') as f:
        string = f.read()
    chain = STAP(string, "chain")
    
    """ Baseline 1: LCD-like algorithm """ 
    lcd = chain.baseline_lcdlike(data, graph)      
    TP_lcd, FN_lcd, FP_lcd, TN_lcd, TPR_lcd, FPR_lcd, ACC_lcd, SHD_lcd = lcd[0], lcd[1], lcd[2], lcd[3], lcd[4], lcd[5], lcd[6], lcd[7]
    print('LCD accuracy from R: TP', TP_lcd, 'FN', FN_lcd, 'FP', FP_lcd, 'TN', TN_lcd, 'TPR', TPR_lcd, 'FPR', FPR_lcd, 'ACC', ACC_lcd, 'SHD', SHD_lcd)
    
    adj_lcd2 = chain.baseline_lcdlike2(data, graph)
    
    fdr_lcd, tpr_lcd, fpr_lcd, shd_lcd, pred_size_lcd = count_accuracy(nx.DiGraph(graph), nx.DiGraph(adj_lcd2))
    print('LCD_like Graph Accuracy: fdr', fdr_lcd, ' tpr ', tpr_lcd, ' fpr ', fpr_lcd, 'shd', shd_lcd, 'pred_size', pred_size_lcd)     


"""
    Function: Estimate conditional covariance between ancestors and descendants
"""
def est_cond_cov(data, descendants, ancestors):
    from pygam import PoissonGAM, s, te
    
    """
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
    from sklearn.feature_extraction import DictVectorizer
    cond_cov = np.zeros(descendants)
    
    """Convert dictionary elements into """
    descendants = np.expand_dims(descendants, axis=0)

    
    for i in range(len(descendants)):
        current_node = descendants[i]
        #ancestors    = ancestors
        #print('current_node', current_node)
        # Compute non-parametric regression using GAM(Generalized Additive Model with integrated smoothness)
        #X = data[:, ancestors]
        #y = data[:, current_node]
        
        with open('utils.R', 'r') as f:
            string = f.read()
        chain = STAP(string, "chain")    
        est = chain.calculate_np(data, ancestors, current_node)

    return np.array(est)
                   
    
"""
    Function: Our Non-parametric equal determiant of covariance matrix algorithm
"""

def NPCOV(data, dag):
    
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
    
    n = data.shape[0]
    p = data.shape[1]
    node_index = list(range(p))
    
    """
    Step1: Find the first source by minimum marginal covariance
    """
    membership, c_size, num = is_chain_graph(cc)
    dic_membership = Convert(membership)   
    print('membership = ', membership)    
    
    # Notes: Chain components drops into cc_section
    cc_section = np.cumsum(num)
    cc_section = np.append(0, cc_section)
    c = len(cc_section) - 1
    
    """Get the determinant of all chain components """
    det_cg = {}  # determinant of chain graph
    
    for idx, item in enumerate(membership):
        
        data_temp = []

        for col in (item):
            data_temp.append(data[:, col])

        data_temp = np.array(data_temp)
        data_temp = np.expand_dims(data_temp, axis=0)
        if data_temp.shape[0] == 1:
            det = np.var(data_temp)
        else:
            det = np.linalg.det(np.cov(data_temp))
        det_cg[idx] = det
        
    #print('Determinant all', det_cg)
    
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
    
    #print('source_idx', source_idx)
    
    """ Store the initial results """
    layer_ancestors[0] = source_idx
    ancestors.append(source_idx)
    
    print('ancestor = ', layer_ancestors)    
    
    """ Then remove the source ancestor node from the node list"""
    #del det_cg[source_idx] 
    for key, val in dic_membership.items():  
        if val == source_idx:
            del det_cg[key]
    
    print('New det_cg is ', det_cg)
    
    """ Step 2: Main loop of the algorithm 
            a. Estimate conditional variance for each descendant 
            b. Choose node with minimum conditional variance and set it as source"""
                 
    """ Estimate conditional covariance"""

    l = 1   # Initialize layers (layer-wise recovery)
    while len(ancestors) < (len(num) - 1):
        res_matrix = []
        Cond_cov = []
        descendants_new = []    
        
        """ Update descendants """
        for key, val in det_cg.items():
            descendants_new.append(membership[key])
    
        for i in range(len(descendants_new)):

            for descent_temp in descendants_new[i]:
                """ Calculate conditional covariance """
                cond_covariance = est_cond_cov(data, descent_temp, ancestors)
                res_matrix.append(cond_covariance)
                
            # Get data within determinant
            """ Q1 : term 2 is too small here """

            if np.array(res_matrix).shape[0] == 1 and (data[:, descendants_new[i]].T.shape[0]) == 1:
                cond_cov = np.var(data[:, descendants_new[i]].T) - np.var(np.array(res_matrix))
                
            if np.array(res_matrix).shape[0] == 1 and (data[:, descendants_new[i]].T.shape[0]) > 1:
                cond_cov = np.linalg.det( np.cov(data[:, descendants_new[i]].T)) - np.var(np.array(res_matrix))
                
            if np.array(res_matrix).shape[0] > 1 and (data[:, descendants_new[i]].T.shape[0]) > 1:
                term2 = np.linalg.det( np.cov( np.array(res_matrix)))
                cond_cov = np.linalg.det( np.cov(data[:, descendants_new[i]].T))  - term2

            """ Calculate conditional covariance for all descendants """
            Cond_cov.append(cond_cov)
                
        #print('Conditional covariance ', Cond_cov)
        
        """ Find the minimum conditional covariance"""
        
        min_cond_cov  = min(Cond_cov)
        min_cond_cov_idx = descendants_new[Cond_cov.index(min(Cond_cov))]
        
        ancestors.append(min_cond_cov_idx )
        layer_ancestors[l] = min_cond_cov_idx 
        
        #dic_membership = Convert(membership)
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
    
    print('layer ancestors', layer_ancestors)
    print('topological sort final = ', ancestors )

    """
    Lastly, infer adjacency matrix by significance given by gam
    """
    



if __name__ == "__main__":
    """
    Function: get parameters for chain graph simulation data generation
    Param:
        det = 3  # Equal determinant
        n = 5  # number of nodes
        c = 2  # number of chain components
        a = 2  # average node degree
        d = 10000  # number of dimensions
    """
    det, n, c, a, d = args.det, args.n, args.c, args.a, args.d
    A, A_undirect, D, U, sorted_D, sorted_U = gen_cg_adj(n, c, a, d)
    cc = connected_comp(n, A_undirect)

    if is_acyclic(A) == 1:
        print("Contains cycles in this graph.")

    membership, c_size, num = is_chain_graph(cc)

    data = gen_data(A, num, d)
    #print('data = ', data)
    
    '''
    data2 = gen_rnorm_data(A, num, d)
    #data2 = data2.T

    plot_chain_graph(cc, D, U)
    
    NPCOV(data2,  A)
    
    """
    Performance evaluation with baseline algorithms
    """
    #baseline_cg(data, A)
    #baseline_cg_lcdlike(data2, A)
    #baseline_cg_lcdlike(data,  A)
    '''
    

def NPCOV_draft(data, dag):
    
    membership, c_size, num = is_chain_graph(cc)
   
    
    # Notes: Chain components drops into cc_section
    cc_section = np.cumsum(num)
    cc_section = np.append(0, cc_section)
    n, c = cc_section[-1], len(cc_section) - 1
    """ Apply NPCOV through 3 implementations
    
        Step 1: Natively recover ordering node one by one;
        Step 2: Recover node layer by layer with fixed det(cov);
        Step 3: Recover no layer by layer, determine det(cov) adaptively
    """
    """
    1. Get index of chain components
    """

    det_cg = {}
    for idx, item in enumerate(membership):
        data_temp = []

        for col in (item):
            data_temp.append(data[:, col])

        data_temp = np.array(data_temp)
        data_temp = np.expand_dims(data_temp, axis=0)
        if data_temp.shape[0] == 1:
            det = np.var(data_temp)
        else:
            det = np.linalg.det(np.cov(data_temp))
        det_cg[idx] = det
        
    top_order = []
    for i in range(len(cc) ):

        min_cg  = min(det_cg.values())

        cc_initial = list(det_cg.keys())[list(det_cg.values()).index(min_cg)]

        del det_cg[cc_initial] 
        
        for key, val in det_cg.items():
            val -= cc_initial
            det_cg[key] = val
        top_order.append(cc_initial)
    
    
    result = []
    for i in range(len(cc)):
        result.append(top_order[i])
