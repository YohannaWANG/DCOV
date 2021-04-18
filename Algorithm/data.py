#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 08:35:27 2021

@author: 
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt 

from utils import args


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
    p = float(a) / (n - 1)
    for j in range(n):
        for i in range(j + 1, n):
            r = np.random.uniform(0, 1)
            if r <= p: A[i][j] = 1  
    A = A + A.T

    """Assign chain components"""
    # Select an integer c randomly as from {1,...,n} as the number of chain components;
    # Seperate nindex np.array(len(0, n)) into R sections, forms R chain components;
    R = np.random.choice(n - 1, c, replace=False) + 1
    R = np.append(np.append(0, R), n)
    R.sort()

    C = {}
    V = list(range(n))
    for i in range(len(R) - 1):
        for v in V[R[i]: R[i + 1]]:
            C[v] = i

    # Set A_{ij}=0 for any (i,j) pair such that i\in C(l), j\in C(m) (l>m);
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
    for v in range(n):
        for w in range(n):
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
    n = A.shape[0]
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
    V = list(range(args.n))  # Nodes indices

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
        mu, sigma = [0] * (cc_section[j + 1] - cc_section[j]), B[cc_section[j]:cc_section[j + 1], cc_section[j]:cc_section[j + 1]]
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
        X[cc_section[j]:cc_section[j + 1]] = np.random.multivariate_normal(mean=mu, cov=sigma, size=d).T

    return np.linalg.solve(B, X).T

"""
    Function: generate DAG data and test the performance of our algorithm
"""

def load_dag_data(bug = False):
    #  # configurations
    d, n = args.d, args.n
    graph_type, degree, sem_type, linear_type = args.graph_type, args.graph_degree, args.graph_sem_type, args.graph_linear_type
    x_dims = args.x_dims

    if args.data_type == 'synthetic':
        # generate data
        G = simulate_random_dag(n, degree, graph_type)
        X = simulate_sem(G, d, x_dims, sem_type, linear_type)

    return X, G

#------------------------------------------------------------------------------
def simulate_random_dag(d: int, #d=10
                        degree: float, #degree=3
                        graph_type: str, 
                        w_range: tuple = (0.5, 2.0)) -> nx.DiGraph:
    """Simulate random DAG with some expected degree.
        Link: https://github.com/xunzheng/notears
    Args:
        d: number of nodes
        degree: expected node degree, in + out
        graph_type: {erdos-renyi, barabasi-albert, full}
        w_range: weight range +/- (low, high)

    Returns:
        G: weighted DAG
    """
    if graph_type == 'erdos-renyi':
        prob = float(degree) / (d - 1)
        B = np.tril((np.random.rand(d, d) < prob).astype(float), k=-1)
    elif graph_type == 'barabasi-albert':
        m = int(round(degree / 2))
        B = np.zeros([d, d])
        bag = [0]
        for ii in range(1, d):
            dest = np.random.choice(bag, size=m)
            for jj in dest:
                B[ii, jj] = 1
            bag.append(ii)
            bag.extend(dest)
    elif graph_type == 'full':  # ignore degree, only for experimental use
        B = np.tril(np.ones([d, d]), k=-1)
    else:
        raise ValueError('unknown graph type')
    # random permutation
    P = np.random.permutation(np.eye(d, d))  # permutes first axis only
    B_perm = P.T.dot(B).dot(P)
    U = np.random.uniform(low=w_range[0], high=w_range[1], size=[d, d])
    U[np.random.rand(d, d) < 0.5] *= -1
    W = (B_perm != 0).astype(float) * U

    G1 = nx.DiGraph(W)
    G2 = nx.DiGraph(B_perm)
    """
    The returned value is a binary adjacency matrix
    """
    return G2


def simulate_sem(G: nx.DiGraph,
                 n: int, x_dims: int,
                 sem_type: str,
                 linear_type: str,
                 noise_scale: float = 1.0) -> np.ndarray:
    """Simulate samples from SEM with specified type of noise.

    Args:
        G: weigthed DAG
        n: number of samples
        sem_type: {linear-gauss,linear-exp,linear-gumbel}
        noise_scale: scale parameter of noise distribution in linear SEM

    Returns:
        X: [n,d] sample matrix
    """
    W = nx.to_numpy_array(G)
    d = W.shape[0]
    X = np.zeros([n, d, x_dims])
    ordered_vertices = list(nx.topological_sort(G))
    assert len(ordered_vertices) == d
    for j in ordered_vertices:
        parents = list(G.predecessors(j))
        if linear_type == 'linear':
            eta = X[:, parents, 0].dot(W[parents, j])
        elif linear_type == 'nonlinear_1':
            eta = np.cos(X[:, parents, 0] + 1).dot(W[parents, j])
        elif linear_type == 'nonlinear_2':
            eta = (X[:, parents, 0]+0.5).dot(W[parents, j])
        else:
            raise ValueError('unknown linear data type')

        if sem_type == 'linear-gauss':
            if linear_type == 'linear':
                X[:, j, 0] = eta + np.random.normal(scale=noise_scale, size=n)
            elif linear_type == 'nonlinear_1':
                X[:, j, 0] = eta + np.random.normal(scale=noise_scale, size=n)
            elif linear_type == 'nonlinear_2':
                X[:, j, 0] = 2.*np.sin(eta) + eta + np.random.normal(scale=noise_scale, size=n)
        elif sem_type == 'linear-exp':
            X[:, j, 0] = eta + np.random.exponential(scale=noise_scale, size=n)
        elif sem_type == 'linear-gumbel':
            X[:, j, 0] = eta + np.random.gumbel(scale=noise_scale, size=n)
        else:
            raise ValueError('unknown sem type')
    if x_dims > 1 :
        for i in range(x_dims-1):
            X[:, :, i+1] = np.random.normal(scale=noise_scale, size=1)*X[:, :, 0] + np.random.normal(scale=noise_scale, size=1) + np.random.normal(scale=noise_scale, size=(n, d))
        X[:, :, 0] = np.random.normal(scale=noise_scale, size=1) * X[:, :, 0] + np.random.normal(scale=noise_scale, size=1) + np.random.normal(scale=noise_scale, size=(n, d))
    X = np.squeeze(X)
    return X


