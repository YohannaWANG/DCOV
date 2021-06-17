"""
@author: yohanna
@email: yohanna.wang0924@gmail.com
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from utils import args

def generate_adj_un():
    '''
    Data generation for unknown chain components
        Step 1: Generate undirected graph. Chain components are (1,2)(3,4),...,(n-1,n)
        Step 2: Generate directed graph over U]
    
    Notes:    Generate each chain components with size 2;
              their covariance matrix is ((a, a*squre(1-sigma)), (a*square(1-sigma), a))      
    '''
    
    n = args.n
    degree = args.a
    
    def _init_directed_adj(d: int,
                   degree: float):
       
        prob = float(degree) / (d - 1) 
        B = np.tril((np.random.rand(d, d) < prob).astype(float), k=-1)
    
        return B
    
    def _init_unidrected_adj():
        U = np.zeros((n,n))
        u_index = np.arange(start=1, stop=n+1, step=2)
    
        for i in u_index:
            j = i -1
            U[i][j] = 1
        U = U + U.T  
        return U

    '''Expand D by inserting 0 between each column and row'''
    n_d = np.int(n/2)   
    D_init = _init_directed_adj(n_d, degree)
    U_init = _init_unidrected_adj()
    
    for i in range(n_d):
        column = np.zeros(n_d)
        D_init = np.insert(D_init, i+1, column, axis=1)
        
    for i in range(n_d):
        rows = np.zeros(n)
        D_init = np.insert(D_init, i+1, rows, axis=0)         

    G = U_init + D_init

    for i in range(n):
        for j in range(n):
            if G[i][j] > 0:
                G[i][j] = 1

    '''Set A_{ij}=0 for any (i,j) pair such that i\in C(l), j\in C(m) (l>m); 
        D: Directed edges
        U: Undirected edges
    '''
    D, U = [], []  
   
    for i in range(n):
        for j in range(n):
            if G[i][j] == 1 and G[j][i] == 1:
                    U.append((j, i))
    G_temp = G - U_init
    for i in range(n):
        for j in range(n):
            if G_temp[i][j] == 1 and G_temp[j][i] == 0:
                D.append((i, j))

    return G, U_init, D, U


def generate_adj():
    '''
    Function: generated chain graph adjacency matrix.

    Arguments:
        
    Parameters:
        n = 5  # number of nodes
        c = 2  # number of chain components
        a = 2  # average node degree
        d = 10000  # number of dimensions

    Returns: 
        A, A_undirect, D, U, sorted_D, sorted_U
    '''
    
    n, c, a = args.n, args.c, args.a 

    def _init_adj():
        '''Initialize chain graph adj matrix'''
        A = np.zeros((n, n))
    
        '''Fill lower triangle'''
        p = float(a) / (n - 1)
        for j in range(n):
            for i in range(j + 1, n):
                r = np.random.uniform(0, 1)
                if r <= p: A[i][j] = 1 
        A = A + A.T
        return A
    
    def _assign_chain_comp():
        '''
        Assign chain components
            1. Select an integer c randomly as from {1,...,n} as the number of chain components;
            2. Seperate nindex np.array(len(0, n)) into R sections, forms R chain components;
        '''
        R = np.random.choice(n - 1, c, replace=False) + 1
        R = np.append(np.append(0, R), n)
        R.sort()
    
        C = {}
        V = list(range(n))
        for i in range(len(R) - 1):
            for v in V[R[i]: R[i + 1]]:
                C[v] = i
        return C
    
    A = _init_adj()
    C = _assign_chain_comp()
    
    '''
    Set A_{ij} = 0 for any (i,j) pair such that i\in C(l), j\in C(m) (l>m)
            D : Directed edges
            U : Undirected edges
    '''
    D, U = [],  []  
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

    '''
    A_undirect: undirected edges from chain graph A
        Remove directed connections and determine all undirected clusters
    '''
    A_undirect = A.copy()  
    
    for i in range(len(D)):
        j, k = D[i]
        A_undirect[(j, k)] = 0
            
    return A, A_undirect, D, U


def generate_data(G_binary, cc):
    '''
    Generate synthetic AMP chain graph data
    
    Arguments:
        G_binary : binary adj matrix
        cc       : connected chain components
    Parameters:
        n        : node number
        d        : number of samples
        det      : determinant of the covariant matrix of a chain component
    Returns:
        X        : synthetic chain graph dara
        top_sort : topological sort
    '''
    
    n, d, det = args.n, args.d, args.det
    
    def _get_weighted_adj(G_binary):
        ''' From binary G to weighted G'''
        U = np.random.uniform(low=0.5, high=1.5, size=[n, n])
        U[np.random.rand(n, n) < 0.5] *= -1
        
        G_weighted = (G_binary != 0).astype(float) * U
        
        return G_weighted 
    
    def _equal_det(n):
        import scipy.linalg as la
        '''
        Notes: Generate size n eigenvector. The product of n elements equals to 'det';
            To achieve this, we first generate 'n-1' random variables (runif) from a uniform distribution;
            The last element equals to 'det/(prod(runif))';
            Eigenvector is formed by the above n elements;
        '''

        if n > 1:
            eigen_val = np.random.uniform(low=float(0.5 / n), high=float((1.5) / n), size=n-1)
            last_val = det / np.prod(eigen_val)
            eigen_val = np.append(eigen_val, last_val)
            
        if n == 1:
            eigen_val = np.expand_dims(det, axis=0)
            
        eigen_vec = np.diag(eigen_val)
        '''
        Notes: eigen_vec = diag(λ1, ..., λn) is a diagonal matrix. 
               Q is a nxn orthogonal matrix (Compute QR decomposition of a matrix.)
        '''
        s = np.random.choice([-1, 1], size=(n, n), replace=True)
        q, _ = la.qr(np.random.rand(n, n) * s)
        semidef = q.T.dot(eigen_vec).dot(q)
    
        return semidef


    def _topological_sort(cc, G):
        '''
        get topological sort from the weighte graph
        '''
        G_weighted = G.copy()
        top_sort = []
        
        for i in cc:
            if len(i) == 1:
                if G_weighted[:,i].sum() == 0:
                    top_sort.append(i)
                    G_weighted[i,:] = 0
            else:
                if len(i) > 1:
                    j = np.setdiff1d( range(n), i)
                    if G_weighted[j , :][:,i].sum() == 0:
                        top_sort.append(i)
                        G_weighted[i,:] = 0

        if len(top_sort) != len(cc):
            print('len_top', len(top_sort))
            others = [item for item in cc if item not in top_sort]
            top_sort.extend(others)
            
        return top_sort
    
    def _predecessors(G, cc):
        '''
        get predecessors from G
        '''
        parents = []
        for nodes in cc:
            parent = np.setdiff1d((np.nonzero(G[:,nodes])), cc)
            parents.append(parent)
        return parents
    
    def _unknown_chain_sigma(size):
        '''
        sigma for unknown chain components
        eg. cc size euqas 2, cov = ((b, b*root(1-sigma)), (b*root(1-sigma), b))     
        '''
        #sigma = args.sigma
        #un_cov = np.matrix([[b, b*np.square(1-args.sigma)], [b*np.square(1-args.sigma), b]])
        un_cov = np.matrix([[2, 1.5], [1.5, 2]])
        return un_cov
        
    G_weighted = _get_weighted_adj(G_binary)
    top_sort = _topological_sort(cc, G_weighted)
    
    X = np.zeros((d, n))
    
    for cc in top_sort:
        mu = [0] * len(cc)
        
        if args.task == "known_cc":
            sigma = _equal_det(len(cc))
            inv_sigma = np.linalg.inv(sigma) 
            
        if args.task == "unknown_cc":
            sigma = _unknown_chain_sigma(len(cc))
            inv_sigma = sigma
            
        Z = np.random.multivariate_normal(mean=mu, cov=inv_sigma, size=d)
        X[:, cc] = Z

        for nodes in cc:
            a1 = np.nonzero(G_binary[:,nodes])[0]
            a2 = np.array(cc)
            parent = np.setdiff1d(a1, a2)

            if len(parent) > 0:
                M_pa = np.array((X[:, parent].dot(G_weighted[parent, nodes])))
                X[:, nodes] += M_pa

    return X, top_sort


def plot_chain_graph(cc, D, U):
    '''
    Plot chain graph
    
    Arguments:
        cc: clustered chain components
        D:  directed edges in chain graph
        U:  undirected edges in chain graph
    '''
    
    C = {}
    with open("colours.txt") as f:
        lines = f.readlines()

    for i in range(len(cc)):
        r, g, b = lines[i].strip().split(",")
        r, g, b = int(r), int(g), int(b)
        C[i] = (r, g, b)

    '''Convert to hex '''
    def clamp(x):
        return max(0, min(x, 255))

    for i in C.keys():
        r, g, b = C[i]
        h = "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))
        C[i] = h

    '''Add network vertices'''
    V = list(range(args.n)) 

    '''Assign colors to each nodes in different chain components'''
    dic_cc = {}
    for i in range(len(cc)):
        for j in cc[i]:
            dic_cc[j] = i
    nc = [C[dic_cc[v]] for v in V]

    '''Add networkx edges'''
    G = nx.DiGraph(directed=True)  
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


class Graph:
    '''
    1. Determine connected components in an undirected graph;
    2. Determine if graph contains cycles or not;
    3. Determine connected chain components
    '''
    
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]

    def DFSUtil(self, temp, v, visited):
        '''Mark the current vertex as visited'''
        visited[v] = True

        '''Store the vertex to list'''
        temp.append(v)

        for i in self.adj[v]:
            if visited[i] == False:
                temp = self.DFSUtil(temp, i, visited)
        return temp

    def addEdge(self, v, w):
        '''add an undirected edge'''
        self.adj[v].append(w)

    def connectedComponents(self):
        '''retrieve connected components in an undirected graph'''
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if visited[v] == False:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc

    def isCyclicUtil(self, v, visited, recStack):
        '''Check is the input graph contains cycles'''
        visited[v] = True
        recStack[v] = True

        for neighbour in self.adj[v]:
            if visited[neighbour] == False:
                if self.isCyclicUtil(neighbour, visited, recStack) == True:
                    return True
            elif recStack[neighbour] == True:
                return True

        recStack[v] = False
        return False

    def isCyclic(self):
        '''Returns true if graph is cyclic else false'''
        visited = [False] * self.V
        recStack = [False] * self.V
        for node in range(self.V):
            if visited[node] == False:
                if self.isCyclicUtil(node, visited, recStack) == True:
                    return True
        return False

def is_acyclic(A):
    '''
    Arguments:
        A: chain graph adjacency matrix
    Returns:
        Bool: True/ False- Whether graph contains cycles

    '''
    
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


def chain_graph_features(A_undirect):
    '''
    Arguments:
        A_undirect  : undirected adjancency matrix from chain graph A.
        
    Returns:
        cc        : nodes in each chain components;
        num_in_cc : number of nodes in each chain components.
    '''
    
    n = args.n
    g = Graph(n)
    for v in range(n):
        for w in range(n):
            if A_undirect[v, w] == 1:
                g.addEdge(v, w)
                
    '''cc: connected chain components'''           
    cc = g.connectedComponents()
    
    num_in_cc = []
    for i in range(len(cc)):
        num_in_cc.append(len(cc[i]))
        
    return cc, num_in_cc


def simulate_random_dag(d: int,
                        degree: float, 
                        graph_type: str, 
                        w_range: tuple = (0.5, 2.0)) -> nx.DiGraph:
    '''
    Simulate random DAG with some expected degree.

    Arguments:
        d: number of nodes
        degree: expected node degree, in + out
        graph_type: {erdos-renyi, barabasi-albert, full}
        w_range: weight range +/- (low, high)

    Returns:
        G: weighted DAG
    '''
    
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
    elif graph_type == 'full': 
        B = np.tril(np.ones([d, d]), k=-1)
    else:
        raise ValueError('unknown graph type')

    P = np.random.permutation(np.eye(d, d))  
    B_perm = P.T.dot(B).dot(P)
    U = np.random.uniform(low=w_range[0], high=w_range[1], size=[d, d])
    U[np.random.rand(d, d) < 0.5] *= -1
    W = (B_perm != 0).astype(float) * U

    G1 = nx.DiGraph(W)
    G2 = nx.DiGraph(B_perm)
    ''' The returned value is a binary adjacency matrix '''
    return G2


def simulate_sem(G: nx.DiGraph,
                 n: int, x_dims: int,
                 sem_type: str,
                 linear_type: str,
                 noise_scale: float = 1.0) -> np.ndarray:
    '''
    Simulate samples from SEM with specified type of noise.

    Arguments:
        G: weigthed DAG
        n: number of samples
        sem_type: {linear-gauss,linear-exp,linear-gumbel}
        noise_scale: scale parameter of noise distribution in linear SEM

    Returns:
        simu: [n,d] sample matrix
    '''
    
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


def load_dag_data(bug = False):
    '''
    load synthetic DAG data (as a special case of chain graph models)
    '''
    
    d, n = args.d, args.n
    graph_type, degree, sem_type, linear_type = args.graph_type, args.graph_degree, args.graph_sem_type, args.graph_linear_type
    x_dims = args.x_dims

    if args.data_type == 'synthetic':
        G = simulate_random_dag(n, degree, graph_type)
        simu = simulate_sem(G, d, x_dims, sem_type, linear_type)

    if args.plot:
        nx.draw_circular(G,arrowsize=15, with_labels=True)
        plt.show()    
        
    return simu, G


def load_data():
    '''
    Load synthetic chain graph data 
    
    Arguments:
        
    Returns:
        cond_cov    : Conditional covariance. COV(descendants | ancestors)    
    '''
    
    if args.graph == "known":
        '''
        chain graph data with known chain components
        '''
        print("Known Chain components")
        
        G, G_undirect, D, U = generate_adj()

    if args.graph == "unknown":
        '''
        chain graph data with unknown chain components
        '''
        print("Unkown Chain components")
        G, G_undirect, D, U = generate_adj_un()
       
    cc, num = chain_graph_features(G_undirect)
    simu, top_sort  = generate_data(G, cc)

    if args.plot:    
        plot_chain_graph(cc, D, U)    
        
    if is_acyclic(G) == True:
        raise ValueError("Contains cycles in this graph.")
             
    return simu, G, num, cc, top_sort



