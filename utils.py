"""
@author: yohanna
@email: yohanna.wang0924@gmail.com
"""

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--graph', type=str, default="unknown", 
                    choices=['unknown', 'known'], 
                    help="generate chain graph data or DAG data")
parser.add_argument('--task', type=str, default="unknown_cc", choices=['known_cc', 'unknown_cc'], 
                    help="known_cc: Known cc with equal determinant; \
                        unknown_cc: Unknown cc with specialized cov")

parser.add_argument('--algorithm', type=str, default="unknown_dcov", choices=['known_dcov','unknown_dcov'], 
                    help="known_dcov: algorithm for known cc with equal determinant;\
                          unknown_dcov: unknown cc with specialized cov")
"""
Notes: Graph: chain_graph_known_cc; Task: eql_det;
       Graph: chain_graph_un_cc;    Task: unknown_cc
"""
parser.add_argument('--regress_method', type=str, default="mgcv",
                    choices=['np', 'mgcv'], 
                    help="regression method (calculate regression residual OR prune)")
parser.add_argument('--sfo', type=str, default='sfo-python', choices = ['sfo-matlab', 'sfo-python'],
                    help="Use either Matlab or Python implementation of SUbmodular function minimization")
parser.add_argument('--plot', action='store_true', default=True,
                    help='plot graph')

parser.add_argument('--det', type=int, default=1,
                    help='equal determinant for chain components')
parser.add_argument('--n', type=int, default=20, help='number of nodes')
parser.add_argument('--c', type=int, default=2, help='number of chain components')
parser.add_argument('--a', type=int, default=2, help='average node degree')
parser.add_argument('--d', type=int, default=1000, help="number of dimensions")


"""Parameters for identifiability of unknown chain components"""
parser.add_argument('--b', type=int, default=1, help='correlations of chain components')
parser.add_argument('--sigma', type=float, default=0.1, help='sigma for cov of chain comp')


parser.add_argument('--operator', type=str, default='det', 
help="(TODO)super-additive function on positive semidefinite matrices. 1.det, 2.trace")
parser.add_argument('--eta', type=float, default=0.01, help="residual error variance")


parser.add_argument('--graph_type', type=str, default='erdos-renyi',
                    help='the type of DAG graph by generation method')
parser.add_argument('--graph_degree', type=int, default=2,
                    help='the number of degree in generated DAG graph')
parser.add_argument('--graph_sem_type', type=str, default='linear-gauss',
                    help='the structure equation model (SEM) parameter type')
parser.add_argument('--graph_linear_type', type=str, default='linear',
                    help='the synthetic data type: linear -> linear SEM, nonlinear_1 -> x=Acos(x+1)+z, nonlinear_2 -> x=2sin(A(x+0.5))+A(x+0.5)+z')
args = parser.parse_args()


