#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:09:39 2021

@author: 
"""

import argparse

parser = argparse.ArgumentParser()
# ----------- data parameters --------
parser.add_argument('--graph', type=str, default="chain_graph", 
                    choices=['dag_graph', 'chain_graph', 'npvar_dag'], 
                    help="generate chain graph data or DAG data")
parser.add_argument('--regress_method', type=str, default="np",
                    choices=['np', 'mgcv'], 
                    help="regression method (calculate regression residual OR prune)")
parser.add_argument('--plot', action='store_true', default=True,
                    help='plot graph')

parser.add_argument('--det', type=int, default=1,
                    help='equal determinant for chain components')
parser.add_argument('--n', type=int, default=20, help='number of nodes')
parser.add_argument('--c', type=int, default=6, help='number of chain components')
parser.add_argument('--a', type=int, default=2, help='average node degree')
parser.add_argument('--d', type=int, default=100, help="number of dimensions")

parser.add_argument('--operator', type=str, default='det', 
help="super-additive function on positive semidefinite matrices. 1.det, 2.trace")
parser.add_argument('--eta', type=float, default=0.01, help="residual error variance")



parser.add_argument('--data_type', type=str, default= 'synthetic',
                    choices=['synthetic', 'discrete', 'real', 'missing'],
                    help='choosing which experiment to do.')
parser.add_argument('--graph_type', type=str, default='erdos-renyi',
                    help='the type of DAG graph by generation method')
parser.add_argument('--graph_degree', type=int, default=2,
                    help='the number of degree in generated DAG graph')
parser.add_argument('--x_dims', type=int, default=1, 
                    help='The number of input dimensions: default 1.')
parser.add_argument('--graph_sem_type', type=str, default='linear-gauss',
                    help='the structure equation model (SEM) parameter type')
parser.add_argument('--graph_linear_type', type=str, default='linear',
                    help='the synthetic data type: linear -> linear SEM, nonlinear_1 -> x=Acos(x+1)+z, nonlinear_2 -> x=2sin(A(x+0.5))+A(x+0.5)+z')
args = parser.parse_args()