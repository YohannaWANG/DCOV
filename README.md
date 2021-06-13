
<img align="left" src="docs/images/dcov.png"> &nbsp; &nbsp;


                                                               
# Identifiability of AMP Chain Graph Models

## Background
Chain graph models contain both <u>directed and undirected edges</u> and can be used to represent both **association** and **causation** in real-world applications

### Chain graph models

<img width="" alt="characterization" src="/docs/images/chain_graph.png"/>         

### Identifiability

<p align="center">
<img width="500" alt="characterization" src="/docs/images/Identifiability.png"/>         
</p>

### Video explanation 
- About chain graph models.
- About this work

## Introduction
We study **identifiability** of **Andersson-Madigan-Perlman (AMP)** chain graph models, which are a common generalization of linear structural equation models and Gaussian graphical models. <u>AMP models are described by DAGs on chain components which themselves are undirected graphs.</u> 

For a known chain component decomposition, we show that the DAG on the chain components is identifiable if the determinants of the residual covariance matrices of the chain components are monotone non-decreasing in topological order. This condition extends the equal variance identifiability criterion for Bayes nets, and it can be generalized from determinants to any super-additive function on positive semidefinite matrices. When the component decomposition is  unknown, we describe conditions that allow recovery of the full structure using a polynomial time algorithm based on submodular function minimization. We also conduct experiments comparing our algorithm's performance against existing baselines.                                       

                                                            
## Prerequisites

- **Python 3.6+**
  - `networkx`
  - `argpase`
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `torch`: Optional, only used for nonlinear model.
- **R 4.0.0**
  - `rpy2`: Python interface, enables calling R from Python. Install [rpy2](https://pypi.org/project/rpy2/) first.
  - `pcalg` : [Methods for Graphical Models and Causal Inference](https://cran.r-project.org/web/packages/pcalg/index.html)
  - `mgvc` : [Mixed GAM Computation Vehicle with Automatic Smoothness Estimation](https://cran.r-project.org/web/packages/mgcv/index.html)
  - `ggm` : [Graphical Markov Models with Mixed Graphs](https://cran.r-project.org/web/packages/ggm/index.html)
  - `lcd` :[Learn Chain graphs via Decomposition](http://www2.uaem.mx/r-mirror/web/packages/lcd/index.html) Notes: [lcd](http://www2.uaem.mx/r-mirror/web/packages/lcd/index.html) is a relatively old package and can be installed from [here](http://www2.uaem.mx/r-mirror/src/contrib/lcd_0.7-3.tar.gz).
- **Matlab R2020b**
  - `Python-Matlab` : Calling Matlab from Python. Install [Python-Matlab](https://www.mathworks.com/help/matlab/matlab-engine-for-python.html) first.
  - `SFO` : A Matlab toolbox for Submodular Function Optimization [SFO(v 2.0)](https://www.mathworks.com/matlabcentral/fileexchange/20504-submodular-function-optimization).


## Contents

- `data.py` - generate synthetic chain graph data, including graph simulation and data simulation
- `evaluate.py` - algorithm accuracy evaluation 
- `utils.py` - simulation parameters, such as selecte graph type, node number, data type, graph degree, etc.  
- `utils.R` - wrapper for scipy's LBFGS-B
- `utils.py` - graph simulation, data simulation, and accuracy evaluation 


## Running a simple demo

The simplest way to try out DCOV is to run a simple example:
```bash
$ git clone https://github.com/YohannaWANG/DCOV.git
$ cd DCOV/
$ python DCOV/demo.py
```

## Runing as a command

Alternatively, if you have a CSV data file `X.csv`, you can install the package and run the algorithm as a command:
```bash
$ pip install git+git://github.com/YohannaWANG/DCOV
$ cd DCOV
$ python main.py --graph chain_graph_known_cc --task eql_det --algorithm known_npcov --regress_method mgcv --n 50 --s 1000 --d 4 --operator det
```

## Algorithms

- ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)  **Algorithm 1** AMP Chain graph identification with known chain components.
   <img width="800" align ="center" alt="characterization" src="/docs/images/algo1.png" >
- ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)  **Algorithm 2**  AMP Chain graph identification from unknown chain components. 
   <img width="800" align ="center" alt="characterization" src="/docs/images/algo2.png">    
   
## Performance


## Citation

## Contacts

Please feel free to contact me if you meet any problem when using this code. I'm glad to hear other advise and update our work. 
I am also open to collaboration if you think that you are working on a problem that I might be interested in it.
Please do not hestitate to contact us!


