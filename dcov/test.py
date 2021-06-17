#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 04:54:48 2021

@author: yohanna
"""

import matlab
import matlab.engine
import numpy as np
eng = matlab.engine.start_matlab()
Sigma = np.array([[1, 0.9], [0.9, 1]])
mat_a = matlab.double(Sigma.tolist())
a = eng.our_experiment(mat_a)
print(a)