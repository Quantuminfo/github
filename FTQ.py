# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import cvxpy as cvx
import numpy as np

x1,x2=cvx.Variable(),cvx.Variable()
constrains=[2*x1+x2>=1,x1+ 3*x2>=1,x1>=0,x2>=0]
f=2*x1 + 3*x2
objective=cvx.Minimize(f)
problem=cvx.Problem(objective,constrains)
problem.solve()
print(f.value,x1.value,x2.value)
