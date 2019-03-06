#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  python/geom_form.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 03.06.2019

# See if we can get a geometric formulation of the problem.
import cvxpy as cp
import numpy as np

N = 3

ex = cp.Variable(N, pos = True, name = "exp(x)")
et = cp.Variable(pos = True, name = "exp(t)")

obj = pow(et, -1)

# Create linearizing constraints
constr = []
for n1 in range(N-1):
    for n2 in range(n1+1, N):
        constr += [et * pow(ex[n1], -1) * ex[n2] <= 1]
        constr += [pow(et, -1) * ex[n1] * pow(ex[n2], -1) <= 1]

# Create bound constraints
for n in range(N):
    constr += [ex[n] >= 1, ex[n] <= np.e]

maximin = cp.Problem(cp.Minimize(obj), constr)

result = maximin.solve(gp = True)

et.value
ex.value

################################################

ex1 = cp.Variable(pos = True, name = "x1")
ex2 = cp.Variable(pos = True, name = "x2")
et = cp.Variable(pos = True, name = "t")

obj = pow(et, -1)
constr1 = et * pow(ex1, -1) * ex2 <= 1
constr2 = pow(et, -1) * ex1 * pow(ex2, -1) <= 1

lb1 = ex1 >= 1
lb2 = ex2 >= 1

ub1 = ex1 <= np.e
ub2 = ex2 <= np.e


print(obj.is_dgp())
print(constr1.is_dgp())
print(constr2.is_dgp())
print(lb1.is_dgp())
print(lb2.is_dgp())
print(ub1.is_dgp())
print(ub2.is_dgp())

constraints = [constr1, constr2, lb1, lb2, ub1, ub2]

maximin = cp.Problem(cp.Minimize(obj), constraints)

maximin.is_dgp()

result = maximin.solve(gp = True)

et.value
