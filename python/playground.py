import cvxpy as cvx
import dccp

N = 2
P = 2
X = cvx.Variable(N*P)
t = cvx.Variable()

# Add maximin constraints
constr = []
for n1 in range(N-1):
    for n2 in range(n1+1, N):
        #constr += [t <= cvx.norm(X[n1,:] - X[n2,:], 2)]
        constr += [t <= cvx.norm(X[(n1*P):((n1+1)*P)] - X[(n2*P):((n2+1)*P)], 2)]
        #constr += [t <= X[n1,:]]

# Add bounding box
for n in range(N):
    for p in range(P):
        constr += [X[n*P + p] >= 0, X[n*P + p] <= 1]
        pass
constr += [t >= 0, t <= 1]

maximin = cvx.Problem(cvx.Maximize(t), constr)

print("problem is DCP:", maximin.is_dcp())   # false
print("problem is DCCP:", dccp.is_dccp(maximin))  # true

result = maximin.solve(method = 'dccp')

X.value
t.value
