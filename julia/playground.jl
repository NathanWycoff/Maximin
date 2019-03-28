using JuMP
using Gurobi

# Define domain
N = 100;
M = 10;#precision of estimate
P = 2;

# Define distance matrices
D = [];
for n1 = 1:(N-1)
    push!(D, []);
    for n2 = (n1+1):N
        Dij = zeros(N*P, N*P);

        start1 = (n1-1)*P + 1;
        start2 = (n2-1)*P + 1;
        Dij[(start1):(start1+P-1),(start1):(start1+P-1)] = eye(P);
        Dij[(start2):(start2+P-1),(start2):(start2+P-1)] = eye(P);
        Dij[(start1):(start1+P-1),(start2):(start2+P-1)] = -eye(P);
        Dij[(start2):(start2+P-1),(start1):(start1+P-1)] = -eye(P);

        push!(D[n1], Dij);
    end
end

# Set up optimization problem.
maximin = Model(solver=GurobiSolver())

@variable(maximin, 0 <= x[n=1:(N*P)] <= 1);
@variable(maximin, b[n=1:(N*P),1:M], Bin);
@variable(maximin, t >= 0);

@objective(maximin, Max, t);

#for n1 = 1:(N-1)
#    for n2 = (n1+1):N
#        @constraint(maximin, t <= x'*D[n1][n2-n1]*x)
#    end
#end

for n1 = 1:(N-1)
    for n2 = (n1+1):N
        @constraint(maximin, t <= sum(0.5^m * b[n1,m] for m=1:M)^2 - 2 * sum(0.5^m * b[n1,m] for m=1:M) *sum(0.5^m * b[n2,m] for m=1:M) + sum(0.5^m * b[n2,m] for m=1:M)^2)
    end
end

for i = 1:(N*P)
    @constraint(maximin, sum(0.5^m * b[i,m] for m=1:M) == x[i]) 
end

solve(maximin)
getvalue(x)
getvalue(t)
