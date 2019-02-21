using JuMP
using Ipopt

# Define domain
N = 5;
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
maximin = Model(solver=IpoptSolver())

@variable(maximin, 0 <= x[i=1:(N*P)] <= 1);
@variable(maximin, t >= 0);

@objective(maximin, Max, t);

for n1 = 1:(N-1)
    for n2 = (n1+1):N
        @constraint(maximin, t <= x'*D[n1][n2-n1]*x)
    end
end

solve(maximin)
