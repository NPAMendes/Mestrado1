function [K, sol] = opt_control(A, B, eps)
%OPT_CONTROL Summary of this function goes here
%   Detailed explanation goes here

nx = size(A, 1);
nu = size(B,2);

P = sdpvar(nx, nx, 'symmetric');
X = sdpvar(nu, nx, 'full');

LMI1 = [P           (A*P + B*X)';
        A*P + B*X   P            ] >= eps*eye(2*nx);
     
LMI2 = P >= eps*eye(nx);

constraints = [LMI1, LMI2];
options = sdpsettings('solver', 'mosek', 'verbose', 0);

sol = optimize(constraints, [], options);

if sol.problem == 0
    K = value(X)/value(P);
else
    disp('Optimization failed');
    exit;
end

end