function [K,sol, alpha] = Control(A, B)
    nx = size(A,1);
    
    P = sdpvar(nx, nx, 'symmetric');
    Y = sdpvar(nx, nx, 'full');
    alpha = sdpvar(1,1,'full');

    Q = [P-alpha*eye(nx), A'*P+Y';
         P*A + Y,         P];

    LMI_K = [Q >= 1e-6*eye(2*nx), P>=1e-6*eye(nx), trace(P)==1, alpha>=0];

    Objective = -alpha;
    options = sdpsettings('solver','mosek','verbose',0);
    sol = optimize(LMI_K, Objective, options);

    alpha = value(alpha);
    P = value(P);
    Y = value(Y);
    K = inv(B'*B)*B'*Y*inv(P);
end