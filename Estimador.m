function [L,sol, gamma, alpha] = Estimador(A, C, D)
    nx = size(A,1);
    ny = size(C, 1);
    ns = size(D, 2);

    Aa = [A,             zeros(nx, ns);
          zeros(ns, nx), eye(ns)];

    Ca = [C, D];

    Ma = [zeros(nx, ns);
          eye(ns)];

    P = sdpvar(nx+ns, nx+ns, 'symmetric');
    Y = sdpvar(nx+ns, ny, 'full');

    gamma = sdpvar(1,1,'full');
    s = sdpvar(1,1,'full');

    Q = [P - eye(nx+ns),  zeros(nx+ns,ns), Aa'*P - Ca'*Y';
         zeros(ns,nx+ns), gamma*eye(ns),   Ma'*P;
         P*Aa - Y*Ca,     P*Ma,            P];

    constraints = [Q >= 1e-6*eye(2*nx+3*ns), P >= 1e-6*eye(nx+ns), gamma >= 0, s >= 0.1, trace(P)==s];
    obj = gamma;
    options = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(constraints, obj, options);

    s = value(s);
    gamma = value(gamma)/s;
    alpha = 1/s;
    P = value(P)/s;
    Y = value(Y)/s;
    L = inv(P)*Y;
    
end