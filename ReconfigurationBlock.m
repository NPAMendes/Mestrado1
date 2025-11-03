function [R1, R2, R3, R4, sol, gamma, alpha] = ReconfigurationBlock(A, B, C, D, K, L, sigma)
    % Parâmetros do sistema
    nx = size(A, 1);
    nu = size(B, 2);
    ny = size(C, 1);
    ns = size(D,2);

    % Matrizes do sistema com falha
    Bf = B*diag(1-sigma);

    % Matrizes do controlador    
    L1 = L(1:nx, :);
    L2 = L(nx+1:end, :);

    Ac = [A + B*K - L1*C, -L1*D;
          -L2*C,          eye(ns) - L2*D];

    Bc = L;

    Cc = [K, zeros(nu, ns)];

    % Matrizes aumentadas para reconfiguração
    Ar = [A,               zeros(nx,nx+ns);
          zeros(nx+ns,nx), Ac];

    Br = [zeros(nx,ny), Bf;
          Bc,           zeros(nx+ns,nu)];

    Cr = [C,            zeros(ny, nx+ns);
          zeros(nu,nx), Cc];

    Dr = [D;
          zeros(nu,ns)];
    
    % Matrizes da LMI
    R1 = sdpvar(ny,ny, 'full');
    R2 = sdpvar(ny,nu, 'full');
    R3 = sdpvar(nu,ny,'full');
    R4 = sdpvar(nu,nu,'full');

    R = [R1, R2;
         R3, R4];
    
    P = sdpvar(2*nx+ns,2*nx+ns,'symmetric');
    
    gamma = sdpvar(1,1,'full');
    alpha = sdpvar(1,1,'full');

    B_cal = [Ar,                   -eye(2*nx+ns),        zeros(2*nx+ns,ny+nu), Br,                 zeros(2*nx+ns,ns);
             Cr,                   zeros(ny+nu,2*nx+ns), -eye(ny+nu),          zeros(ny+nu,ny+nu), Dr;
             zeros(ny+nu,2*nx+ns), zeros(ny+nu,2*nx+ns), R,                    -eye(ny+nu),        zeros(ny+nu,ns)];

    Q = [-P + alpha*eye(2*nx+ns), zeros(2*nx+ns,2*nx+ns), zeros(2*nx+ns,ny+nu), zeros(2*nx+ns,ny+nu), zeros(2*nx+ns,ns);
         zeros(2*nx+ns,2*nx+ns),  P,                      zeros(2*nx+ns,ny+nu), zeros(2*nx+ns,ny+nu), zeros(2*nx+ns,ns);
         zeros(ny+nu,2*nx+ns),    zeros(ny+nu,2*nx+ns),   zeros(ny+nu,ny+nu),   zeros(ny+nu,ny+nu),   zeros(ny+nu,ns);
         zeros(ny+nu,2*nx+ns),    zeros(ny+nu,2*nx+ns),   zeros(ny+nu,ny+nu),   zeros(ny+nu,ny+nu),   zeros(ny+nu,ns);
         zeros(ns,2*nx+ns),       zeros(ns,2*nx+ns),      zeros(ns,ny+nu),      zeros(ns,ny+nu),      -gamma*eye(ns)];

    F1 = sdpvar(2*nx+ns,2*nx+ns, 'full');
    F2 = sdpvar(2*nx+ns,2*nx+ns, 'full');
    F3 = sdpvar(ny+nu,2*nx+ns, 'full');
    F4 = sdpvar(ny+nu,2*nx+ns, 'full');
    F5 = sdpvar(ns,2*nx+ns, 'full');

    H1 = sdpvar(2*nx+ns,ny+nu, 'full');
    H2 = sdpvar(2*nx+ns,ny+nu, 'full');
    H3 = sdpvar(ny+nu,ny+nu, 'full');
    H4 = sdpvar(ny+nu,ny+nu, 'full');
    H5 = sdpvar(ns,ny+nu, 'full');

    X_cal = [F1, H1, zeros(2*nx+ns,ny+nu);
             F2, H2, zeros(2*nx+ns,ny+nu);
             F3, H3, eye(ny+nu);
             F4, H4, eye(ny+nu);
             F5, H5, zeros(ns,ny+nu);];

    LMI = Q + X_cal*B_cal + B_cal'*X_cal';

    constraints = [LMI <= -1e-6*eye(4*nx+3*ns+2*ny+2*nu), P>=1e-6*eye(2*nx+ns), alpha>=0, gamma>=0, trace(P) == 1];
    options = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug', 1);
    sol = optimize(constraints, [], options);
    
    alpha = value(alpha);
    gamma = value(gamma);
    R1 = value(R1);
    R2 = value(R2);
    R3 = value(R3);
    R4 = value(R4);
end