clear; clc;

A = [1.2, 0.5; 
     0,   0.8];
B = [1; 0];
C = [1, 0];
D = 0;
%%
syms z

nx = size(A,1);

M = [A - z*eye(nx), B; 
     C,             zeros(size(C,1), size(B,2))];

condicao = simplify(det(M'*M));

zeros_invariantes = solve(condicao == 0, z);

%%
nu = size(B,2);
ny = size(C, 1);
ns = size(D,2);

Q = sdpvar(nx, nx, 'symmetric');
Y = sdpvar(nu, nx, 'full');

LMI1 = [Q,        (A*Q + B*Y)';
        A*Q + B*Y, Q ] >= 1e-6*eye(2*nx);
     
LMI2 = Q >= 1e-6*eye(nx);

constraints = [LMI1, LMI2];
options = sdpsettings('solver', 'mosek', 'verbose', 0);
sol = optimize(constraints, [], options);

K = Y*(C*value(Q))'*inv(C*value(Q)*(C*value(Q))');
%%
Aa = [[A,            zeros(nx,ns)];
      [zeros(ns,nx), eye(ns)]];
Ba = [B;
      zeros(ns,nu)];
Ca = [C, eye(ns)];
Ma = [[zeros(nx, ns)];
    [eye(ns)]];

mu = 1;

sigma = 1;
tau = 1;
kapa = 1;

E = [eye(nx+ns),     zeros(nx+ns,nu), zeros(nx+ns,nu);
    zeros(nu,nx+ns), eye(nu),         zeros(nu,nu);
     -mu*Ba',        -eye(nu),        eye(nu)];

X = sdpvar(nx+ns, nx+ns, 'full');
L = sdpvar(nx+ns, ns, 'full');

P = sdpvar(nx+ns+2*nu, nx+ns+2*nu, 'symmetric');

Q = [E'*P*E,                        zeros(nx+ns+2*nu, nx+ns+2*nu), zeros(nx+ns+2*nu, nu+ns);
     zeros(nx+ns+2*nu, nx+ns+2*nu), -sigma*E'*P*E,                 zeros(nx+ns+2*nu, nu+ns);
     zeros(nu+ns, nx+ns+2*nu),      zeros(nu+ns, nx+ns+2*nu),      -tau*eye(nu+ns)];

F1 = sdpvar(nx+ns, nu, 'full'); H1 = sdpvar(nx+ns, nu, 'full');
F2 = sdpvar(nu, nu, 'full');    H2 = sdpvar(nu, nu, 'full');
F3 = sdpvar(nu, nu, 'full');    H3 = sdpvar(nu, nu, 'full');
F4 = sdpvar(nx+ns, nu, 'full'); H4 = sdpvar(nx+ns, nu, 'full');
F5 = sdpvar(nu, nu, 'full');    H5 = sdpvar(nu, nu, 'full');
F6 = sdpvar(nu, nu, 'full');    H6 = sdpvar(nu, nu, 'full');
F7 = sdpvar(nu, nu, 'full');    H7 = sdpvar(nu, nu, 'full');
F8 = sdpvar(nu, nu, 'full');    H8 = sdpvar(nu, nu, 'full');

X_cal = [kapa*eye(nx+ns), F1, H1;
         zeros(nu,nx+ns), F2, H2;
         zeros(nu,nx+ns), F3, H3;
         eye(nx+ns),      F4, H4;
         zeros(nu,nx+ns), F5, H5;
         zeros(nu,nx+ns), F6, H6;
         zeros(nu,nx+ns), F7, H7;
         zeros(nu,nx+ns), F8, H8];

B_cal = [-X',              zeros(nx+ns,nu), mu*Ba;
         zeros(nu,nx+ns),  -eye(nu),        eye(nu);
         zeros(nu,nx+ns),  zeros(nu,nu),    -eye(nu);
         (X*Aa - L*Ca)',   -mu*Aa'*Ba,      zeros(nx+ns, nu);
         zeros(nu, nx+ns), zeros(nu,nu),    zeros(nu, nu);
         Ba'*X',           -mu*(Ba'*Ba),    eye(nu);
         zeros(nu,nx+ns),  eye(nu),         zeros(nu, nu);
         Ma'*X',           zeros(ns, nu),   zeros(ns, nu)]';

LMI = Q + X_cal * B_cal + B_cal' * X_cal';

constraints = [LMI <= -1e-6*eye(2*nx+3*ns+5*nu), P >= 1e-6*eye(nx+ns+2*nu)];

Objective = sigma + tau + kapa;
options = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1);
sol = optimize(constraints, [], options);

X = value(X);
L = value(L);

L = X*L