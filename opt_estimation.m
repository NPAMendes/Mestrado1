function [Lt, L, sol] = opt_estimation(A,B,C,D,sigma,kappa,mu,eps)

%% Definição das dimensões do sistema
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);
ns = size(D,2);

Aa = [A, zeros(nx, ns); zeros(ns, nx), eye(ns)];
Ba = [B; zeros(ns, nu)];
Ca = [C, D];
Ma = [zeros(nx, ns); eye(ns)];

%% Declaração de variáveis sdp
tau = sdpvar(1,1,'full');

X = sdpvar(nx+ns, nx+ns, 'full');
Lt = sdpvar(nx+ns, ny, 'full');
P = sdpvar(nx+ns+2*nu, nx+ns+2*nu, 'symmetric');

%% Definição das matrizes do Finsler

E = [eye(nx+ns),     zeros(nx+ns,nu), zeros(nx+ns,nu);
    zeros(nu,nx+ns), eye(nu),         zeros(nu,nu);
     -mu*Ba',        -eye(nu),        eye(nu)];

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

X_cal = [kappa*eye(nx+ns), F1, H1;
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
         (X*Aa - Lt*Ca)',   -mu*Aa'*Ba,      zeros(nx+ns, nu);
         zeros(nu, nx+ns), zeros(nu,nu),    zeros(nu, nu);
         Ba'*X',           -mu*(Ba'*Ba),    eye(nu);
         zeros(nu,nx+ns),  eye(nu),         zeros(nu, nu);
         Ma'*X',           zeros(ns, nu),   zeros(ns, nu)]';

LMI = Q + X_cal * B_cal + B_cal' * X_cal';

constraints = [LMI <= -eps*eye(2*nx+3*ns+5*nu), P >= eps*eye(nx+ns+2*nu), tau >= eps];

obj = tau;
options = sdpsettings('solver', 'mosek', 'verbose', 0);
sol = optimize(constraints, obj, options);

if sol.problem == 0
    Lt = value(Lt);
    L = value(X)\value(Lt);
else
    disp('Optimization failed');
    exit;
end

end