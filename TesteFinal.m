clear; clc;

A = [1.2, 0.5; 
     0,   0.8];
B = [1; 0];
C = [1, 1];

[n, m] = size(B);
p = size(C, 1);

X = sdpvar(n, n, 'symmetric');
Y = sdpvar(m, p, 'full');

LMI1 = [ -X, (A*X + B*Y*C)';
          A*X + B*Y*C, -X ] <= -1e-6*eye(2*n);
     
LMI2 = X >= 1e-6*eye(n);

constraints = [LMI1, LMI2];
options = sdpsettings('solver', 'mosek', 'verbose', 0);
sol = optimize(constraints, [], options);

K = value(Y);
%%
Aa = [[A,          zeros(n,m)];
      [zeros(m,n), eye(m)]];
Ba = [[B];
      [zeros(m,m)]];
Ca = [C, eye(m)];

mu = 1.0;

sigma = 0.99;
tau = 1; %sdpvar(1,1);
kapa = 1; %sdpvar(1,1);

E = [eye(n+m),    zeros(n+m,m), zeros(n+m,m);
    zeros(m,n+m), eye(m),       zeros(m,m);
     -mu*Ba',  -eye(m),      eye(m)];

X = sdpvar(n+p, n+p, 'symmetric');
L = sdpvar(n+p, p, 'full');

P = sdpvar((n+3*m), (n+3*m), 'symmetric');

dims = n + p + 2*m;

Q = [E'*P*E,          zeros(dims, dims),    zeros(dims, m);
     zeros(dims, dims), -sigma*E'*P*E,       zeros(dims, m);
     zeros(m, dims),    zeros(m, dims),      -tau*eye(m)];

F1 = sdpvar(n+m, m, 'full'); H1 = sdpvar(n+m, m, 'full');
F2 = sdpvar(m, m, 'full');   H2 = sdpvar(m, m, 'full');
F3 = sdpvar(m, m, 'full');   H3 = sdpvar(m, m, 'full');
F4 = sdpvar(n+m, m, 'full'); H4 = sdpvar(n+m, m, 'full');
F5 = sdpvar(m, m, 'full');   H5 = sdpvar(m, m, 'full');
F6 = sdpvar(m, m, 'full');   H6 = sdpvar(m, m, 'full');
F7 = sdpvar(m, m, 'full');   H7 = sdpvar(m, m, 'full');

X_cal = [kapa*eye(n+m), F1, H1;
         zeros(m,n+m),     F2, H2;
         zeros(m,n+m),     F3, H3;
         eye(n+m),         F4, H4;
         zeros(m,n+m),     F5, H5;
         zeros(m,n+m),     F6, H6;
         zeros(m,n+m),     F7, H7];

B_cal = [-X',          zeros(n+m,m),  mu*Ba;
         zeros(m,n+m),   -eye(m),     eye(m);
         zeros(m,n+m),   zeros(m,m),  -eye(m);
         (X*Aa - L*Ca)', -mu*Aa'*Ba, zeros(n+m, m);
         zeros(m, n+m),  zeros(m,m),  zeros(m, m);
         Ba'*X',        -mu*(Ba'*Ba), eye(m);
         zeros(m,n+m),     eye(m),    zeros(m, m)]';

LMI = Q + X_cal * B_cal + B_cal' * X_cal';

constraints = [LMI <= -1e-6*eye(2*(n+m)+5*m), P >= 1e-6*eye(n+p+2*m)];  % X definida positiva

Objective = sigma + tau + kapa;
options = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1);
sol = optimize(constraints, [], options);

X = value(X);
L = value(L);

L = X*L