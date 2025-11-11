function [R3, R4, sol, gamma, alpha] = ReconfigurationBlock(A, B, C, D, K, sigma)
    % Parâmetros do sistema
    nx = size(A, 1);
    nu = size(B, 2);
    ny = size(C, 1);
    ns = size(D,2);

    % Matrizes do sistema com falha
    Bf = B*diag(1-sigma);

    % Matrizes aumentadas
    Aa = [A,            zeros(nx,ns);
          zeros(ns,nx), eye(ns)];

    Baf = [Bf;
           zeros(ns,nu)];

    Cs = [C, zeros(ny,ns)];

    Ka = [K, zeros(nu,ns)];

    Ma = [zeros(nx,ns);
          eye(ns)];

    % Variáveis a serem encontradas
    P = sdpvar(nx+ns, nx+ns, 'symmetric');
    Y = sdpvar(nx+ns, ny, 'full');
    G = sdpvar(nx+ns, nu, 'full');
    gamma = sdpvar(1, 1, 'full');
    s = sdpvar(1, 1, 'full');

    LMI = [P + eye(nx+ns),     zeros(nx+ns, ns), Aa'*P + Cs'*Y' + Ka'*G';
           zeros(ns, nx+ns),   gamma*eye(ns),    Ma'*P;
           P*Aa + Y*Cs + G*Ka, P*Ma,             P];

    constraints = [LMI >= 1e-6*eye(2*nx+3*ns), P >= 1e-6*eye(nx+ns), gamma >= 0, s >= 0.1, trace(P)==s];
    obj = gamma;
    options = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(constraints, obj, options);
    
    alpha = 1/value(s);
    gamma = value(gamma)*alpha;
    P = value(P)*alpha;
    Y = value(Y)*alpha;
    G = value(G)*alpha;

    R3 = inv(Baf'*Baf)*Baf'*inv(P)*Y;
    R4 = inv(Baf'*Baf)*Baf'*inv(P)*G;
end