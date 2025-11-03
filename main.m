clear; clc; close('all');

% Definição do sistema
A = [ 1.178   0.001   0.511  -0.403;
     -0.051   0.661  -0.011   0.061;
      0.076   0.335   0.560   0.382;
      0.000   0.335   0.089   0.849 ];

B = [ 0.004  -0.087;
      0.467   0.001;
      0.213  -0.235;
      0.213  -0.016 ];

C = eye(4);

D = [1;
     0;
     0;
     0];

% Definição do controlador (ralimentação de saída)
[K, solC, alphaC] = Control(A, B);

if solC.problem ~= 0 || alphaC < 0
    disp('Não foi possível projetar o controlador.')
else 
    writematrix(K, 'K.csv');
end

% Definição do estimador 
[L, solL, gammaL, alphaL] = Estimador(A, C, D);

if solL.problem ~= 0
    disp('Não foi possível projetar o estimador.')
else 
    writematrix(L, 'L.csv');
end

% Definição do bloco de reconfiguração
sigma = [0.12, 0.13];

[R1, R2, R3, R4, solR, gammaR, alphaR] = ReconfigurationBlock(A, B, C, D, K, L, sigma);

if solR.problem ~= 0 || alphaR<0 || gammaR<0
    disp('Não foi possível projetar o bloco de reconfiguração.')
else
    writematrix(R1, 'R1.csv');
    writematrix(R2, 'R2.csv');
    writematrix(R3, 'R3.csv');
    writematrix(R4, 'R4.csv');
end

Bf = B*diag(1-sigma);

nx = size(A,1);
ns = size(D,2);
nu = size(B,2);

L1 = L(1:nx, :);
    L2 = L(nx+1:end, :);

    Ac = [A + B*K - L1*C, -L1*D;
          -L2*C,          eye(ns) - L2*D];

    Bc = L;

    Cc = [K, zeros(nu, ns)];

Acl = [A+Bf*R3*C, Bf*R4*Cc;
       Bc*R1*C,   Ac+Bc*R2*Cc];

if all(abs(eig(Acl))<1)
    disp('Sistema estável')
end
