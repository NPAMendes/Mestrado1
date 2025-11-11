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
sigma = [0, 0.8];

[R3, R4, solR, gammaR, alphaR] = ReconfigurationBlock(A, B, C, D, K, sigma);

if solR.problem ~= 0 || alphaR<0 || gammaR<0
    disp('Não foi possível projetar o bloco de reconfiguração.')
else
    writematrix(R3, 'R3.csv');
    writematrix(R4, 'R4.csv');
end