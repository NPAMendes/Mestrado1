clear; clc;

%% Definição das Matrizes do Sistema
A = [1.2, 0.5; 
     0,   0.8];
B = [1; 0.0];
C = [1, 0; 0, 1];
D = [0; 1];

%% Verificação do Espaço Nulo do Sistema
syms z

nx = size(A,1);

M = [A - z*eye(nx), B; 
     C,             zeros(size(C,1), size(B,2))];

condicao = simplify(det(M'*M));

zeros_invariantes = solve(condicao == 0, z);

%% Obtem o valor do ganho do controlador
[K, sol1] = opt_control(A,B,1e-6);

%% Definição dos parâmetros do estimador

sigma = 0.99;
kappa = 10;
mu = 0.1;

[Lt, L, sol2] = opt_estimation(A,B,C,D,sigma,kappa,mu,eps);

%% Simulação

% Definição das dimensões do sistema
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);
ns = size(D,2);

Aa = [A, zeros(nx, ns); zeros(ns, nx), eye(ns)];
Ba = [B; zeros(ns, nu)];
Ca = [C, D];
Ma = [zeros(nx, ns); eye(ns)];
Cs = Ma';

t = 0:1:2000;

% Inicializações
Tx = length(t);

xa = zeros(nx+ns,Tx+1);
y = zeros(ny,Tx);
u = zeros(nu,Tx);
au = zeros(nu,Tx);
ay = zeros(ns,Tx+1);

ym = zeros(ny,Tx);
aum = zeros(nu,Tx);
aym = zeros(ns,Tx);
taum = zeros(nu,Tx+1);
xam = zeros(nx+ns,Tx+1);

eau = zeros(nu,Tx);
eay = zeros(ns,Tx);

% Initial Values (optional)
xa(1:(nx+ns),1) = [-0.1; 0.2; 0];
xam(1:(nx+ns),1) = [0; 0; 0];

% falhas em rampa, degrau e senoide
au(:,1:Tx) = attack_au(t,nu);
ay(:,1:Tx+1) = attack_ay(t,ns);

for k = 1:Tx
    % System equations
    if k > 2000
        u(:,k) = K*xa(1:end-1,k)-aum(:,k);
    else
        u(:,k) = K*xa(1:end-1,k);
    end

    xa(:,k+1) = Aa*xa(:,k) + Ba*u(:,k) + Ba*au(:,k) + Ma*(ay(:,k+1)-ay(:,k));
    y(:,k) = Ca*xa(:,k);

    % Estimator equations
    ym(:,k) = Ca*xam(:,k);
    aym(:,k) = Cs*xam(:,k);
    taum(:,k+1) = -mu*(Ba'*Ba)*aum(:,k) - mu*Ba'*(Aa*xam(:,k)+Ba*u(:,k));
    xam(:,k+1) = Aa*xam(:,k) + Ba*u(:,k) + Ba*aum(:,k) + L*(y(:,k) - ym(:,k));
    aum(:,k+1) = aum(:,k) + mu*Ba'*xam(:,k+1) + taum(:,k+1);

    % Error equations
    eau(:,k) = au(:,k) - aum(:,k);
    eay(:,k) = ay(:,k) - aym(:,k);

end

%% Plotting
figure(1)

subplot(2,2,1)
stairs(t,xa(1,1:end-1),LineWidth=1.5,LineStyle="-",Color='b')
hold on
stairs(t,xam(1,1:end-1),LineWidth=1.5,LineStyle="--",Color='r')
ylim([-1,1]);
% xlabel('$Sample(k)$','Interpreter','latex','FontSize',12);
ylabel('$x_1,\hat{x}_1$','Interpreter','latex','FontSize',16);
lgd1 = legend('$x_1$','$\hat{x}_1$');
set(lgd1,'Interpreter','latex','FontSize',16);
legend('boxon')
legend('Location','southeast')
title('(a) Real and estimated states.','Interpreter','latex','FontSize',14)

subplot(2,2,2)
stairs(t,xa(2,1:end-1),LineWidth=1.5,LineStyle="-",Color='b')
hold on
stairs(t,xam(2,1:end-1),LineWidth=1.5,LineStyle="--",Color='r')
%ylim([-0.3,0.2]);
% xlabel('$Sample(k)$','Interpreter','latex','FontSize',12);
ylabel('$x_2,\hat{x}_2$','Interpreter','latex','FontSize',16);
lgd2 = legend('$x_2$','$\hat{x}_2$');
set(lgd2,'Interpreter','latex','FontSize',16);
legend('boxon')
legend('Location','southeast')
title('(b) Real and estimated states.','Interpreter','latex','FontSize',14)

subplot(2,2,3);
stairs(t,au,LineWidth=1.5,LineStyle="-",Color='b')
hold on
stairs(t,aum(1,1:end-1),LineWidth=1.5,LineStyle="--",Color='r')
%ylim([-1.5 3]);
% xlabel('$Sample(k)$','Interpreter','latex','FontSize',12);
ylabel('$a_{u}, \hat{a}_{u}$','Interpreter','latex','FontSize',16);
lgd3 = legend('$a_{u}$','$\hat{a}_{u}$');
set(lgd3,'Interpreter','latex','FontSize',16);
legend('boxon')
legend('Location','northeast')
title('(c) Actuator attack and its estimation.','Interpreter','latex','FontSize',14)

subplot(2,2,4);

stairs(t,ay(1,1:end-1),LineWidth=1.5,LineStyle="-",Color='b')
hold on
stairs(t,aym,LineWidth=1.5,LineStyle="--",Color='r')
%ylim([-2 2.5]);
ylabel('$a_{y}, \hat{a}_{y}$','Interpreter','latex','FontSize',16);
lgd4 = legend('$a_{y}$','$\hat{a}_{y}$');
set(lgd4,'Interpreter','latex','FontSize',16);
legend('boxon')
legend('Location','southeast')
title('(d) Sensor attack and its estimation.','Interpreter','latex','FontSize',14)