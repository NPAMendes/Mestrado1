import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all')

# Definição do sistema
A = np.array([[ 1.178,  0.001,  0.511, -0.403],
              [-0.051,  0.661, -0.011,  0.061],
              [ 0.076,  0.335,  0.560,  0.382],
              [ 0.000,  0.335,  0.089,  0.849]])

B = np.array([[ 0.004, -0.087],
              [ 0.467,  0.001],
              [ 0.213, -0.235],
              [ 0.213, -0.016]])

C = np.eye(4)

D = np.array([[1],
              [0],
              [0],
              [0]])
#%% Simulação malha aberta
nx = A.shape[0] # Estados
nu = B.shape[1] # Entradas
ny = C.shape[0] # Saída
ns = D.shape[1] # Ataque de sensor

N = 200
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

def system(x, u):
    return (A @ x + B@u)

for i in range (N-1):
    X.T[i+1] = system(X.T[i], np.zeros(B.shape[1]))
    
plt.figure('Malha Aberta')
plt.subplot(221)
plt.plot(k, X[0])
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1])
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2])
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 30)

plt.subplot(224)
plt.plot(k, X[3])
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 30)

#%% Simualçaõ com controlador u = Kx
base_path = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(base_path, "K.csv")

K = np.loadtxt(file_path, delimiter=",")

def control(x):
    return K@x

N = 100
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

u = np.zeros((nu,N))

for i in range (N-1):
    u.T[i] = control(X.T[i])
    X.T[i+1] = system(X.T[i], u.T[i])
    
plt.figure('Sistema Controlador')
plt.subplot(221)
plt.plot(k, X[0])
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1])
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2])
plt.ylabel('$x_3$', fontsize = 30)
plt.xlabel('Amostras', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(224)
plt.plot(k, X[3])
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Sinal de controle do sistema controlado')
plt.subplot(121)
plt.plot(k, u[0])
plt.ylabel('$u_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(122)
plt.plot(k, u[1])
plt.ylabel('$u_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

#%% Simulações com controlador e observador

file_path = os.path.join(base_path, "L.csv")

L = np.loadtxt(file_path, delimiter=",")

Aa = np.block([[A,               np.zeros((nx,ns))],
               [np.zeros((ns,nx)), np.eye(ns)]])

Ba = np.block([[B],
               [np.zeros((ns,nu))]])

Ca = np.block([[C, D]])

def observer(Xo,Y,Yo, u):
    Xom = Aa@Xo + Ba@u +L@(Y - Yo)
    return Xom

N = 300
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

Xo = np.zeros((nx+ns, N))
Xo.T[0] = [0,0,0,0,0]

y = np.zeros((ny,N))
y.T[0] = [1,1,1,1]

yo = np.zeros((ny,N))
yo.T[0] = [0,0,0,0]

u = np.zeros((nu,N))

for i in range (N-1):
    
    u.T[i] = control(Xo.T[i, :nx])
    
    X.T[i+1] = system(X.T[i], u.T[i])
    y.T[i+1] = C@X.T[i+1]
    
    Xo.T[i+1] = observer(Xo.T[i], y.T[i], yo.T[i], u.T[i])
    yo.T[i+1] = Ca@Xo.T[i+1]
    
u.T[-1] = u.T[-2]

plt.figure('Sistema-Controlador-Estimador')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xo[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xo[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xo[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xo[3], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Sinal de controle-Controlador-Estimador')
plt.subplot(121)
plt.plot(k, u[0])
plt.ylabel('$u_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(122)
plt.plot(k, u[1])
plt.ylabel('$u_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Ataque durante estabilização')
plt.plot(k, Xo[4])
plt.ylabel('$\\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

#%% Ataque de sensor
N = 2000
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

Xo = np.zeros((nx+ns, N))
Xo.T[0] = [0,0,0,0,0]

y = np.zeros((ny,N))
y.T[0] = [1,1,1,1]

yo = np.zeros((ny,N))
yo.T[0] = [0,0,0,0]

u = np.zeros((nu,N))

ay = np.zeros((1,N))

for i in range (N-1):
    
    if i >=150 and i <550: ay[0][i+1] = 1.5
    elif i >= 650 and i < 1300: ay[0][i+1] = np.sin(np.pi*(i-1300)/1300)
    elif i >= 1400 and i < 1900: ay[0][i+1] = i/380 - 3
    
    u.T[i] = control(Xo.T[i, :nx])
    
    X.T[i+1] = system(X.T[i], u.T[i])
    y.T[i+1] = C@X.T[i+1] + D@ay.T[i+1]
    
    Xo.T[i+1] = observer(Xo.T[i], y.T[i], yo.T[i], u.T[i])
    yo.T[i+1] = Ca@Xo.T[i+1]
    
plt.figure('Estimativa de ataque no sensor')
plt.plot(k, ay[0], label = 'Real')
plt.plot(k, Xo[4], label = 'Estimado')
plt.legend(fontsize = 20)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)
plt.ylabel('$a_y, \\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)

plt.figure('Sistema-Controlador-Estimador Sob ataque')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xo[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xo[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xo[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xo[3], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Sinal de controle-Controlador-Estimador Sob ataque')
plt.subplot(121)
plt.plot(k, u[0])
plt.ylabel('$u_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(122)
plt.plot(k, u[1])
plt.ylabel('$u_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

#%% Sistema com falha sem bloco de reconfiguração
plt.close('all')
sigma = np.array([0, 0.8])
Bf = B@np.diag(1-sigma)

def systemF(x, u):
    return (A @ x + Bf@u)

Baf = np.block([[Bf],
                [np.zeros((ns,nu))]])

def observerF(Xo,Y,Yo, u):
    Xom = Aa@Xo + Baf@u +L@(Y - Yo)
    return Xom

N = 1000
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

Xo = np.zeros((nx+ns, N))
Xo.T[0] = [0,0,0,0,0]

y = np.zeros((ny,N))
y.T[0] = [1,1,1,1]

yo = np.zeros((ny,N))
yo.T[0] = [0,0,0,0]

u = np.zeros((nu,N))

for i in range (N-1):
    
    u.T[i] = control(Xo.T[i, :nx])
    
    if i < 300:
        X.T[i+1] = system(X.T[i], u.T[i])
        y.T[i+1] = C@X.T[i+1]
        
        Xo.T[i+1] = observer(Xo.T[i], y.T[i], yo.T[i], u.T[i])
        yo.T[i+1] = Ca@Xo.T[i+1]
        
    else:
        X.T[i+1] = systemF(X.T[i], u.T[i])
        y.T[i+1] = C@X.T[i+1]
        
        Xo.T[i+1] = observerF(Xo.T[i], y.T[i], yo.T[i], u.T[i])
        yo.T[i+1] = Ca@Xo.T[i+1]
    
u.T[-1] = u.T[-2]

plt.figure('Sistema-Controlador-Estimador-Falha')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xo[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xo[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xo[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xo[3], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Sinal de controle-Controlador-Estimador')
plt.subplot(121)
plt.plot(k, u[0])
plt.ylabel('$u_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(122)
plt.plot(k, u[1])
plt.ylabel('$u_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Ataque durante estabilização')
plt.plot(k, Xo[4])
plt.ylabel('$\\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

#%% Sistema com falha
plt.close('all')
file_path = os.path.join(base_path, "R3.csv")
R3 = np.loadtxt(file_path, delimiter=",")

file_path = os.path.join(base_path, "R4.csv")
R4 = np.loadtxt(file_path, delimiter=",")

Cs = np.block([[C, np.zeros((ny,ns))]])

def RBur (y,u):
    ur = R3@y + R4@u
    return ur

N = 1000
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

Xo = np.zeros((nx+ns, N))
Xo.T[0] = [0,0,0,0,0]

y = np.zeros((ny,N))
y.T[0] = [1,1,1,1]

yo = np.zeros((ny,N))
yo.T[0] = [0,0,0,0]

u = np.zeros((nu,N))

for i in range (N-1):
    
    u.T[i] = control(Xo.T[i, :nx])
    
    if i < 400:
        X.T[i+1] = system(X.T[i], u.T[i])
        y.T[i+1] = C@X.T[i+1]
        
        Xo.T[i+1] = observer(Xo.T[i], y.T[i], yo.T[i], u.T[i])
        yo.T[i+1] = Ca@Xo.T[i+1]
        
    else:
        u.T[i] = RBur(yo.T[i], u.T[i])
        
        X.T[i+1] = systemF(X.T[i], u.T[i])
        y.T[i+1] = C@X.T[i+1]
        
        Xo.T[i+1] = observerF(Xo.T[i], y.T[i], yo.T[i], u.T[i])
        yo.T[i+1] = Cs@Xo.T[i+1]
    
u.T[-1] = u.T[-2]

plt.figure('Sistema-Controlador-Estimador-Falha')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xo[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xo[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xo[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xo[3], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Sinal de controle-Controlador-Estimador')
plt.subplot(121)
plt.plot(k, u[0])
plt.ylabel('$u_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(122)
plt.plot(k, u[1])
plt.ylabel('$u_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Ataque durante estabilização')
plt.plot(k, Xo[4])
plt.ylabel('$\\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

#%%
N = 2000
k = np.arange(0,N,1)

X = np.zeros((nx, N))

Xo = np.zeros((nx+ns, N))

y = np.zeros((ny,N))

yo = np.zeros((ny,N))

yos = np.zeros((ny,N))

u = np.zeros((nu,N))

ay = np.zeros((1,N))

for i in range (N-1):
    
    if i >=150 and i <550: ay[0][i+1] = 1.5
    elif i >= 650 and i < 1300: ay[0][i+1] = np.sin(np.pi*(i-1300)/1300)
    elif i >= 1400 and i < 1900: ay[0][i+1] = i/380 - 3
    
    u.T[i] = control(Xo.T[i, :nx])
    
    if i < 50:
        X.T[i+1] = system(X.T[i], u.T[i])
        y.T[i+1] = C@X.T[i+1]+ D@ay.T[i+1]
        
        Xo.T[i+1] = observer(Xo.T[i], y.T[i], yo.T[i], u.T[i])
        yo.T[i+1] = Ca@Xo.T[i+1]
        
    else:
        u.T[i] = RBur(yos.T[i], u.T[i])
        
        X.T[i+1] = systemF(X.T[i], u.T[i])
        y.T[i+1] = C@X.T[i+1]+ D@ay.T[i+1]
        
        Xo.T[i+1] = observerF(Xo.T[i], y.T[i], yo.T[i], u.T[i])
        yos.T[i+1] = Cs@Xo.T[i+1]
        yo.T[i+1] = Ca@Xo.T[i+1]
    
plt.figure('Estimativa de ataque no sensor')
plt.plot(k, ay[0], label = 'Real')
plt.plot(k, Xo[4], label = 'Estimado')
plt.legend(fontsize = 20)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)
plt.ylabel('$a_y, \\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)

plt.figure('Sistema-Controlador-Estimador Sob ataque')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xo[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xo[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xo[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xo[3], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Sinal de controle-Controlador-Estimador Sob ataque')
plt.subplot(121)
plt.plot(k, u[0])
plt.ylabel('$u_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(122)
plt.plot(k, u[1])
plt.ylabel('$u_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)
