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

N = 100
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

for i in range (N-1):
    X.T[i+1] = system(X.T[i], K@X.T[i])
    
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

u = K@X

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

L1 = L[:nx, :]
L2 = L[nx:, :]

Ac = np.block([[A + B @ K - L1 @ C,  -L1 @ D],
               [-L2 @ C,             np.eye(ns) - L2 @ D]])

Bc = np.block([[L1],
               [L2]])

Cc = np.block([[K, np.zeros((nu, ns))]])

def control (y, xc):
    xcm = Ac@xc + Bc@y
    u = Cc@xc
    return (xcm, u)

N = 300
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [1,1,1,1]

Xc = np.zeros((nx+ns, N))
Xc.T[0] = [0,0,0,0,0]

u = np.zeros((nu,N))

for i in range (N-1):
    
    Xc.T[i+1], u.T[i] = control(C@X.T[i], Xc.T[i])
        
    X.T[i+1] = system(X.T[i], u.T[i])
    
u.T[-1] = u.T[-2]

plt.figure('Sistema-Controlador-Estimador')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xc[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xc[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xc[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xc[3], label = 'Estimado')
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
plt.plot(k, Xc[4])
plt.ylabel('$\\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

#%% Ataque de sensor
N = 2000
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [0,0,0,0]

Xc = np.zeros((nx+ns, N))
Xc.T[0] = [0,0,0,0,0]

u = np.zeros((nu, N))

ay = np.zeros((1,N))

for i in range (N-1):
    
    if i >=150 and i <550: ay[0][i] = 1.5
    elif i >= 650 and i < 1300: ay[0][i] = np.sin(np.pi*(i-1300)/1300)
    elif i >= 1400 and i < 1900: ay[0][i] = i/380 - 3
    
    Xc.T[i+1], u.T[i] = control(C@X.T[i] + D@ay.T[i], Xc.T[i])
    
    X.T[i+1] = system(X.T[i], u.T[i])
    
plt.figure('Estimativa de ataque no sensor')
plt.plot(k, ay[0], label = 'Real')
plt.plot(k, Xc[4], label = 'Estimado')
plt.legend(fontsize = 20)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)
plt.ylabel('$a_y, \\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)

plt.figure('Sistema-Controlador-Estimador Sob ataque')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xc[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xc[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xc[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xc[3], label = 'Estimado')
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
sigma = np.array([0.12, 0.13])
Bf = B@np.diag(1-sigma)

def systemF(x, u):
    return (A @ x + Bf@u)

N = 2300
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [0,0,0,0]

Xc = np.zeros((nx+ns, N))
Xc.T[0] = [0,0,0,0,0]

u = np.zeros((nu,N))

ay = np.zeros((1,N))

for i in range (N-1):
    
    if i >=450 and i <850: ay[0][i] = 1.5
    elif i >= 950 and i < 1600: ay[0][i] = np.sin(np.pi*(i-1300)/1300)
    elif i >= 1700 and i < 2200: ay[0][i] = i/380 - 3
    
    Xc.T[i+1], u.T[i] = control(C@X.T[i] + D@ay.T[i], Xc.T[i])
    
    if k[i]<=300:
        X.T[i+1] = system(X.T[i], u.T[i])
    else:
        X.T[i+1] = systemF(X.T[i], u.T[i])
u.T[-1] = u.T[-2]

plt.figure('Sistema-Controlador-Estimador-Falha')
plt.subplot(221)
plt.plot(k, X[0], label = 'Real')
plt.plot(k, Xc[0], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_1$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(222)
plt.plot(k, X[1], label = 'Real')
plt.plot(k, Xc[1], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_2$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()

plt.subplot(223)
plt.plot(k, X[2], label = 'Real')
plt.plot(k, Xc[2], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_3$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.subplot(224)
plt.plot(k, X[3], label = 'Real')
plt.plot(k, Xc[3], label = 'Estimado')
plt.legend(fontsize = 15)
plt.ylabel('$x_4$', fontsize = 30)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

plt.figure('Sinal de controle-Controlador-Estimador-Falha')
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

plt.figure('Ataque durante estabilização - Falha')
plt.plot(k, Xc[4])
plt.ylabel('$\\hat{a}_y$', fontsize = 20)
plt.tick_params(axis = 'both', labelsize = 20)
plt.gca().yaxis.get_offset_text().set_size(15)
plt.grid()
plt.xlabel('Amostras', fontsize = 20)

#%% Sistema com falha
plt.close('all')
file_path = os.path.join(base_path, "R1.csv")
R1 = np.loadtxt(file_path, delimiter=",")

file_path = os.path.join(base_path, "R2.csv")
R2 = np.loadtxt(file_path, delimiter=",")

file_path = os.path.join(base_path, "R3.csv")
R3 = np.loadtxt(file_path, delimiter=",")

file_path = os.path.join(base_path, "R4.csv")
R4 = np.loadtxt(file_path, delimiter=",")

def RByr (y,u):
    yr = R1@y + R2@u
    return yr

def RBur (y,u):
    ur = R3@y + R4@u
    return ur

N = 2300
k = np.arange(0,N,1)

X = np.zeros((nx, N))
X.T[0] = [0,0,0,0]

Xc = np.zeros((nx+ns, N))
Xc.T[0] = [0,0,0,0,0]

u = np.zeros((nu,N))

y = np.zeros((ny, N))

ay = np.zeros((1,N))

for index in range(len(k)-1):
    if index >=450 and index <850: ay[0][index] = 1.5
    elif index >= 950 and index < 1600: ay[0][index] = np.sin(np.pi*(index-1300)/1300)
    elif index >= 1700 and index < 2200: ay[0][index] = index/380 - 3
    
    if k[index] < 300:
        y.T[index] = C @ X.T[index] + D@ay.T[index]
        
        Xc.T[index+1], u.T[index] = control(y.T[index], Xc.T[index])
            
        X.T[index+1] = system(X.T[index], u.T[index])
        
    else:
        y.T[index] = C @ X.T[index] + D@ay.T[index]
        
        yr = RByr(y.T[index], Cc@Xc.T[index])
        
        Xc.T[index+1], u.T[index] = control(yr, Xc.T[index])
        
        ur = RBur(y.T[index], u.T[index])
            
        X.T[index+1] = systemF(X.T[index], ur)
        
plt.figure()
plt.subplot(221)
plt.plot(k,X[0])

plt.subplot(222)
plt.plot(k,X[1])

plt.subplot(223)
plt.plot(k,X[2])

plt.subplot(224)
plt.plot(k,X[3])

plt.figure()
plt.plot(k,Xc[4])
plt.plot(k, ay[0])


