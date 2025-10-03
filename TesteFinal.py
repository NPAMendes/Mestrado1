import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

A = np.array([[1.2, 0.5],
              [0,   0.8]])

B = np.array([[1],
              [0]])

C = np.array([[1, 0]])

D = np.array([[0]])

n = 2 # Estados
m = 1 # Entradas
p = 1 # Saída

N = 200
k = np.arange(0,N,1)

X = np.zeros((n, N))
X.T[0] = [1,1]

def system(x, u):
    return (A @ x + (B.flatten() * u))

for i in range (N-1):
    X.T[i+1] = system(X.T[i], 0)
    
plt.figure()
plt.subplot(211)
plt.plot(k, X[0])

plt.subplot(212)
plt.plot(k, X[1])

#%%

Ky = np.array([[-1.2]])

N = 5000
k = np.arange(0,N,1)

X = np.zeros((n, N))
X.T[0] = [1,1]

for i in range (N-1):
    X.T[i+1] = system(X.T[i], Ky@C@X.T[i].reshape(-1,1))
    
plt.figure()
plt.subplot(211)
plt.plot(k, X[0])

plt.subplot(212)
plt.plot(k, X[1])

u = Ky@(C@X)

plt.figure()
plt.plot(k, u[0])

#%%
mu = 1

Da = np.eye(p)

L1 = np.array([[-29.3273],
              [208591778.1075]])
L2 = np.array([[271.5512]])
                  
Ac = np.block([[A-L1@C, -(L1+B@Ky)@Da,  np.zeros((n,m))],
               [-L2@C,  np.eye(p)-L2@D, np.zeros((p,m))],
               [-mu*B.T@L1@C, -mu*B.T@L1@Da, np.eye(m)]])

Bc = np.block([[B@Ky + L1],
               [L2],
               [mu*B.T@L1]])

Cc = np.block([[np.zeros((m,n)), -Ky@Da, -np.eye(m)]])

Dc = Ky

def control (y, xc):
    xcm = Ac@xc + Bc@y
    u = Cc@xc + Dc@y
    return (xcm, u)

#%%
plt.close('all')
N = 200
k = np.arange(0,N,1)

X = np.zeros((n, N))
X.T[0] = [0.1,0.1]

Xc = np.zeros((n+p+m, N))
Xc.T[0] = [0.1,0.1,0,0]

u = np.zeros(N)

for i in range (N-1):
    
    Xc.T[i+1], u[i] = control(C@X.T[i], Xc.T[i])
        
    X.T[i+1] = system(X.T[i], u[i])
    
plt.figure()
plt.subplot(211)
plt.plot(k, X[0])

plt.subplot(212)
plt.plot(k, X[1])

plt.figure()
plt.plot(k,u)

plt.figure()
plt.subplot(211)
plt.plot(k, Xc[2])

plt.subplot(212)
plt.plot(k, Xc[3])

#%%
plt.close('all')

N = 320
k = np.arange(0,N,1)

X = np.zeros((n, N))
X.T[0] = [1,1]

Xc = np.zeros((n+p+m, N))
Xc.T[0] = [1,1,0,0]

u = np.zeros(N)

au = np.zeros(N)
ay = np.zeros(N)

for i in range (N-1):
    
    if i >=150 and i <550: ay[i] = 1.5
    elif i >= 650 and i < 1300: ay[i] = np.sin(np.pi*(i-1300)/1300)
    elif i >= 1400 and i < 1900: ay[i] = i/380 - 3
    
    Xc.T[i+1], u[i] = control(C@X.T[i] + ay[i], Xc.T[i])
    
    if i >=150 and i <500: au[i] = 4*i/1000
    elif i >= 750 and i < 1200: au[i] = 0.75
    elif i >= 1350 and i < 2000: au[i] = np.sin(5*np.pi*(i-1350)/1000)
        
    X.T[i+1] = system(X.T[i], u[i]+au[i])
    
plt.figure()
plt.subplot(211)
plt.plot(k, X[0])

plt.subplot(212)
plt.plot(k, X[1])

plt.figure()
plt.plot(k,u)