import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

A = np.array([[1.2, 0.5],
              [0,   0.8]])

B = np.array([[1],
              [0]])

C = np.array([[1, 1]])

D = np.array([[0]])

n = 2
m = 1
p = 1

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

Ky = np.array([[-0.6834]])

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

Aa = np.block([[A,           np.zeros((n,p))],
               [np.zeros((p,n)), np.eye(p)]])
Ba = np.block([[B],
               [np.zeros((p,p))]])
Ca = np.block([C, np.eye(p)])

Da = np.eye(p)

Ma = np.vstack([np.zeros((n,1)), np.eye(p)])

L = np.array([[1.3073e-10],
              [3.4473e-10],
              [3.8979e-10]])
                  
Ac = np.block([[Aa - Ba@Ky@Da@Ma.T - L@Ca, np.zeros((n+p,m))],
               [-mu*Ba.T@L@Ca,           np.eye(m)]])

Bc = np.block([[Ba@Ky + L],
               [mu*Ba.T@L]])

Cc = np.block([[-Ky@Da@Ma.T, -np.eye(m)]])

Dc = Ky

def control (y, xc):
    xcm = Ac@xc + Bc@y
    u = Cc@xc + Dc@y
    return (xcm, u)

#%%
N = 320
k = np.arange(0,N,1)

X = np.zeros((n, N))
X.T[0] = [1,1]

Xc = np.zeros((n+p+m, N))
Xc.T[0] = [1,1,0,0]

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

#%%

N = 2500
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