import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

lam_1 = 1e4
k = 2.4e-8
d_1 = .01
sig_1 = .01
gam_r = 5e-6
sig_2 = .0002
lam_2 = 5e5
d_2 = 1
f = .95
d_3 = .1 # only this number was changed
p_v1 = 1000
d_4 = 23
D_1 = .1
p_v2 = 2000
D_2 = .2
N_C = 15
d_5 = 6.6


T_B_0 = 1e6
T_L_0 = 5e7
T_Bx_0 = 0
T_Lx_0 = 0
Mx_0 = 0
V1_0 = 1e-3
V2_0 = 1e-3
C_0 = 0



def model(X,t):
    T_B = X[0]
    T_L = X[1]
    T_Bx = X[2]
    T_Lx = X[3]
    Mx = X[4]
    V1 = X[5]
    V2 = X[6]
    C = X[7]

    dT_B = lam_1 - k*V1*T_B -d_1*T_B - sig_1*(1+gam_r*C)*T_B + sig_2*T_L
    dT_L = lam_2 + sig_1*(1+gam_r*C)*T_B - k*V2*T_L - sig_2*T_L - d_1*T_L
    dT_Bx = k*V1*T_B - d_2*T_Bx
    dT_Lx = (1-f)*k*V2*T_L - d_2*T_Lx
    dMx = f*k*V2*T_L - d_3*Mx
    dV1 = p_v1*T_Bx - d_4*V1 + D_1*(V2 - V1)
    dV2 = p_v2*T_Lx - d_4*V2 + D_2*(V1 - V2)
    dC = N_C*d_3*Mx - d_5*C

    return [dT_B,dT_L,dT_Bx,dT_Lx,dMx,dV1,dV2,dC]

### Simulation
x0 = [
T_B_0,
T_L_0,
T_Bx_0,
T_Lx_0,
Mx_0,
V1_0,
V2_0,
C_0]

t = np.linspace(0,365*10,1e6)
sol,info = odeint(model, x0, t,full_output=True)


T_B = sol[:,0]
T_L = sol[:,1]
T_Bx = sol[:,2]
T_Lx = sol[:,3]
Mx = sol[:,4]
V1 = sol[:,5]
V2 = sol[:,6]
C = sol[:,7]

plt.figure(1)

plt.plot(t,T_B,'k-')
plt.plot(t,T_Bx,'r-')



plt.show()

plt.figure(1)

plt.subplot(1,2,1)
plt.plot(t,T_L,'k-')
plt.plot(t,T_Lx,'r-')
plt.plot(t,Mx,'b-')

plt.subplot(1,2,2)
plt.plot(t,V2)


plt.show()
plt.figure(3)

plt.plot(t,C)

plt.show()

