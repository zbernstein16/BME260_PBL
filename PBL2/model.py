import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

Volume = 5  # Liters


####### Blood #######

#Healthy T Cell properties
T_B_0 = 5e9 # Initial blood T cell count [count]
lambda_T_B = 5e7 # Generation of CDT cell in blood
sig_BL = 0.01 # Blood to lymph T cell transfer [1/day]
d_T = 0.01
# Infected T Cell properties
T_Bx_0 = 0

# Healthy Monocyte properties
M_B_0 = 1.8e9 # Initial blood Monocyte cell count [count]
d_M = 0.1  # Monocyte death rate [1/day]
lambda_M = d_M * M_B_0  # Steady state of dM/dt = lambda_M - d_M * M_B --> 0 = lambda_M - d_M * M_B_0
k_M = 1.19  # 1/day
rho_M = 1  # virions/cell Viruses required to infect monocyte

# Infected monocyte properties
M_Bx_0 = 0
d_Mx = 0.087  # 1/day Infected monocyte death rate

# Virus properties
V_B_0 = 1e12
d_V = 3  # 1/day Virus natural death rate
rho_V = 34  # viruses produced per infected cell
c_1 = 3.083e9  # michaelsen menton virions half saturation

####### Lymph #######

#Healthy T Cell properties
lambda_T_L = 50*lambda_T_B
T_L_0 = lambda_T_L/d_T
sig_LB = .0002

# Infected T Cell properties
T_Lx_0 = 0

# Virus properties
V_L_0 = 0

def model(X, t):

    # Blood Variables
    T_B = X[0] # Healthy Blood CD4 T
    T_Bx = X[1] # Infected Blood CD4 T
    M_B = X[2] # Healthy BLood monocyte
    M_Bx = X[3] # Infected blood monocyte
    V_B = X[4] # Blood virion
    # Lymph Variables
    T_L = X[5] # Healthy lymph CD4 T
    T_Lx = X[6] # Infected Lymph CD4 T
    V_L = X[7] # Lymph Virion

    # Blood Variables
    dT_B = lambda_T_B - d_T*T_B - sig_BL*T_B + sig_LB*T_L
    dT_Bx =  0
    dM_B = lambda_M - d_M * M_B
    dM_Bx =  0
    dV_B =   0
    # Lymph Variables
    dT_L =  lambda_T_L - d_T*T_L + sig_BL*T_B - sig_LB*T_L
    dT_Lx = 0
    dV_L = 0

    #
    # dV_B = rho_V * M_Bx - d_V * V_B
    # # produced - death - used to infect
    # dM_B = lambda_M - d_M * M_B - (2.4e-12) * V_B * M_B
    # dM_Bx = (2.4e-12) * V_B * M_B - d_Mx * M_Bx

    return [dT_B,
    dT_Bx,
    dM_B,
    dM_Bx,
    dV_B,
    dT_L,
    dT_Lx,
    dV_L,]


x0 = [T_B_0,T_Bx_0,M_B_0,M_Bx_0,V_B_0,T_L_0,T_Lx_0,V_L_0]
t = np.linspace(0, 80, 1e5)
sol = odeint(model, x0, t)

# Blood Variables
T_B = sol[:,0]  # Healthy Blood CD4 T
T_Bx = sol[:,1]  # Infected Blood CD4 T
M_B = sol[:,2]  # Healthy BLood monocyte
M_Bx = sol[:,3] # Infected blood monocyte
V_B = sol[:,4]  # Blood virion
# Lymph Variables
T_L = sol[:,5]  # Healthy lymph CD4 T
T_Lx = sol[:,6]  # Infected Lymph CD4 T
V_L = sol[:,7]  # Lymph Virion


plt.figure(1)

plt.subplot(1,3,1)
plt.plot(t, V_B)
plt.ylabel('Virions')


plt.subplot(1,3,2)
plt.plot(t,T_B,'k-')
plt.plot(t,T_Bx,'r-')
plt.legend(['Healthy T','Infected T'])


plt.show()
plt.figure(2)
plt.plot(t, M_B, 'k-')
plt.plot(t, M_Bx, 'r-')
plt.ylabel('Monocyte')
plt.legend(['Healthy', 'Infected'])
plt.show()

plt.figure(3)

plt.show()
