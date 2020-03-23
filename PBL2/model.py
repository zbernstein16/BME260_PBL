import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


Volume = 5  # Liters. Total Blood volume

################################ Properties ################################

##################### T Cells #################

#%%
########## Overall System

## Healthy
d_T = 0.01 # Death Rate [1/day]
k_T = 2.4e-10 # T infection Rate [1/virion-day]

## Infected
d_Tx = 1 # Death Rate [1/day]

#%%
##########  Blood

## Healthy
T_B_0 = 1e6 # Initial blood T cell count [cell/mL]
lambda_T_B = d_T*T_B_0 # [cell/mL/day] Gen. of blood T cell 
#                                     (from Steady state of dT_B/dt = lambda_T_B - d_T * T_B_0 = 0)
sig_BL = 0.01 # Blood to lymph T cell transfer [1/day]

## Infected
T_Bx_0 = 0
rho_T_Bx = 2.0e-3 # Blood T cell virus production [virion/cell/day]

#%%
########## Lymph

## Healthy
lambda_T_L = 50*lambda_T_B # Gen. of lymph T cell [cell/mL/day]
T_L_0 = lambda_T_L/d_T # Initial lymph T cell count  [cell/mL]
sig_LB = .0002

## Infected
T_Lx_0 = 0
rho_T_Lx = 2000 # Blood T cell virus production [virion/cell/day]


########## END T CELL PROPERTIES
#%%
##################### Monocytes #################

## Healthy
M_B_0 = 1.8e9 # Initial blood Monocyte cell count [count]
d_M = 0.1  # Monocyte death rate [1/day]
lambda_M = d_M * M_B_0  # Steady state of dM/dt = lambda_M - d_M * M_B --> 0 = lambda_M - d_M * M_B_0
k_M = 1.19  # 1/day
rho_M = 1  # virions/cell Viruses required to infect monocyte

## Infected
M_Bx_0 = 0
d_Mx = 0.087  # 1/day Infected monocyte death rate

#%%
##################### Virus #################


##########  Overall
d_V = 3  # 1/day Virus natural death rate
rho_V = 34  # viruses produced per infected cell
c_1 = 3.083e9  # michaelis menton virions half saturation

##########  Blood
V_B_0 = 1

##########  Lymph

V_L_0 = 1





#%%
################################ Model ################################

def model(X,t):

    ## Get current state values

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

    ## Update States

    # Blood Variables
            #SS Gen     SS death   B->L trans L->B trans  - Infection
    dT_B = lambda_T_B - d_T*T_B - sig_BL*T_B + sig_LB*T_L - k_T*V_B*T_B
            #Infection     Inf Death
    dT_Bx = k_T*V_B*T_B - d_Tx*T_Bx
            
            #TODO: Add monocytes
    dM_B = 0
    dM_Bx =  0
    
            #V from inf T  nat V Death
    dV_B = rho_T_Bx*T_Bx - d_V*V_B
    
    # Lymph Variables
            #SS Gen      SS death  B->L trans   L->B trans      - Infection
    dT_L =  lambda_T_L - d_T*T_L + sig_BL*T_B - sig_LB*T_L - k_T*V_L*T_L
            #Infection - Inf Death
    dT_Lx = k_T*V_L*T_L - d_Tx*T_Lx
            #V from inf T  nat V Death 
    dV_L = rho_T_Lx*T_Lx - d_V*V_L

    return [dT_B,
    dT_Bx,
    dM_B,
    dM_Bx,
    dV_B,
    dT_L,
    dT_Lx,
    dV_L,]

#%%
################################ Simulation ################################
x0 = [T_B_0,T_Bx_0,M_B_0,M_Bx_0,V_B_0,T_L_0,T_Lx_0,V_L_0]
t = np.linspace(0,500,5e3)
sol,info = odeint(model, x0, t,full_output=True)

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

#%%
################################ Plotting ################################
plt.figure(1)

plt.subplot(1,3,1)
plt.plot(t, V_B,'g-')
plt.ylabel('Virions')
plt.xlabel("Days")

plt.subplot(1,3,2)
plt.plot(t,T_B,'k-')
plt.plot(t,T_Bx,'r-')
plt.legend(['Healthy T','Infected T'])


plt.title("BLOOD")

plt.show()


plt.figure(2)

plt.subplot(1,3,1)
plt.plot(t, V_L,'g-')
plt.ylabel('Virions')
plt.xlabel("Days")

plt.subplot(1,3,2)
plt.plot(t,T_L,'k-')
plt.plot(t,T_Lx,'r-')
plt.legend(['Healthy T','Infected T'])

plt.title("LYMPH")

plt.show()

