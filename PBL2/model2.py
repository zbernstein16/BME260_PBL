import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


## Sources
# A: Modeling of HIV-1 Infection: Insights to the Role of
# Monocytes/Macrophages, Latently Infected T4 Cells, and
# HAART Regimes
## https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0046026&type=printable

# B: Modeling the Slow CD4+ T Cell Decline in HIV-Infected Individuals
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4692447/#pcbi.1004665.ref031

# C :
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5502436/

# Initial Conditions
#TODO: Changing these breaks model, but is technically correct... unclear
T_B_0 = 1e6 #  cells/mL (B)
T_L_0 = 5e7 #  cells/mL (B)
T_Bx_0 = 0
T_Lprodx_0 = 0
T_Labortx_0 = 0
V1_0 = 1e-3
V2_0 = 1e-6
C_0 = 0
M_B_0 = 1e5 # Monocytes in Blood [cell/mL] (A)
M_Bx_0 = 0

############# Overall
k_T = 8e-9# Fit, OLD: 2.4e-8 # T cell Infection rate (B)
k_M = k_T*(1.19/.089) # Monocyte Infection rate (A)
d_T = 0.01 # Death rate of healthy T cell (B). Model EXTREMELY sensitive to this
d_Tprodx = 1 # Death rate of productively infected T Cells (B)
d_M = np.log(2)/1 # Death rate of healthy monocytes [day^-1] (C, gives half-life ~1 day)
d_Mx = 0.087 # (A)
d_C = 6.6 # Decay rate of cytok_Tine (B)
d_V = 23 # Death rate of virus (B)

############# Blood

# T Cell
lam_T_B = 1e4

# Monocyte
D_MB2CNS = 0.1
lam_M_B = d_M*M_B_0 + D_MB2CNS*M_B_0

############# Lymph

# T Cell
lam_T_L = 5e5 # ORIG 5e5, 50x blood
d_Tabortx = 0.01 #ORIG 0.001 Death rate of latently infected T cells (B)
                 #CHANGED to 0.01. Should be same as death of health or greater

############# CNS

alpha_abort2prod = .0001 # Conversion of latently infected to productively
alpha_prod2abort = alpha_abort2prod/100 # Conversion of productively infected to latent



sig_1 = .01
gam_r = 5e-6
sig_2 = .0002


f = .95
p_v1 = 1000
D_1 = .1
p_v2 = 2000
D_2 = .2
N_C = 15

def cytoGamma(C):
    #All from Source (B)
    threshhold = 4000#molecules/mL
    if C < threshhold:
        return 0    
    return gam_r #mL/molecule


def model(X,t):
    T_B = X[0] #Healthy Blood
    T_L = X[1] #Healthy Lymph
    T_Bx = X[2] #Inf Blood
    T_Lprodx = X[3] #Inf Lymph Productive
    T_Labortx = X[4] #Inf Lymph Abortive (pyroptosis)
    V1 = X[5] #Virus in blood
    V2 = X[6] #Virus in lymph
    C = X[7]  #Cytokine
    M_B = X[8] #Monocytes in Blood
    M_Bx = X[9] #Inf blood monocytes
    
    k_Tenh = k_T * (1+t/2000) # ENHANCED
    #      CD4 GEN   Virus->Inf  Nat Death      Pyroptosis Attr    Lymph->Blood
    dT_B = lam_T_B - k_Tenh*V1*T_B -d_T*T_B - sig_1*(1+cytoGamma(C)*C)*T_B + sig_2*T_L
   
    #       Virus->Inf   DeathInf (100x natural)  
    dT_Bx = k_Tenh*V1*T_B - d_Tprodx*T_Bx
    
    #Blood
    #   Viral Prod.    Viral Death    Transfer Rate
    dV1 = p_v1 * T_Bx - d_V * V1 + D_1 * (V2 - V1)
    
    #     Mono GEN   Nat Death   Bl -> CNS     Virus-> Inf
    dM_B = lam_M_B - d_M*M_B - D_MB2CNS*M_B - k_M*V1*M_B
    
    #        Virus->Inf    Inf Death
    dM_Bx = k_M*V1*M_B - d_Mx*M_Bx
    
    #      CD4 GEN     Pyroptosis Attraction            Virus->Inf    Lymph->Blood   Nat Death
    dT_L = lam_T_L + sig_1 * (1 + cytoGamma(C) * C) * T_B - k_Tenh * V2 * T_L - sig_2 * T_L - d_T * T_L
    
    #               prod. Gen           prod. death
    dT_Lprodx = (1-f)*k_Tenh*V2*T_L - d_Tprodx*T_Lprodx
    
    #               abort Gen          abort death
    dT_Labortx = f*k_Tenh*V2*T_L - d_Tabortx*T_Labortx
    
    #Notice that prod. + abort GEN = total Inf GEN
    #TODO Try estimating productive and abortive as constant ratio
    
    #     Prod.Cell     Vir death   virus Transfer?
    dV2 = p_v2*T_Lprodx - d_V*V2 + D_2*(V1 - V2)
    
    #     #C * death * #abortive    Cyto Death             
    dC = N_C*d_Tabortx*T_Labortx - d_C*C

    return [dT_B,dT_L,dT_Bx,dT_Lprodx,dT_Labortx,dV1,dV2,dC,dM_B,dM_Bx]


### Simulation
x0 = [
T_B_0,
T_L_0,
T_Bx_0,
T_Lprodx_0,
T_Labortx_0,
V1_0,
V2_0,
C_0,
M_B_0,
M_Bx_0]

t = np.linspace(0,365*10,1e5)
sol,info = odeint(model, x0, t,full_output=True)


T_B = sol[:,0]
T_L = sol[:,1]
T_Bx = sol[:,2]
T_Lprodx = sol[:,3]
T_Labortx = sol[:,4]
V1 = sol[:,5]
V2 = sol[:,6]
C = sol[:,7]
M_B = sol[:,8]
M_Bx = sol[:,9]


plt.figure(1)
plt.title("Blood T Cells")
plt.plot(t,T_B,'k-')
#plt.plot(t,T_Bx,'r--')

plt.legend(['Healthy'])
plt.yscale("log")
plt.show()

plt.figure(2)
plt.title("Lymph T Cells")
plt.plot(t,T_L,'k-')
plt.plot(t,T_Lprodx,'r-')
plt.plot(t,T_Labortx,'b-')

plt.legend(['T Lymph','T product','T abortive'])
plt.show()

# Blood monocytes
#plt.figure(3)
#plt.title("Blood Monocyte Profile")
#
#plt.plot(t,M_B,'k-')
#plt.plot(t,M_Bx,'r-')
#
#plt.legend(['Healthy','Infected'])
#plt.show()

# Virus
plt.figure(4)
plt.title("Virions")
plt.plot(t,V1,'k-')
plt.plot(t,V2,'b-')

plt.legend(["Blood","Lymph"])
plt.yscale("log")
plt.show()