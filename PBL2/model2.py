import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


## Sources
# A: Modeling of HIV-1 Infection: Insights to the Role of
# Monocytes/Macrophages, Latently Infected T4 Cells, and
# HAART Regimes
## https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0046026&type=printable

# B: Modeling the Slow CD4+ T Cell Decline in HIV-Infected Individuals
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4692447/#pcbi.1004665.ref031

# C :
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5502436/

# Physiological Properties
Vol_CNS = 1273.6 + 150 # mL
Vol_blood = 5000 # mL

# Initial Conditions
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
M_CNSx_0 = 0
G_CNS_0 = .05*(1.7e9+100e9)/Vol_CNS # [cell/mL]
G_CNSx_0 = 0
V_CNS_0 = 0


############# Overall
k_T = 8e-9# Fit, OLD: 2.4e-8 # T cell Infection rate (B)
k_M = k_T*(1.19/.089) # Monocyte Infection rate (A)
k_G = k_M # unsure of this value! will study its variations
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
p_MB = 34 # Monocyte Virus production rate [particles/cell/day]
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
    M_CNSx = X[10]
    G_CNS = X[11]
    G_CNSx = X[12]
    V_CNS = X[13]


    
    k_Tenh = k_T * (1+t/2000) # ENHANCED
    #      CD4 GEN   Virus->Inf  Nat Death      Pyroptosis Attr    Lymph->Blood
    dT_B = lam_T_B - k_Tenh*V1*T_B -d_T*T_B - sig_1*(1+cytoGamma(C)*C)*T_B + sig_2*T_L
   
    #       Virus->Inf   DeathInf (100x natural)  
    dT_Bx = k_Tenh*V1*T_B - d_Tprodx*T_Bx
    
    #Blood
    #   Viral Prod  T. Prod. M   Viral Death    Transfer Rate
    dV1 = p_v1 * T_Bx + p_MB*M_Bx - d_V * V1 + D_1 * (V2 - V1)
    
    #     Mono GEN   Nat Death   Bl -> CNS     Virus-> Inf
    dM_B = lam_M_B - d_M*M_B - D_MB2CNS*M_B - k_M*V1*M_B
    
    #        Virus->Inf    Inf Death    Leaving
    dM_Bx = k_M*V1*M_B - d_Mx*M_Bx - D_MB2CNS*M_Bx
    
    #      CD4 GEN     Pyroptosis Attraction            Virus->Inf    Lymph->Blood   Nat Death
    dT_L = lam_T_L + sig_1 * (1 + cytoGamma(C) * C) * T_B - k_Tenh * V2 * T_L - sig_2 * T_L - d_T * T_L
    
    #               prod. Gen           prod. death
    dT_Lprodx = (1-f)*k_Tenh*V2*T_L - d_Tprodx*T_Lprodx
    
    #               abort Gen          abort death
    dT_Labortx = f*k_Tenh*V2*T_L - d_Tabortx*T_Labortx
    
    #     Prod.Cell     Vir death   virus Transfer?
    dV2 = p_v2*T_Lprodx - d_V*V2 + D_2*(V1 - V2)
    
    #     #C * death * #abortive    Cyto Death             
    dC = N_C*d_Tabortx*T_Labortx - d_C*C


    # CNS
    dM_CNSx = D_MB2CNS*M_Bx*Vol_blood/Vol_CNS - d_Mx*M_CNSx
    dG_CNS = - k_G*G_CNS*V_CNS
    dG_CNSx = k_G*G_CNS*V_CNS - d_Mx * G_CNSx
    dV_CNS = p_MB*M_CNSx - d_V*V_CNS

    return [dT_B,dT_L,dT_Bx,dT_Lprodx,dT_Labortx,dV1,dV2,dC,dM_B,dM_Bx,dM_CNSx,dG_CNS,dG_CNSx,dV_CNS]


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
M_Bx_0,
M_CNSx_0,
G_CNS_0,
G_CNSx_0,
V_CNS_0]



t = np.linspace(0,365*10,1e5)
aids = []
for ts in t:
    aids.append(200000)
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
M_CNSx = sol[:,10]
G_CNS = sol[:,11]
G_CNSx = sol[:,12]
V_CNS = sol[:,13]

axes_size = "15"
plt.figure(0).clf()
fig0, ax0 = plt.subplots(num=0,clear=True)
ax0.plot(t/365, T_B, 'k-', label = "Healthy")
ax0.plot(t/365, aids, 'r--', label = "AIDS Threshold")
ax0.set_xlabel("Time (Years)", size = axes_size)
ax0.set_ylabel("Concentration (Cells/mL)", size = axes_size)
ax0.set_title("CD4 Concentration in the Blood")
plt.xticks(np.arange(0,11,1))
ax0.legend(loc = "upper right")
fig0.tight_layout()
#fig0.savefig("CD4_Blood.png")



plt.figure(1).clf()
fig1, ax1 = plt.subplots(num=1,clear=True)
ax1.plot(t/365,T_L,'k-')
ax1.plot(t/365,T_Lprodx,'r-')
ax1.plot(t/365,T_Labortx,'b-')
ax1.set_xlabel("Time (Years)", size = axes_size)
ax1.set_ylabel("Concentration (Cells/mL)", size = axes_size)
ax1.set_title("CD4 Concentration in Lymph")
ax1.legend(['Healthy','Productively Inf.','Abortively Inf.'], loc = "upper center", prop={"size":15})
fig1.tight_layout()
#fig1.savefig("CD4_Lymph.png")




plt.figure(2).clf()
fig2, ax2 = plt.subplots(num=2,clear=True)
ax2.plot(t/365,M_B,'k-')
ax2.plot(t/365,M_Bx,'r-')
ax2.set_xlabel("Time (Years)", size = axes_size)
ax2.set_ylabel("Concentration (Cells/mL)", size = axes_size)
ax2.set_title("Monocyte Concentration in Blood")
ax2.legend(['Healthy','Infected'], loc = "lower right", prop={"size":15})
plt.yscale('log')
fig2.tight_layout()
#fig2.savefig("Monocyte_Blood.png")





# Virus


plt.figure(3).clf()
fig3, ax3 = plt.subplots(num=3,clear=True)
ax3.plot(t/365,V1,'k-')
ax3.plot(t/365,V2,'r-')
ax3.plot(t/365,V_CNS,'b-')
ax3.set_xlabel("Time (Years)", size = axes_size)
ax3.set_ylabel("Concentration (Cells/mL)", size = axes_size)
ax3.set_title("Virus Concentrations")
plt.yscale("log")
ax3.legend(["Blood","Lymph","Brain"], prop={"size":15})
fig3.tight_layout()
#fig3.savefig("Virions.png")





# Monocytes in Brain


plt.figure(4).clf()
fig4, ax4 = plt.subplots(num=4,clear=True)
ax4.plot(t/365,G_CNS,'k-')
#ax4.plot(t/365,M_CNSx,'r-')
ax4.plot(t/365,G_CNSx,'r-')
ax4.set_xlabel("Time (Years)", size = axes_size)
ax4.set_ylabel("Concentration (Cells/mL)", size = axes_size)
ax4.set_title("HIV Infection in the Brain")
#plt.yscale("log")
plt.legend(["Healthy Microglia", "Infected Microglia"], prop={"size":15})
fig4.tight_layout()
#fig4.savefig("Monocyte_Brain.png")




###### Parameter variation Tests
# Some colors to use
Ncolors = 5

rmap = cm.get_cmap('winter', Ncolors)
reds = rmap(np.linspace(0, 1, Ncolors))

bmap = cm.get_cmap('Blues', Ncolors)
blues = bmap(np.linspace(0, 1, Ncolors))

gmap = cm.get_cmap('Greens', Ncolors)
greens = gmap(np.linspace(0, 1, Ncolors))


## e.g. using blues
# plt.figure(100)
# x = np.linspace(1,10)
# k_test = np.linspace(10,20) # some parameter we are going to vary to plot
# norm = mpl.colors.Normalize(vmin=k_test.min(),vmax=k_test.max())
# cmap = mpl.cm.ScalarMappable(norm=norm, cmap=bmap)
# cmap.set_array([])
# fig, ax = plt.subplots(figsize=(6, 6))
# for indx,k in enumerate(np.linspace(1,5,Ncolors)):
#     y = k*(x**2)
#     ax.plot(x,y,c=blues[indx])
#
# fig.colorbar(cmap, ticks=k_test[::5])
# fig.show()

def model_varyparam(X,t,param):

    global k_G # change this to the parameter you are varying
    k_G = param # change LHS
    return model(X,t)

param_values = np.linspace(1e-9,1e-8,Ncolors)
norm = mpl.colors.Normalize(vmin=param_values.min(),vmax=param_values.max())

figA,axA=plt.subplots()
for indx,param in enumerate(param_values):
    t = np.linspace(0,365*10,1e5)
    sol,info = odeint(model_varyparam, x0, t,full_output=True,args=(param,))


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
    M_CNSx = sol[:,10]
    G_CNS = sol[:,11]
    G_CNSx = sol[:,12]
    V_CNS = sol[:,13]


#    axA.plot(t, G_CNS, c=greens[indx])
    axA.plot(t,G_CNSx,c=reds[indx])

# Creating fake legend items
red_line = mpl.lines.Line2D([0], [0], color=rmap(.6), lw=4)
blue_line =  mpl.lines.Line2D([0], [0], color=bmap(.6), lw=4)
green_line = mpl.lines.Line2D([0], [0], color=gmap(.6), lw=4)


## Setting up figure A

axA.set_yscale("log")
axA.set_ylim([1e1,1e7])
# Fake legends
axA.legend([green_line,red_line],['Healthy Microglia','Infected Microglia'])
# Set colorbar for figure A
colormap4colorbar = rmap # use red as most useful
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=colormap4colorbar)
cmap.set_array([])

param_values2label = param_values[0:Ncolors:4] # maybe too many labels, so this just chooses a selection of them
cbar = figA.colorbar(cmap,label="k_G",ticks=param_values2label)
cbar.ax.set_yticklabels(["{:.2e}".format(x) for x in param_values2label])


axA.set_title("Microglia")


## Setting up figure B




figA.show()
figA.savefig("Spread Plot")

#f= open("No_Pyroptosis.txt","w+")
#for num in range(0, 100000):
#    f.write("{} {} {} {}\n".format(T_B[num], T_L[num], V1[num], G_CNS[num]))
#f.close()



