"""
Fichier principal d'execution du devoir 2
"""


#############################
#Importation de bibliothèques
#############################
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.stats import linregress
import config
import functions
import errors


############################
# Importation des constantes
############################
D_EFF = config.D_EFF            # en m^2/s
C_E = config.C_E                # en mol/m^3
D = config.D                    # en m
R = config.R                    # en m
Tf = config.Tf                  # Temps caractéristique
k = config.k                    # Coefficient de réaction


########################################
# Définition de la solution manufacturée
########################################
r, t = sp.symbols('r t')

# Définition de la solution manufacturée
C_chapeau = C_E * (r**2 / R**2) + sp.exp(t / Tf) * (1 - r**2 / R**2)

# Transformation de la solution manufacturée sympy en fonction traçable par matplotlib
C_chapeau = sp.lambdify([r,t], C_chapeau, modules=['numpy'])

# Tracé de la solution manufacturée (question b)
x_values = np.linspace(0, R, 500)
plt.figure()

for t in ([0,Tf/2,Tf]):
    y_values = C_chapeau(x_values,t)
    plt.plot(x_values, y_values, label=f't={t} s')

plt.title("Solution manufacturée fonction de r à différents instants")
plt.xlabel('r (m)')
plt.ylabel('Solution manufacturée')
plt.legend()
plt.grid()
plt.savefig('Devoir 2/results/solution_manufacturée.png')
plt.show()

##############
# Terme source
##############
# Tracé du terme source (question c)
r_values = np.linspace(0, R, 100)
plt.figure()

for t in ([0,Tf/2,Tf]):
    S_values = np.zeros(len(r_values))
    for i in range(len(r_values)):
        S_values[i] = functions.terme_source(r_values[i],t)
    plt.plot(r_values, S_values, label=f't={t} s')

plt.title("Terme source en fonction de r à différents instants")
plt.xlabel('r (m)')
plt.ylabel('Terme source')
plt.legend()
plt.grid()
plt.savefig('Devoir 2/results/terme_source.png')
plt.show()


##############################################################################
# Analyse de convergence spatiale avec delta temporel fixé à 3,1536e7 s (1 an)
##############################################################################
list_N_spatiaux=[3, 4, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]

delta_t = 365*24*3600      # 1 an en s
N_temporel = (Tf/delta_t)+1
N_temporel = int(N_temporel)

vectDelta_r = np.zeros(len(list_N_spatiaux))
vectL1 = np.zeros(len(list_N_spatiaux))
vectL2 = np.zeros(len(list_N_spatiaux))
vectLinf = np.zeros(len(list_N_spatiaux))

for i in (range(len(list_N_spatiaux))):
    delta_r = R/(list_N_spatiaux[i]-1)
    N_spatial = list_N_spatiaux[i]
    
    C_i_n, r_i, t_i=functions.Concentrations(delta_r, delta_t)
    C_spatiaux = C_i_n[-1,:]
    C_temporels = C_i_n[:,0]
    C_exact_spatial = C_chapeau(r_i, Tf)
    C_exacte_temporel = C_chapeau(0, t_i)
    vectDelta_r[i]=delta_r

    L1=errors.ErreurL1(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel,N_spatial,N_temporel)
    vectL1[i]=L1
    L2=errors.ErreurL2(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel,N_spatial,N_temporel)
    vectL2[i]=L2
    Linf=errors.ErreurLinf(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel)
    vectLinf[i]=Linf

plt.figure()

# L1
plt.plot(vectDelta_r, vectL1, 'o', color='red', label='L1')
slope_L1, intercept_L1, r_value_L1, p_value_L1, std_err_L1 = linregress(np.log(vectDelta_r[6:]), np.log(vectL1[6:]))
y_pred_L1 =  np.exp(intercept_L1) * vectDelta_r[6:]**slope_L1
plt.plot(vectDelta_r[6:], y_pred_L1, '--', color='red', label=f'Ordre de convergence L1: {slope_L1}')

# L2
plt.plot(vectDelta_r, vectL2, 'o', color='green', label='L2')
slope_L2, intercept_L2, r_value_L2, p_value_L2, std_err_L2 = linregress(np.log(vectDelta_r[6:]), np.log(vectL2[6:]))
y_pred_L2 = np.exp(intercept_L2) * vectDelta_r[6:]**slope_L2
plt.plot(vectDelta_r[6:], y_pred_L2, '--', color='green', label=f'Ordre de convergence L2: {slope_L2}')

# Linf
plt.plot(vectDelta_r, vectLinf, 'o', color='blue', label='Linf')
slope_L3, intercept_L3, r_value_L3, p_value_L3, std_err_L3 = linregress(np.log(vectDelta_r[6:]), np.log(vectLinf[6:]))
y_pred_Linf = np.exp(intercept_L3) * vectDelta_r[6:]**slope_L3
plt.plot(vectDelta_r[6:], y_pred_Linf, '--', color='blue', label=f'Ordre de convergence Linf: {slope_L3}')

plt.title('Erreurs L1, L2 et Linf en fonction du nombre de points N spatiaux')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta r (m)')
plt.ylabel('Erreur (mol/m^3)')
plt.grid(True, which="both", ls="--")
plt.savefig('Devoir 2/results/erreurs_spatiales.png')
plt.show()


##############################################################################
# Analyse de convergence temporelle avec delta spatial fixé à 0,125 m (5 points)
##############################################################################
list_N_temporel=[100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000]

delta_r = R/(20-1)      # en m
N_spatial = 20

vectDelta_t = np.zeros(len(list_N_temporel))
vectL1 = np.zeros(len(list_N_temporel))
vectL2 = np.zeros(len(list_N_temporel))
vectLinf = np.zeros(len(list_N_temporel))

for i in (range(len(list_N_temporel))):
    delta_t = Tf/(list_N_temporel[i]-1)
    N_temporel = list_N_temporel[i]
    
    C_i_n, r_i, t_i=functions.Concentrations(delta_r, delta_t)
    C_spatiaux = C_i_n[-1,:]
    C_temporels = C_i_n[:,0]
    C_exact_spatial = C_chapeau(r_i, Tf)
    C_exacte_temporel = C_chapeau(0, t_i)
    vectDelta_t[i]=delta_t

    L1=errors.ErreurL1(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel,N_spatial,N_temporel)
    vectL1[i]=L1
    L2=errors.ErreurL2(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel,N_spatial,N_temporel)
    vectL2[i]=L2
    Linf=errors.ErreurLinf(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel)
    vectLinf[i]=Linf

plt.figure()

# L1
plt.plot(vectDelta_t, vectL1, 'o', color='red', label='L1')
slope_L1, intercept_L1, r_value_L1, p_value_L1, std_err_L1 = linregress(np.log(vectDelta_t[:3]), np.log(vectL1[:3]))
y_pred_L1 =  np.exp(intercept_L1) * vectDelta_t[:3]**slope_L1
plt.plot(vectDelta_t[:3], y_pred_L1, '--', color='red', label=f'Ordre de convergence L1: {slope_L1}')

# L2
plt.plot(vectDelta_t, vectL2, 'o', color='green', label='L2')
slope_L2, intercept_L2, r_value_L2, p_value_L2, std_err_L2 = linregress(np.log(vectDelta_t[:3]), np.log(vectL2[:3]))
y_pred_L2 = np.exp(intercept_L2) * vectDelta_t[:3]**slope_L2
plt.plot(vectDelta_t[:3], y_pred_L2, '--', color='green', label=f'Ordre de convergence L2: {slope_L2}')

# Linf
plt.plot(vectDelta_t, vectLinf, 'o', color='blue', label='Linf')
slope_L3, intercept_L3, r_value_L3, p_value_L3, std_err_L3 = linregress(np.log(vectDelta_t[:3]), np.log(vectLinf[:3]))
y_pred_Linf = np.exp(intercept_L3) * vectDelta_t[:3]**slope_L3
plt.plot(vectDelta_t[:3], y_pred_Linf, '--', color='blue', label=f'Ordre de convergence Linf: {slope_L3}')

plt.title('Erreurs L1, L2 et Linf en fonction du nombre de points N temporels')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta t (s)')
plt.ylabel('Erreur (mol/m^3)')
plt.grid(True, which="both", ls="--")
plt.savefig('Devoir 2/results/erreurs_temporels.png')
plt.show()


###################################################
# Tracer de la concentration pour N=11 (question f)
###################################################
N_spatial = 11
N_temporel = 100
delta_r = R/(N_spatial-1)
delta_t = Tf/(N_temporel-1)

C_i_n, r_i, t_i=functions.Concentrations(delta_r, delta_t)
plt.figure()
plt.plot(r_i, C_i_n[-1,:], label='Concentration à t=4e9 s')
plt.title("Concentration en fonction de r pour N=11")
plt.xlabel('r (m)')
plt.ylabel('Concentration (mol/m^3)')
plt.legend()
plt.grid()
plt.savefig('Devoir 2/results/trace_question_f.png')
plt.show()
