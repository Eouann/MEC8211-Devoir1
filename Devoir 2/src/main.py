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
N_spatial = config.N_Spatial    # Nombre de points spatiaux pour la méthode des différences finies
N_Temporel = config.N_Temporel  # Nombre de points temporels pour la méthode des différences finies
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
list_delta_r=[0.25, 0.166666667, 0.125, 0.055555556, 
              0.026315789, 0.017241379, 0.010204082, 
              0.005050505, 0.002512563, 0.001672241, 
              0.001002004, 0.000500501]
delta_t = (3.1536e7)      # 1 an en s
vectDelta_r = np.zeros(len(list_delta_r))
vectL1 = np.zeros(len(list_delta_r))
vectL2 = np.zeros(len(list_delta_r))
vectLinf = np.zeros(len(list_delta_r))

for i in (range(len(list_delta_r))):
    delta_r = list_delta_r[i]
    N = R / delta_r + 1
    N = int(N)
    
    C_i, r_i, t_i=functions.Concentrations(delta_r, delta_t)
    y_values = C_chapeau(r_i, Tf)
    vectDelta_r[i]=delta_r

    L1=errors.ErreurL1(C_i,y_values,N)
    vectL1[i]=L1
    L2=errors.ErreurL2(C_i,y_values,N)
    vectL2[i]=L2
    Linf=errors.ErreurLinf(C_i,y_values)
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
list_delta_t=[935000000, 415555555.55555600, 128965517.24137900, 76326530.61224490,
37777777.77777780, 12508361.20401340, 7494989.97995992, 3743743.74374374]
delta_r = 0.125      # en m
vectDelta_t = np.zeros(len(list_delta_t))
vectL1 = np.zeros(len(list_delta_t))
vectL2 = np.zeros(len(list_delta_t))
vectLinf = np.zeros(len(list_delta_t))

for i in (range(len(list_delta_t))):
    delta_t = list_delta_t[i]
    N = Tf / delta_t + 1
    N = int(N)
    
    C_i, r_i, t_i=functions.Concentrations(delta_r, delta_t)
    y_values = C_chapeau(0.25, t_i)
    vectDelta_t[i]=delta_t

    L1=errors.ErreurL1(C_i[2],y_values,N)
    vectL1[i]=L1
    L2=errors.ErreurL2(C_i[2],y_values,N)
    vectL2[i]=L2
    Linf=errors.ErreurLinf(C_i[2],y_values)
    vectLinf[i]=Linf

plt.figure()

# L1
plt.plot(vectDelta_t, vectL1, 'o', color='red', label='L1')
slope_L1, intercept_L1, r_value_L1, p_value_L1, std_err_L1 = linregress(np.log(vectDelta_t[4:]), np.log(vectL1[4:]))
y_pred_L1 =  np.exp(intercept_L1) * vectDelta_t[4:]**slope_L1
plt.plot(vectDelta_t[4:], y_pred_L1, '--', color='red', label=f'Ordre de convergence L1: {slope_L1}')

# L2
plt.plot(vectDelta_t, vectL2, 'o', color='green', label='L2')
slope_L2, intercept_L2, r_value_L2, p_value_L2, std_err_L2 = linregress(np.log(vectDelta_t[4:]), np.log(vectL2[4:]))
y_pred_L2 = np.exp(intercept_L2) * vectDelta_t[4:]**slope_L2
plt.plot(vectDelta_t[4:], y_pred_L2, '--', color='green', label=f'Ordre de convergence L2: {slope_L2}')

# Linf
plt.plot(vectDelta_t, vectLinf, 'o', color='blue', label='Linf')
slope_L3, intercept_L3, r_value_L3, p_value_L3, std_err_L3 = linregress(np.log(vectDelta_t[4:]), np.log(vectLinf[4:]))
y_pred_Linf = np.exp(intercept_L3) * vectDelta_t[4:]**slope_L3
plt.plot(vectDelta_t[4:], y_pred_Linf, '--', color='blue', label=f'Ordre de convergence Linf: {slope_L3}')

plt.title('Erreurs L1, L2 et Linf en fonction du nombre de points N temporels')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta t (s)')
plt.ylabel('Erreur (mol/m^3)')
plt.grid(True, which="both", ls="--")
plt.savefig('Devoir 2/results/erreurs_temporels.png')
plt.show()
