"""
Fichier principal d'execution du devoir 1
"""


#############################
#Importation de bibliothèques
#############################
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import config
import functions
import errors


############################
# Importation des constantes
############################
N_TOT = config.N_TOT    # nombre de point pour la méthodes des différences finies
D_EFF = config.D_EFF    # en m^2/s
C_E = config.C_E        # en mol/m^3
D = config.D            # en m
R = config.R            # en m
S = config.S            # en mol/m^3/s


# Calcul des concentrations pour les deux cas et pour le cas analytique
Concentrations_CAS1,delta_r,r_i=functions.Concentrations(N_TOT,1)
Concentrations_CAS2,delta_r,r_i=functions.Concentrations(N_TOT,2)
x_values = functions.x_values
y_values = functions.y_values


# Tracé des fonctions analytiques et différences finies

plt.figure()
plt.plot(x_values, y_values, color='blue', label='Analytique')
plt.plot(r_i, Concentrations_CAS1, 'o', color='red', label='Cas 1 - numérique')
plt.plot(r_i, Concentrations_CAS2, 'o', color='green',label='Cas 2 - numérique')
plt.title('Concentration en fonction de r')
plt.xlabel('r (m)')
plt.ylabel('Concentration (mol/m^3)')
plt.grid()
plt.legend()
plt.savefig('Devoir 1/results/fonctions.png')
plt.show()

# Analyse de convergence cas 1
listN=[3, 4, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]
vectDelta_r = np.zeros(len(listN))
vectL1 = np.zeros(len(listN))
vectL2 = np.zeros(len(listN))
vectLinf = np.zeros(len(listN))

for i in (range(len(listN))):
    N = listN[i]
    Concentrations_CAS1,delta_r,r_i=functions.Concentrations(N,1)
    vectDelta_r[i]=delta_r
    x_Nvalues,y_Nvalues=functions.C_analytique_N(N)
    L1=errors.ErreurL1(Concentrations_CAS1,y_Nvalues,N)
    vectL1[i]=L1
    L2=errors.ErreurL2(Concentrations_CAS1,y_Nvalues,N)
    vectL2[i]=L2
    Linf=errors.ErreurLinf(Concentrations_CAS1,y_Nvalues)
    vectLinf[i]=Linf

plt.figure()

# L1
plt.plot(vectDelta_r, vectL1, 'o', color='red', label='L1')
slope_L1, intercept_L1, r_value_L1, p_value_L1, std_err_L1 = linregress(vectDelta_r[6:], vectL1[6:])
y_pred_L1 = slope_L1 * vectDelta_r[6:] + intercept_L1
plt.plot(vectDelta_r[6:], y_pred_L1, '--', color='red', label=f'Ordre de convergence L1: {np.log(slope_L1)}')

# L2
plt.plot(vectDelta_r, vectL2, 'o', color='green', label='L2')
slope_L2, intercept_L2, r_value_L2, p_value_L2, std_err_L2 = linregress(vectDelta_r[6:], vectL2[6:])
y_pred_L2 = slope_L2 * vectDelta_r[6:] + intercept_L2
plt.plot(vectDelta_r[6:], y_pred_L2, '--', color='green', label=f'Ordre de convergence L2: {np.log(slope_L2)}')

# Linf
plt.plot(vectDelta_r, vectLinf, 'o', color='blue', label='Linf')
slope_L3, intercept_L3, r_value_L3, p_value_L3, std_err_L3 = linregress(vectDelta_r[6:], vectLinf[6:])
y_pred_Linf = slope_L3 * vectDelta_r[6:] + intercept_L3
plt.plot(vectDelta_r[6:], y_pred_Linf, '--', color='blue', label=f'Ordre de convergence Linf: {np.log(slope_L3)}')

plt.title('Erreurs L1, L2 et Linf en fonction du nombre de points N')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta r (m)')
plt.ylabel('Erreur (mol/m^3)')
plt.grid(True, which="both", ls="--")
plt.savefig('Devoir 1/results/erreurs.png')
plt.show()
