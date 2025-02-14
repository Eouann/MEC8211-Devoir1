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
plt.plot(x_values, y_values, label='Analytique')
plt.plot(r_i, Concentrations_CAS1, 'o', color='red', label='Cas 1 - numérique')
plt.plot(r_i, Concentrations_CAS2, 'o', color='green',label='Cas 2 - numérique')
plt.title('Concentration en fonction de r')
plt.xlabel('r (m)')
plt.ylabel('Concentration (mol/m^3)')
plt.grid()
plt.legend()
plt.savefig('Devoir 1/results/fonctions.png')
plt.show()

# Tracé des erreurs L1, L2 et Linf CAS1
vectDelta_r = []
vectL1 = []
vectL2 = []
vectLinf = []

for N in ([3, 4, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]):
    Concentrations_CAS1,delta_r,r_i=functions.Concentrations(N,1)
    vectDelta_r.append(delta_r)
    x_Nvalues,y_Nvalues=functions.C_analytique_N(N)
    L1=errors.ErreurL1(Concentrations_CAS1,x_Nvalues,N)
    vectL1.append(L1)
    L2=errors.ErreurL2(Concentrations_CAS1,x_Nvalues,N)
    vectL2.append(L2)
    Linf=errors.ErreurLinf(Concentrations_CAS1,x_Nvalues)
    vectLinf.append(Linf)

plt.figure()

# L1
plt.plot(vectDelta_r, vectL1, 'o', color='red', label='L1')
plt.title('Erreur L1 en fonction du pas')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta r (m)')
plt.ylabel('Erreur (mol/m^3)')
plt.grid(True, which="both", ls="--")
'''
# L2
plt.plot(vectDelta_r, vectL2, 'o', color='green', label='L2')
plt.title('Erreur L2 en fonction du pas')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta r (m)')
plt.grid(True, which="both", ls="--")

# Linf
plt.plot(vectDelta_r, vectLinf, 'o', color='blue', label='Linf')
plt.title('Erreur Linf en fonction du pas')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta r (m)')
plt.grid(True, which="both", ls="--")
'''
plt.savefig('Devoir 1/results/erreurs.png')
plt.show()
