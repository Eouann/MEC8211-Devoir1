"""
Fichier principal d'execution du devoir 1
"""


#############################
#Importation de bibliothèques
#############################
import numpy as np
import matplotlib.pyplot as plt
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


# Importation des coeficients pour Ntot points CAS 1 et 2
a1,b1,c1,a2,b2,c2,r_i,delta_r = functions.Coefficients(N_TOT)


# Calcul des concentrations pour les deux cas et pour le cas analytique
Concentrations_CAS1=functions.Concentrations(a1,b1,c1,N_TOT,1)
Concentrations_CAS2=functions.Concentrations(a2,b2,c2,N_TOT,2)
x_values = functions.x_values
y_values = functions.y_values


# Tracé des fonctions analytiques et différences finies

plt.figure()
plt.plot(x_values, y_values, label='Analytique')
plt.plot(r_i, Concentrations_CAS1, 'o', color='red', label='Cas 1 - numérique')
plt.plot(r_i, Concentrations_CAS2, 'o', color='green',label='Cas 2 - numérique')
plt.xlabel('r (m)')
plt.ylabel('Concentration (mol/m^3)')
plt.grid()
plt.legend()
plt.savefig("Devoir 1/results/fonctions.png")
plt.show()

# Tracé des erreurs L1, L2 et Linf CAS1
vectDelta_r = []
vectL1 = []
vectL2 = []
vectLinf = []

for N in range(1000, 10000, 800):
    a1,b1,c1,a2,b2,c2,r_i,delta_r = functions.Coefficients(N)
    vectDelta_r.append(delta_r)
    Concentrations_CAS1=functions.Concentrations(a1,b1,c1,N,1)
    x_Nvalues,y_Nvalues=functions.C_analytique_N(N)
    L1=errors.ErreurL1(Concentrations_CAS1,x_Nvalues,N)
    vectL1.append(L1)
    L2=errors.ErreurL2(Concentrations_CAS1,x_Nvalues,N)
    vectL2.append(L2)
    Linf=errors.ErreurLinf(Concentrations_CAS1,x_Nvalues)
    vectLinf.append(Linf)

plt.figure()
#plt.plot(vectDelta_r, vectL1, 'o', color='red', label='L1')
plt.plot(vectDelta_r, vectL2, 'o', color='green', label='L2')
#plt.plot(vectDelta_r, vectLinf, 'o', color='blue', label='Linf')
# Regression lineaire de L1
#a,b = np.polyfit(np.log(vectDelta_r[8:]), np.log(vectL1[8:]), 1)
#plt.plot(vectDelta_r[8:], np.exp(a*np.log(vectDelta_r[8:])+b), color='red', label='Regression L1')
# Regression lineaire de L2
#c,d = np.polyfit(np.log(vectDelta_r[8:]), np.log(vectL2[8:]), 1)
#plt.plot(vectDelta_r[8:], np.exp(c*np.log(vectDelta_r[8:])+d), color='green', label='Regression L2')
# Regression lineaire de Linf
#e,f = np.polyfit(np.log(vectDelta_r), np.log(vectLinf), 1)
#plt.plot(vectDelta_r, np.exp(e*np.log(vectDelta_r)+f), label='Regression Linf')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta r (m)')
plt.ylabel('Erreur (mol/m^3)')
plt.grid()
plt.savefig("Devoir 1/results/erreurs.png")
plt.show()
