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
Ntot = config.Ntot      # nombre de point pour la méthodes des différences finies
Deff = config.Deff      # en m^2/s
Ce = config.Ce          # en mol/m^3
D = config.D            # en m
R = config.R            # en m
S = config.S            # en mol/m^3/s


# Importation des coeficients pour Ntot points CAS 1 et 2
a1,b1,c1,a2,b2,c2,r_i,delta_r = functions.Coefficients(Ntot)


# Calcul des concentrations pour les deux cas et pour le cas analytique
Concentrations_CAS1=functions.Concentrations(a1,b1,c1)
Concentrations_CAS2=functions.Concentrations(a2,b2,c2)
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
plt.show()

# Tracé des erreurs L1, L2 et Linf CAS1
vectDelta_r = []
vectL1 = []
vectL2 = []
vectLinf = []

for N in range(10, 1010, 100):
    a1,b1,c1,a2,b2,c2,r_i,delta_r = functions.Coefficients(N)
    vectDelta_r.append(delta_r)
    Concentrations_CAS1=functions.Concentrations(a1,b1,c1)
    L1=errors.ErreurL1(Concentrations_CAS1,functions.x_Ntotvalues)
    vectL1.append(L1)
    L2=errors.ErreurL2(Concentrations_CAS1,functions.x_Ntotvalues)
    vectL2.append(L2)
    Linf=errors.ErreurLinf(Concentrations_CAS1,functions.x_Ntotvalues)
    vectLinf.append(Linf)



plt.figure()
plt.plot(vectDelta_r, vectL1, 'o', color='red', label='L1')
plt.plot(vectDelta_r, vectL2, 'o', color='green', label='L2')
plt.plot(vectDelta_r, vectLinf, 'o', color='blue', label='Linf')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta r (m)')
plt.ylabel('Erreur (mol/m^3)')
plt.grid()
plt.show()