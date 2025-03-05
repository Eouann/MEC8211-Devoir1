"""
Fichier de calcul des fonction analytique et des concentrations avec la methode des differences finies
"""
"""
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import config2

###########################
# Définition des constantes
###########################
N_spatial = config2.N_Spatial    # nombre de points pour la méthode des différences finies
N_Temporel = config2.N_Temporel 
D_EFF = config2.D_EFF        # en m^2/s
C_E = config2.C_E            # en mol/m^3
D = config2.D             # en m
R = config2.R            # en m
S = config2.S             # en mol/m^3/s
Tf = config2.Tf        # Temps caractéristique
k = config2.k            # Coefficient de réaction
"""
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import config2
###########################
# Définition des constantes
###########################
N_spatial = config2.N_Spatial    # nombre de points pour la méthode des différences finies
N_Temporel = config2.N_Temporel 
D_EFF = config2.D_EFF        # en m^2/s
C_E = config2.C_E            # en mol/m^3
D = config2.D             # en m
R = config2.R            # en m
S = config2.S             # en mol/m^3/s
Tf = config2.Tf        # Temps caractéristique
k = config2.k            # Coefficient de réaction

# Définition des symboles
r, t = sp.symbols('r t')

# Définition de la solution manufacturée correcte
C_r = C_E * (r**2 / R**2) + sp.exp(t / Tf) * (1 - r**2 / R**2)

# Fonction pour calculer le terme source
def Terme_Source(r,t):
   S = np.exp(t/Tf)*(1/Tf)*(1-r**2/R**2) - D_EFF* (2*C_E/R**2-2*np.exp(t/Tf)/R**2+ 1/r*(2*C_E*r/R**2-np.exp(t/Tf)*2*r/R**2))+k*(C_E*r**2/R**2+np.exp(t/Tf)*(1-r**2/R**2))
   return S

# Discrétisation
delta_r = R / (N_spatial - 1)
delta_t = Tf / (N_Temporel - 1)
r_i = np.linspace(0, R, N_spatial)
t_i = np.linspace(0, Tf, N_Temporel)

# Initialisation des matrices et vecteurs
C_i = np.zeros(N_spatial)
matA = np.zeros((N_spatial, N_spatial))
vectB = np.zeros(N_spatial)

# Conditions aux limites
# Condition de Neumann à r = 0 avec Gear
matA[0, 0] = -3
matA[0, 1] = 4
matA[0, 2] = -1

# Condition de Dirichlet à r = R
matA[-1, -1] = 1
vectB[-1] = C_E

# Construction de la matrice A (Coefficients de l'équation aux différences finies)
for i in range(1, N_spatial - 1):
    matA[i, i - 1] = D_EFF * (1 / delta_r**2 - 1 / (2 * r_i[i] * delta_r))
    matA[i, i] = -(k + (1 / delta_t) + 2 * (D_EFF / delta_r**2))
    matA[i, i + 1] = D_EFF * (1 / delta_r**2 + 1 / (2 * r_i[i] * delta_r))

# Boucle temporelle pour résoudre C_i à chaque pas de temps
for t in t_i:
    for j in range(1, N_spatial - 1):
        vectB[j] = -C_i[j] / delta_t + Terme_Source(r_i[j], t)

    # Résolution du système matriciel A * C_new = B
    C_new = np.linalg.solve(matA, vectB)
    C_i = C_new  # Mise à jour des concentrations pour l'itération suivante

# Tracer C_i en fonction de r pour le dernier pas de temps, en affichant les points
plt.figure(figsize=(10, 6))
plt.scatter(r_i, C_i, label=f'Temps = {t_i[-1]:.1e} s', color='b', marker='o')

plt.title("Concentration en fonction de r au dernier temps")
plt.xlabel('r (m)')
plt.ylabel('Concentration (mol/m^3)')
plt.legend()
plt.grid(True)
plt.show()
