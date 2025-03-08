"""
Fichier de calcul du terme source et des concentrations avec la methode des differences finies
"""


#############################
#Importation de bibliothèques
#############################
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import config


###########################
# Définition des constantes
###########################
D_EFF = config.D_EFF            # en m^2/s
C_E = config.C_E                # en mol/m^3
D = config.D                    # en m
R = config.R                    # en m
Tf = config.Tf                  # Temps caractéristique
k = config.k                    # Coefficient de réaction


########################
# Calcul du terme source
########################
def terme_source(r,t):
    """Définition de la fon,ctopn de terme source en implicite"""
    S = np.exp(t/Tf) * (1/Tf) * (1 - r**2/R**2) \
       - D_EFF * (2 * C_E / R**2 - 2 * np.exp(t/Tf) / R**2 + 1 / r \
           * ( 2 *C_E * r / R**2 - np.exp(t/Tf) * 2 * r/R**2))+k\
           * (C_E * r**2 / R**2 + np.exp(t/Tf) * (1 - r**2/R**2))
    return S


###########################
# Calcul des concentrations
###########################
def Concentrations(delta_r, delta_t):
    # Discrétisation
    N_spatial = (R / delta_r) + 1
    N_spatial = int(N_spatial)
    N_temporel = (Tf / delta_t) + 1
    N_temporel = int(N_temporel)
    r_i = np.linspace(0, R, N_spatial)
    t_i = np.linspace(0, Tf, N_temporel)

    # Initialisation des matrices et vecteurs
    C_i = np.zeros(N_spatial)
    matA = np.zeros((N_spatial, N_spatial))
    vectB = np.zeros(N_spatial)

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
        matA[i, i] = - (k + (1 / delta_t) + 2 * (D_EFF / delta_r**2))
        matA[i, i + 1] = D_EFF * (1 / delta_r**2 + 1 / (2 * r_i[i] * delta_r))

    # Boucle temporelle pour résoudre C_i à chaque pas de temps
    for t in t_i:
        for j in range(1, N_spatial - 1):
            vectB[j] = - C_i[j] / delta_t + terme_source(r_i[j], t)

        # Résolution du système matriciel A * C_new = B
        C_new = np.linalg.solve(matA, vectB)
        C_i = C_new     # Mise à jour des concentrations pour l'itération suivante
    return C_i, r_i, t_i
