"""
Question A : Determination de u_num grace au GCI (Grid Convergence Index)
"""


# Importation des bibliothèques
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


# Résultats de simulations effectuées avec les fichiers matlab (seed=4852)
r = 2 # Choisi arbitrairement
list_Nx = np.array([100,200,400])
liste_dx = np.array([8e-6,4e-6,2e-6])
liste_k = np.array([49.1695,36.7123,32.2762]) # f_r^2.h, f_r.h, f_h
p_f = 2 # D'après l'énoncé


# Définition de la fonction u_num
def u_num(liste_f,liste_dx,r,p_f):
    """Fonction de calcl de u_num, l'incertitude numérique, grace au GCI (Grid Convergence Index)"""
    p_chapeau = np.log((liste_f[0]-liste_f[1])/(liste_f[1]-liste_f[2]))/np.log(r)
    if np.abs((p_chapeau-p_f)/p_f) <= 0.1:
        GCI = 1.25/(r**p_f-1)*np.abs(liste_f[1]-liste_f[2])
    else:
        p=min(max(0.5,p_chapeau),p_f)
        GCI = 3/(r**p-1)*np.abs(liste_f[1]-liste_f[2])
    u_num = GCI/2
    return u_num


# Affichage de la valeur de u_num
print("La valeur de u_num est : ",u_num(liste_k,liste_dx,r,p_f))