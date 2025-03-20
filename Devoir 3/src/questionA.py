"""
Question A : Determination de u_num grace au GCI (Grid Convergence Index)
"""


# Importation des bibliothèques
import numpy as np


# Résultats de simulations effectuées avec les fichiers matlab (seed=4852)
r = 1.6 # Choisi arbitrairement
list_Nx = [200,320,512]
liste_dx = [2.00e-6,1.25e-6,7.81e-7]
liste_k = [27.9845,27.7591,27.1026] # f_r^2.h, f_r.h, f_h
p_f = 2 # D'après l'énoncé


# Définition de la fonction u_num
def u_num(liste_f,liste_dx,r,p_f):
    """Fonction de calcul de l'incertitude numérique u_num grace au GCI"""
    p_chapeau = np.log((liste_f[2]-liste_f[1])/(liste_f[1]-liste_f[0]))/np.log(r)
    if np.abs((p_chapeau-p_f)/p_f) <= 0.1:
        GCI = 1.25/(r**p_f-1)*np.abs(liste_f[1]-liste_f[0])
    else:
        p=min(max(0.5,p_chapeau),p_f)
        GCI = 3/(r**p-1)*np.abs(liste_f[1]-liste_f[0])
    u_num = GCI/2
    return u_num


# Affichage de la valeur de u_num
print("La valeur de u_num est : ",u_num(liste_k,liste_dx,r,p_f))