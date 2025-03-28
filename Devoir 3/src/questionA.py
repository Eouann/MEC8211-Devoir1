"""
Question A : Calcul de u_num grace au GCI (Grid Convergence Index)
"""


# Importation des bibliothèques
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


# Vérification de la convergence asymptotique
k_1 = np.array([22.6343,30.5561,17.3802,30.5561,23.0971,17.8138,53.9237,38.582,30.4274,40.7066])    # Liste de k simulés pour Nx=133 et dx=1.5e-6
k_2 = np.array([23.7583,39.275,40.2556,23.1456,27.7324,41.4252,30.1901,44.6739,14.4448,13.4319])    # Liste de k simulés pour Nx=167 et dx=1.2e-6
k_3 = np.array([12.8425,34.2975,42.7621,23.937,26.6145,41.8116,21.6885,50.4074,22.0472,13.1])       # Liste de k simulés pour Nx=222 et dx=0.9e-6

k_1_moyenne = np.mean(k_1)
k_2_moyenne = np.mean(k_2)
k_3_moyenne = np.mean(k_3)

liste_k_moyen = np.array([k_1_moyenne,k_2_moyenne,k_3_moyenne])
liste_dx = np.array([1.5e-6,1.2e-6,0.9e-6])

plt.plot(liste_dx, liste_k_moyen, 'o', color='red', label='k moyen')
slope, intercept, r_value, p_value, std_err = linregress(np.log(liste_dx), np.log(liste_k_moyen))
y_pred =  np.exp(intercept) * liste_dx**slope
plt.plot(liste_dx, y_pred, '--', color='red', label=f'Ordre de convergence : {slope}')
plt.title('Vérification de la convergence asymptotique')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delta x (m)')
plt.ylabel('k (μm^3)')
plt.grid(True, which="both", ls="--")
plt.savefig('Devoir 3/results/analyse-de-convergence.png')
plt.show()

print(slope,intercept)

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