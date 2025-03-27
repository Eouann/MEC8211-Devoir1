"""
Question D et E : Calcul de E l'erreur de simulation et delta_modele l'erreur du modèle d'aprés la norme V&V20 
"""


# Importation des bibliothèques
import numpy as np


# Importation des données
S = XXX
D = 80.6 # En µm^2, d'après l'énoncé
u_num = 3.6801030800014876 # En µm^2, d'après la question A
u_input = XXX
u_D = XXX
k = 2 # Car u_num a ete déterminé grace au GCI


# Définition de la fonction E
def Calcul_E(S,D):
    """Fonction de calcul de E l'erreur de simulation"""
    E = S-D
    return E


# Execution du calcul
E=Calcul_E(S,D)
print("La valeur de E est : ",E)


# Calcul de u_val, l'incertitude de validation
def Calcul_u_val(u_num,u_input,u_D):
    """Fonction de calcul de u_val l'incertitude de validation"""
    u_val = np.sqrt(u_num**2+u_input**2+u_D**2)
    return u_val


# Intervalle de delta_model
intervalle_sup=E+k*Calcul_u_val(u_num,u_input,u_D)
intervalle_inf=E-k*Calcul_u_val(u_num,u_input,u_D)

print("Delta_model, l'erreur du modèle se trouve dans l'intervalle [",intervalle_inf,";",intervalle_sup,"] avec une confiance de 95,4 %")