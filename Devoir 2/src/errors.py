"""
Fichiers de calcul des erreurs L1, L2 et Linf discrètes
"""


#############################
#Importation de bibliothèques
#############################
import numpy as np


###################
#Erreur L1 discrète
###################
def ErreurL1(C_numerique,C_exacte,N_spatial,N_temporel):
    """Fonction de calcul de l'erreur L1 discrète"""
    normeL1 = 1/N_spatial*1/N_temporel*np.sum(np.abs(C_numerique-C_exacte))
    return normeL1


###################
#Erreur L2 discrète
###################
def ErreurL2(C_numerique,C_exacte,N_spatial,N_temporel):
    """Fonction de calcul de l'erreur L2 discrète"""
    normeL2 = np.sqrt(1/N_spatial*1/N_temporel*np.sum(np.abs(C_numerique-C_exacte)**2))
    return normeL2


#####################
#Erreur Linf discrète
#####################
def ErreurLinf(C_numerique,C_exacte,):
    """Fonction de calcul de l'erreur Linf discrète"""
    normeLinf = np.max(np.abs(C_numerique-C_exacte))
    return normeLinf
