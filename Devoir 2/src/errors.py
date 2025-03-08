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
def ErreurL1(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel,N_spatial,N_temporel):
    """Fonction de calcul de l'erreur L1 discrète"""
    normeL1 = 1/N_spatial*np.sum(np.abs(C_spatiaux-C_exact_spatial))+1/N_temporel*np.sum(np.abs(C_temporels-C_exacte_temporel))
    return normeL1


###################
#Erreur L2 discrète
###################
def ErreurL2(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel,N_spatial,N_temporel):
    """Fonction de calcul de l'erreur L2 discrète"""
    normeL2 = np.sqrt(1/N_spatial*np.sum(np.abs(C_spatiaux-C_exact_spatial)**2))+np.sqrt(1/N_temporel*np.sum(np.abs(C_temporels-C_exacte_temporel)**2))
    return normeL2


#####################
#Erreur Linf discrète
#####################
def ErreurLinf(C_spatiaux,C_temporels,C_exact_spatial,C_exacte_temporel):
    """Fonction de calcul de l'erreur Linf discrète"""
    normeLinf = np.max(np.abs(C_spatiaux-C_exact_spatial))+np.max(np.abs(C_temporels-C_exacte_temporel))
    return normeLinf
