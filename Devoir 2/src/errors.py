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
def ErreurL1(Ci,Ci_exact,N_spatial,N_temporel):
    """Fonction de calcul de l'erreur L1 discrète"""
    normeL1 = 1/N_spatial*np.sum(np.abs(Ci-Ci_exact))+1/N_temporel*np.sum(np.abs(Ci-Ci_exact))
    return normeL1


###################
#Erreur L2 discrète
###################
def ErreurL2(Ci,Ci_exact,N_spatial,N_temporel):
    """Fonction de calcul de l'erreur L2 discrète"""
    normeL2 = np.sqrt(1/N_spatial*np.sum(np.abs(Ci-Ci_exact)**2))+np.sqrt(1/N_temporel*np.sum(np.abs(Ci-Ci_exact)**2))
    return normeL2


#####################
#Erreur Linf discrète
#####################
def ErreurLinf(Ci,Ci_exact):
    """Fonction de calcul de l'erreur Linf discrète"""
    normeLinf = np.max(np.abs(Ci-Ci_exact))+np.max(np.abs(Ci-Ci_exact))
    return normeLinf
