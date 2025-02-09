#############################
#Importation de bibliothèques
#############################
import numpy as np
import config
import functions


############################
# Importation des constantes
############################
Ntot = config.Ntot      # nombre de point pour la méthodes des différences finies


###################
#Erreur L1 discrète
###################
def ErreurL1(Ci,Ci_exact):
    # Fonction de calcul de l'erreur L1 discrète
    normeL1 = 1/Ntot*np.sum(np.abs(Ci-Ci_exact))
    return normeL1


###################
#Erreur L2 discrète
###################
def ErreurL2(Ci,Ci_exact):
    # Fonction de calcul de l'erreur L2 discrète
    normeL2 = (1/Ntot*np.sum(np.abs(Ci-Ci_exact)**2))**0.5
    return normeL2


#####################
#Erreur Linf discrète
#####################
def ErreurLinf(Ci,Ci_exact):
    # Fonction de calcul de l'erreur Linf discrète
    normeLinf = np.max(np.abs(Ci-Ci_exact))
    return normeLinf