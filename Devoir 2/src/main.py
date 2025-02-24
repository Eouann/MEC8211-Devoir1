"""
Fichier principal d'execution du devoir 1
"""


#############################
#Importation de bibliothèques
#############################
import numpy as np
import matplotlib.pyplot as plt
import config
import functions
import errors


############################
# Importation des constantes
############################
N_TOT = config.N_TOT    # nombre de point pour la méthodes des différences finies
dt = config.dt          # pas de temps en s
t_fin = config.t_fin    # temps final en s
D_EFF = config.D_EFF    # en m^2/s
C_E = config.C_E        # en mol/m^3
D = config.D            # en m
R = config.R            # en m
k = config.k            # s^-1
