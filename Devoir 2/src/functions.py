"""
Fichier de calcul des fonction analytique et des concentrations avec la methode des differences finies
"""


#############################
#Importation de bibliothèques
#############################
import sympy as sp
import numpy as np
import config


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
