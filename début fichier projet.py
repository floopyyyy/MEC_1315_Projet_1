# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 20:11:38 2023

@author: noahd
"""
"test"

j
import numpy as np
from MEC1315_STL import *
from mes_fonctions_complets import *
#créer fichier pur mini_planete


#Anneau
Anneau_p=rep_circulaire2('mini_planete.stl', 120, 135, 0, 0, 1/50) #Répétition circulaire
Anneau_s1=rep_circulaire2('mini_planete.stl', 120, 120, 0, 10, 1/50)
Anneau_s2=rep_circulaire2('mini_planete.stl', 120, 120, 0, -10, 1/50)



#Planète colonisée
drapeau=fonction_drapeau('Cylindre.stl', 'Triangle.stl', 1)
drapeau_minion=minion_drapeau(drapeau, 'minion.stl')
planete_col=planete_colonisee(drapeau_minion, 'planete.stl')
#Satellite
satellite=fonction_satellite('cube.stl', -300, 300,10) #Répétiton lnéaire dans la fonction

#Galaxy 
galaxy=repéperso(6, planete_col, 'planete.stl', 2/3) #Répétition selon la suite de fibonacci

#fusion et export du fichier STL
scene=fusion((Anneau_p,Anneau_s1,Anneau_s2,galaxy,satellite))
f,v,n=scene[0],scene[1],scene[2]
EcrireSTLASCII('scene_final_XX.stl', f, v, n)
