# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 20:32:29 2023

@author: flore
"""

import numpy as np
from MEC1315_STL import *

# fonction pour homothétie
# facteur doit être float

def homothetie(V, facteur):
    V = facteur * V
    return V

# fonction pour translation
# déplacement doit être array de format [1,3]

def translation(V, deplacement):
    V = V + deplacement
    return V

# fonction pour répétition circulaire
# nb_rep correspond au nombre de répétitions souhaité, doit être int

def rep_circulaire(F, V, N, nb_rep):
    F_final, V_final, N_final = np.empty([0,3]), np.empty([0,3]), np.empty([0,3]) # création des arrays F, V et N finaux
    nb_vertex = len(V) # on détermine le nombre de vertex de l'objet original
    
    for i in range(nb_rep):
        theta = 2*np.pi/nb_rep*i # on trouve l'angle pour la répétition i
        R = Rz(theta) # matrice de rotation par rapport à l'axe des z
        
        F_i, V_i, N_i = F, V, N # créations de arrays F, V et N pour l'itération i
        
        F_final = np.vstack((F_final, F_i+nb_vertex*i)) # concaténation et ajout de nb_vertex*i sur F_i
        V_final = np.vstack((V_final, V_i.dot(R))) # rotation et concaténation
        N_final = np.vstack((N_final, N_i.dot(R))) # rotation et concaténation
        
    return F_final, V_final, N_final

def fusion(objets): # objets doit être une liste de listes, avec chaque objet à fusionner comme étant une élément-liste tel que [F, V, N]
    F_tot, V_tot, N_tot =np.empty([0,3]),np.empty([0,3]),np.empty([0,3]) # création des arrays F, V et N finaux
    i = 0 # initialisation du compteur
    
    for objet_individuel in objets: # on prend F, V et N de chaque objet à fusionner
        if i == 0:
            F_tot = np.vstack((F_tot, objet_individuel[0])) # ajout du premier objet à la matrice vide F_tot
        else:
            nb_vertex = len(V_tot) # on détermine le nombre de vertex total des objets déjà ''fusionnés''
            F_tot = np.vstack((F_tot, objet_individuel[0]+nb_vertex)) # concaténation et ajout de nb_vertex sur F_tot
            
        V_tot = np.vstack((V_tot, objet_individuel[1])) # concaténation
        N_tot = np.vstack((N_tot, objet_individuel[2])) # concaténation
        i += 1

    return F_tot, V_tot, N_tot

def rep_circulaire2(fichier, nb_rep, x, y, z, Grandissement):
    F_final, V_final, N_final = np.empty([0,3]), np.empty([0,3]), np.empty([0,3]) # création des arrays F, V et N finaux
    F,V,N=LireSTL(fichier)
    nb_vertex = len(V) # on détermine le nombre de vertex de l'objet original
    V=homothetie(V, Grandissement)
    V=translation(V, np.array([x, y, z]))
    for i in range(nb_rep):
        theta = 2*np.pi/nb_rep*i # on trouve l'angle pour la répétition i
        R = Rz(theta) # matrice de rotation par rapport à l'axe des z
        
        F_i, V_i, N_i = F, V, N # créations de arrays F, V et N pour l'itération i
        
        
        F_final = np.vstack((F_final, F_i+nb_vertex*i)) # concaténation et ajout de nb_vertex*i sur F_i
        V_final = np.vstack((V_final, V_i.dot(R))) # Rotation et concaténation
        N_final = np.vstack((N_final, N_i.dot(R))) # Rotation concaténation
        
    return F_final, V_final, N_final

def repéperso(planète_total,planete_répliquer,planète_centrale,grandissement_centrale): 
    f,v,n=planete_répliquer  
    fc,vc,nc=LireSTL(planète_centrale)
    nbre=planète_total
    vc = homothetie(vc, grandissement_centrale)
    objet=[[fc,vc,nc]]
    #grandissement planète secondaire
    m=nbre+1
    gi=1/m
    
    #suite de fibonacci dans les 4 quadrants
    i=0
    j=1
    a=0
    res=np.array([0])
    
    ppmm=np.array([1,1,-1,-1])
    ppmmtot=np.array([1,1,-1,-1])
    pmmp=np.array([1,-1,-1,1])
    pmmptot=np.array([1,-1,-1,1])
    #Matrice pour créer la rotation dans les 4 quadrants
    while len(pmmptot)<nbre: #Pour créer la suite +--+ pour x
        pmmptot=np.hstack([pmmptot,pmmp])
    while len(ppmmtot)<nbre: #Pour créer la suite ++-- pour y
        ppmmtot=np.hstack([ppmmtot,ppmm])
   #Suite de fibonnacci     
    while a<nbre: 
        k=np.array([i+j])
        res=np.hstack([res,k])
        j=i
        i=res[-1]
        a+=1
        
    positionx=np.array([1])
    positiony=np.array([1])
    positionz=np.zeros(nbre)
    for o in range(1,nbre):
        x=res[o+1]*pmmptot[o]
        positionx=np.hstack([positionx,x])
        y=res[o+1]*ppmmtot[o]
        positiony=np.hstack([positiony,y])
    position=np.vstack([positionx, positiony,positionz]).T
    p=position
    #Placement des planètes sur la suite de fibonacci
    for i in range(nbre):
        f1,n1=f,n
        
        v1=v*[i/m]
        v1=v1+110*p[i] #comme la suite est minime comparer à la taille des planètes, grandissement de la translation de 200
        objet.append([f1,v1,n1])
    F_final,V_final,N_final = fusion(objet)
    return F_final,V_final,N_final

#Change dimension d'un objet en appliquant un coefficient a, b et c sur la composante x, y et z
def affinite_vectorielle(F,V,N,a,b,c):
    F1, V1, N1 = F.copy(), V.copy(), N.copy()
    V1[ :,0]=V1[ :,0]*a
    V1[ :,1]=V1[ :,1]*b
    V1[ :,2]=V1[ :,2]*c
    N1=CalculNormal( F1, V1 )
    return F1, V1, N1


def repetition_rectiligne(objet, repetition, espacement):
    F_final, V_final, N_final = np.empty([0,3]), np.empty([0,3]), np.empty([0,3]) # création des arrays F, V et N finaux
    F, V, N = objet[0], objet[1], objet[2]
    nb_vertex = len(V) # on détermine le nombre de vertex de l'objet original
    
    for i  in range (repetition):
        F_i, V_i, N_i = F, V, N # créations de arrays F, V et N pour l'itération i
        V_i[ :,0]=V_i[ :,0]+espacement
        F_final = np.vstack((F_final, F_i+nb_vertex*i)) # concaténation et ajout de nb_vertex*i sur F_i
        V_final = np.vstack((V_final, V_i)) # concaténation
        N_final = np.vstack((N_final, N_i)) # concaténation
        
    return F_final, V_final, N_final

def centrer(objet):
    F, V, N = objet[0], objet[1], objet[2]
    centre_x=(min(V[:,0])+max(V[:,0]))/2
    centre_y=(min(V[:,1])+max(V[:,1]))/2
    min_z=min(V[:,2])
    V=translation(V, np.array([-centre_x,-centre_y,-min_z]))
    return F, V, N

def fonction_drapeau(cylindre,triangle,grandissement):
    #Cylindre et triangle
    fc,vc,nc = LireSTL(cylindre)
    ft,vt,nt= LireSTL(triangle)
    #dimension cylindre
    vc[:,0]=homothetie(vc[:,0], 0.5)
    vc[:,2]=homothetie(vc[:,2], 3)
    vc[:,1]=homothetie(vc[:,1], 0.5)
    #dimension triangle
    vt[:,1]=homothetie(vt[:,1], 6)
    vt[:,2]=homothetie(vt[:,2], 4)
    vt=homothetie(vt, 0.15)

    #Placer le triangle en haut et à gauche du cylindre
    ycmax=max(vc[:,1])
    vt[:,1]=vt[:,1]+ycmax
    zt=max(vt[:,2])-min(vt[:,2])
    ht=max(vc[:,2])-zt
    vt[:,2]=vt[:,2]+ht


    objet=[[fc,vc,nc],[ft,vt,nt]]
    f,v,n=fusion(objet)
    v=homothetie(v, grandissement)
    return f,v,n

def fonction_satellite(cube,x,y,grandissement):
    f0,v0,n0=LireSTL(cube)

    rectangle1=affinite_vectorielle(f0,v0,n0, 0.5, 1.5, 0.2)
    rectangles = repetition_rectiligne(rectangle1, 4, 1)

    f1,v1,n1=centrer(rectangles) #les 4 panneaux supérieur de l'aile
    v1=translation(v1,[0,0.9,0])
    Ailesup= [f1,v1,n1]
    
    f2,v2,n2=centrer(affinite_vectorielle(f0,v0,n0, 5, 0.4, 0.5)) #centre de l'aile
    v2=translation(v2,[0.5,0,-0.2])
    Ailemid= [f2,v2,n2]

    f3,v3,n3=centrer(rectangles) #les 4 panneaux inférieur de l'aile
    v3=translation(v3,[0,-0.9,0])
    Aileinf=[f3,v3,n3]

    f4,v4,n4= fusion((Ailesup,Ailemid,Aileinf)) #une des ailes du satellite
    aile1=[f4,v4,n4] 

    v5=translation(v4.dot(Rz(np.pi)),[7,0,0])  #la deuxième aile du satellite
    aile2=[f4,v5,n4.dot(Rz(np.pi))]

    f6,v6,n6=centrer(affinite_vectorielle(f0,v0,n0, 2, 2.5,2)) #Centre satellite
    v6=translation(v6, [3.5,0,-0.8])
    centre=[f6,v6,n6]

    f7,v7,n7=fusion((aile1,aile2,centre))
    v7=homothetie(v7, grandissement)
    v7[:,0],v7[:,1]=v7[:,0]+x,v7[:,1]+y
    return f7,v7,n7

def minion_drapeau(fichier_drapeau, fichier_minion, grandissement): # fusion de l'objet drapeau et l'objet minion
    f1, v1, n1 = fichier_drapeau
    nv1 = len(v1)
    
    # coïncidence avec l'origine de l'objet 1
    
    centre_x=(min(v1[:,0])+max(v1[:,0]))/2
    centre_y=(min(v1[:,1])+max(v1[:,1]))/2
    min_z=min(v1[:,2])
    v1=v1+np.array([-centre_x,-centre_y,-min_z])
    
    

    v1 = translation(v1, [0,-2,0]) # translation de 2 unités vers les y négatif

    v1 = translation(v1, [0,0,-1]) # translation de 1 unitées vers les z négatifs



    f2, v2, n2 = LireSTL(fichier_minion) # importation de l'objet (2) minion
    
    nv2 = len(v2)
    v2= homothetie(v2, grandissement)
    #coïcidence avec l'origine de l'objet 2

    centre_x=(min(v2[:,0])+max(v2[:,0]))/2
    centre_y=(min(v2[:,1])+max(v2[:,1]))/2
    min_z=min(v2[:,2])
    v2=v2+np.array([-centre_x,-centre_y,-min_z])
    
    
    # fusion des deux objets et création du fichier "minion_drapeau.stl"
    
    minion_drapeau = fusion([[f1,v1,n1],[f2, v2, n2]])
    ff,vf,nf=minion_drapeau
    
    return ff,vf,nf
    
    
def planete_colonisee(fichier_minion_drapeau, fichier_planete): # fusion de l'objet minion_drapeau et de la planète
                                                                # à coloniser

    f1, v1, n1 = fichier_minion_drapeau # importation de l'objet (1) minion_drapeau
    nv1 = len(v1)
    

    v1 = homothetie(v1, 20) # homotéthie de l'objet 1
    

    f2, v2, n2 = LireSTL(fichier_planete) # importation de l'objet (2) planète
    nv2 = len(v2) 
    

    v2 = translation(v2, [0,0,-(max(v2[:,2])-min(v2[:,2]))/2]) # translation vers z négatif de la longueur du rayon<
    
    # fusion des deux objets et création du fichier "planete_colonisee.stl"
    
    planete_colonisee = fusion([[f1,v1,n1],[f2, v2, n2]])
    ff,vf,nf=planete_colonisee
    return ff,vf,nf