# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 14:46:14 2020
Exemple de script permettant de calculer la valeur du coefficient de friction 
d'un sujet en fonction de la force normale qu'il applique. Concretement, on
calcule les valeurs de k et n dans la formule:
    
    mu=k*NF^(n-1)

Ceci permet de calculer la slip force et la marge de securite

Dans ce script, on genere aussi les graphes montrant les differents moments du
glissement detectes ainsi que les graphes montrant mu en fonction de NF.

@author: fschiltz
"""
# Importation des differentes libraries necessaires
import numpy as np
import matplotlib.pyplot as plt


import glm_data_processing as glm
import get_mu_points as gmp
import get_mu_fit as gmf

# Fermeture des figures ouvertes
plt.close('all')

# Chemins d'acces aux fichiers (A MODIFIER)
filepath_strong = '.\S10_CF_004.glm'
filepath_medium = '.\S10_CF_005.glm'
filepath_weak = '.\S10_CF_006.glm'

# Importation des donnes
glm_strong_df = glm.import_data(filepath_strong)
glm_medium_df = glm.import_data(filepath_medium)
glm_weak_df = glm.import_data(filepath_weak)

#%% Filtrage des donnees
freqAcq=800 # Frequence d'acquisition des donnees
freqFiltForces=20; # Frequence de coupure du filtrage des forces (peut etre modifié)

for i in range(21,43):
    glm_strong_df.iloc[:,i]=glm.filter_signal(glm_strong_df.iloc[:,i],freqAcq,freqFiltForces)
    glm_medium_df.iloc[:,i]=glm.filter_signal(glm_medium_df.iloc[:,i],freqAcq,freqFiltForces)
    glm_weak_df.iloc[:,i]=glm.filter_signal(glm_weak_df.iloc[:,i],freqAcq,freqFiltForces)

# Frequence de coupure du filtrage du COP (peut etre modifié)
# COP = Center Of Pressure. Il s'agit de la résultante des forces appliquées 
# par le doigt sur le capteur. Cela donne donc une idee de la position du doigt
# sur celui-ci
freqFiltCOP=40 
glm_strong_df.loc[:,'OPxgr']=glm.filter_signal(glm_strong_df.loc[:,'OPxgr'],freqAcq,freqFiltCOP)
glm_strong_df.loc[:,'OPxgl']=glm.filter_signal(glm_strong_df.loc[:,'OPxgl'],freqAcq,freqFiltCOP)
glm_medium_df.loc[:,'OPxgr']=glm.filter_signal(glm_medium_df.loc[:,'OPxgr'],freqAcq,freqFiltCOP)
glm_medium_df.loc[:,'OPxgl']=glm.filter_signal(glm_medium_df.loc[:,'OPxgl'],freqAcq,freqFiltCOP)
glm_weak_df.loc[:,'OPxgr']=glm.filter_signal(glm_weak_df.loc[:,'OPxgr'],freqAcq,freqFiltCOP)
glm_weak_df.loc[:,'OPxgl']=glm.filter_signal(glm_weak_df.loc[:,'OPxgl'],freqAcq,freqFiltCOP)

#%% Extraction des signaux et mise a zero des forces en soustrayant la valeur 
# moyenne des 500 premieres ms
baseline = range(0,400)
# Normal Force exerted by the thumb
NF_thumb_strong  = glm_strong_df.loc[:,'Fygl']-np.nanmean(glm_strong_df.loc[baseline,'Fygl'])
# Vertical Tangential Force exerted by the thumb
TFx_thumb_strong  = glm_strong_df.loc[:,'Fxgl']-np.nanmean(glm_strong_df.loc[baseline,'Fxgl'])
#Horizontal Tangential Force exerted by the thumb
TFz_thumb_strong  = glm_strong_df.loc[:,'Fzgl']-np.nanmean(glm_strong_df.loc[baseline,'Fzgl'])
#Vertical component of the Center Of Pressure of the thumb
COP_thumb_strong = -glm_strong_df.loc[:,'OPxgl']


# Normal Force exerted by the index
NF_index_strong  = -(glm_strong_df.loc[:,'Fygr']-np.nanmean(glm_strong_df.loc[baseline,'Fygr']))
# Vertical Tangential Force exerted by the index
TFx_index_strong  = glm_strong_df.loc[:,'Fxgr']-np.nanmean(glm_strong_df.loc[baseline,'Fxgr'])
#Horizontal Tangential Force exerted by the index
TFz_index_strong  = glm_strong_df.loc[:,'Fzgr']-np.nanmean(glm_strong_df.loc[baseline,'Fzgr'])
#Vertical component of the Center Of Pressure of the index
COP_index_strong = -glm_strong_df.loc[:,'OPxgr']


# Normal Force exerted by the thumb
NF_thumb_medium  = glm_medium_df.loc[:,'Fygl']-np.nanmean(glm_medium_df.loc[baseline,'Fygl'])
# Vertical Tangential Force exerted by the thumb
TFx_thumb_medium  = glm_medium_df.loc[:,'Fxgl']-np.nanmean(glm_medium_df.loc[baseline,'Fxgl'])
#Horizontal Tangential Force exerted by the thumb
TFz_thumb_medium  = glm_medium_df.loc[:,'Fzgl']-np.nanmean(glm_medium_df.loc[baseline,'Fzgl'])
#Vertical component of the Center Of Pressure of the thumb
COP_thumb_medium = -glm_medium_df.loc[:,'OPxgl']


# Normal Force exerted by the index
NF_index_medium  = -(glm_medium_df.loc[:,'Fygr']-np.nanmean(glm_medium_df.loc[baseline,'Fygr']))
# Vertical Tangential Force exerted by the index
TFx_index_medium  = glm_medium_df.loc[:,'Fxgr']-np.nanmean(glm_medium_df.loc[baseline,'Fxgr'])
#Horizontal Tangential Force exerted by the index
TFz_index_medium  = glm_medium_df.loc[:,'Fzgr']-np.nanmean(glm_medium_df.loc[baseline,'Fzgr'])
#Vertical component of the Center Of Pressure of the index
COP_index_medium = -glm_medium_df.loc[:,'OPxgr']

# Normal Force exerted by the thumb
NF_thumb_weak  = glm_weak_df.loc[:,'Fygl']-np.nanmean(glm_weak_df.loc[baseline,'Fygl'])
# Vertical Tangential Force exerted by the thumb
TFx_thumb_weak  = glm_weak_df.loc[:,'Fxgl']-np.nanmean(glm_weak_df.loc[baseline,'Fxgl'])
#Horizontal Tangential Force exerted by the thumb
TFz_thumb_weak  = glm_weak_df.loc[:,'Fzgl']-np.nanmean(glm_weak_df.loc[baseline,'Fzgl'])
#Vertical component of the Center Of Pressure of the thumb
COP_thumb_weak = -glm_weak_df.loc[:,'OPxgl']


# Normal Force exerted by the index
NF_index_weak  = -(glm_weak_df.loc[:,'Fygr']-np.nanmean(glm_weak_df.loc[baseline,'Fygr']))
# Vertical Tangential Force exerted by the index
TFx_index_weak  = glm_weak_df.loc[:,'Fxgr']-np.nanmean(glm_weak_df.loc[baseline,'Fxgr'])
#Horizontal Tangential Force exerted by the index
TFz_index_weak  = glm_weak_df.loc[:,'Fzgr']-np.nanmean(glm_weak_df.loc[baseline,'Fzgr'])
#Vertical component of the Center Of Pressure of the index
COP_index_weak = -glm_weak_df.loc[:,'OPxgr']

#%% Calcul du coefficient de friction au moment du glissement. Voir fonction 
# get_mu_points et papier "Barrea et al. 2016" pour plus de détails
mu_thumb_strong,slip_indexes_thumb_strong,start_search_zones_thumb_strong,end_search_zones_thumb_strong = \
gmp.get_mu_points(COP_thumb_strong,TFz_thumb_strong,TFx_thumb_strong,NF_thumb_strong)

mu_index_strong,slip_indexes_index_strong,start_search_zones_index_strong,end_search_zones_index_strong = \
gmp.get_mu_points(COP_index_strong,TFz_index_strong,TFx_index_strong,NF_index_strong)


mu_thumb_medium,slip_indexes_thumb_medium,start_search_zones_thumb_medium,end_search_zones_thumb_medium = \
gmp.get_mu_points(COP_thumb_medium,TFz_thumb_medium,TFx_thumb_medium,NF_thumb_medium)

mu_index_medium,slip_indexes_index_medium,start_search_zones_index_medium,end_search_zones_index_medium = \
gmp.get_mu_points(COP_index_medium,TFz_index_medium,TFx_index_medium,NF_index_medium)

mu_thumb_weak,slip_indexes_thumb_weak,start_search_zones_thumb_weak,end_search_zones_thumb_weak = \
gmp.get_mu_points(COP_thumb_weak,TFz_thumb_weak,TFx_thumb_weak,NF_thumb_weak)

mu_index_weak,slip_indexes_index_weak,start_search_zones_index_weak,end_search_zones_index_weak = \
gmp.get_mu_points(COP_index_weak,TFz_index_weak,TFx_index_weak,NF_index_weak)


#%% Calcul des valeurs de k et n pour l'index
mu_vector_index=np.append(mu_index_strong,np.append(mu_index_medium,mu_index_weak))
NF_vector_index=np.append(NF_index_strong[slip_indexes_index_strong],\
                                         np.append(NF_index_medium[slip_indexes_index_medium],\
                                         NF_index_weak[slip_indexes_index_weak]))
k_index,n_index=gmf.get_mu_fit(mu_vector_index,NF_vector_index)

#%% Calcul des valeurs de k et n pour le pouce
mu_vector_thumb=np.append(mu_thumb_strong,np.append(mu_thumb_medium,mu_thumb_weak))
NF_vector_thumb=np.append(NF_thumb_strong[slip_indexes_thumb_strong],\
                                         np.append(NF_thumb_medium[slip_indexes_thumb_medium],\
                                         NF_thumb_weak[slip_indexes_thumb_weak]))
k_thumb,n_thumb=gmf.get_mu_fit(mu_vector_thumb,NF_vector_thumb)

# Les variables k_index, k_thumb, n_index et n_thumb sont ce que l'on cherche
# a obtenir avec ce script!!! Elles permettent de calculer la slip force et
# donc la marge de securite






# Generating some beautiful graphs to check that everything went well

#%% Graphe de détection du moment du glissement pour l'essai "index strong"
siz= len(start_search_zones_index_strong)
time = glm_strong_df.loc[:,'time'].to_numpy()
fig = plt.figure(figsize = [15,7])
ax  = fig.subplots(3,1)
ax[0].set_title("Index strong", fontsize=14, fontweight="bold")
ax[0].plot(time, COP_index_strong*1000)
ax[0].set_ylabel("COP [mm]", fontsize=13)

ax[1].plot(time,TFx_index_strong, label="TFv")
ax[1].plot(time,NF_index_strong, label="NFv")
ax[1].legend(fontsize=12)
ax[1].set_ylabel("Forces [N]", fontsize=13)

ax[2].plot(time,np.divide(np.hypot(TFx_index_strong,TFz_index_strong),NF_index_strong))
ax[2].plot(time[slip_indexes_index_strong],mu_index_strong,marker='.',linestyle='')
ax[2].set_ylim([0,1.5])
ax[2].set_ylabel("TF/NF [-]", fontsize=13)
ax[2].set_xlabel("Time [s]", fontsize=13)

for i in range(0,siz):
    rect0=plt.Rectangle((time[start_search_zones_index_strong[i]],ax[0].get_ylim()[0]),\
                       time[end_search_zones_index_strong[i]-start_search_zones_index_strong[i]],\
                       ax[0].get_ylim()[1]-ax[0].get_ylim()[0],color='k',alpha=0.3)
    rect1=plt.Rectangle((time[start_search_zones_index_strong[i]],ax[1].get_ylim()[0]),\
                       time[end_search_zones_index_strong[i]-start_search_zones_index_strong[i]],\
                       ax[1].get_ylim()[1]-ax[1].get_ylim()[0],color='k',alpha=0.3)
    rect2=plt.Rectangle((time[start_search_zones_index_strong[i]],ax[2].get_ylim()[0]),\
                       time[end_search_zones_index_strong[i]-start_search_zones_index_strong[i]],\
                       ax[2].get_ylim()[1]-ax[2].get_ylim()[0],color='k',alpha=0.3)
    ax[0].add_patch(rect0)
    ax[1].add_patch(rect1)
    ax[2].add_patch(rect2)
    
#%% Graphe de détection du moment du glissement pour l'essai "index medium"
siz= len(start_search_zones_index_medium)
time = glm_medium_df.loc[:,'time'].to_numpy()
fig = plt.figure(figsize = [15,7])
ax  = fig.subplots(3,1)
ax[0].set_title("Index medium", fontsize=14, fontweight="bold")
ax[0].plot(time, COP_index_medium*1000)
ax[0].set_ylabel("COP [mm]", fontsize=13)

ax[1].plot(time,TFx_index_medium, label="TFv")
ax[1].plot(time,NF_index_medium, label="NFv")
ax[1].legend(fontsize=12)
ax[1].set_ylabel("Forces [N]", fontsize=13)

ax[2].plot(time,np.divide(np.hypot(TFx_index_medium,TFz_index_medium),NF_index_medium))
ax[2].plot(time[slip_indexes_index_medium],mu_index_medium,marker='.',linestyle='')
ax[2].set_ylim([0,1.5])
ax[2].set_ylabel("TF/NF [-]", fontsize=13)
ax[2].set_xlabel("Time [s]", fontsize=13)

for i in range(0,siz):
    rect0=plt.Rectangle((time[start_search_zones_index_medium[i]],ax[0].get_ylim()[0]),\
                       time[end_search_zones_index_medium[i]-start_search_zones_index_medium[i]],\
                       ax[0].get_ylim()[1]-ax[0].get_ylim()[0],color='k',alpha=0.3)
    rect1=plt.Rectangle((time[start_search_zones_index_medium[i]],ax[1].get_ylim()[0]),\
                       time[end_search_zones_index_medium[i]-start_search_zones_index_medium[i]],\
                       ax[1].get_ylim()[1]-ax[1].get_ylim()[0],color='k',alpha=0.3)
    rect2=plt.Rectangle((time[start_search_zones_index_medium[i]],ax[2].get_ylim()[0]),\
                       time[end_search_zones_index_medium[i]-start_search_zones_index_medium[i]],\
                       ax[2].get_ylim()[1]-ax[2].get_ylim()[0],color='k',alpha=0.3)
    ax[0].add_patch(rect0)
    ax[1].add_patch(rect1)
    ax[2].add_patch(rect2)

#%% Graphe de détection du moment du glissement pour l'essai "index weak"
siz= len(start_search_zones_index_weak)
time = glm_weak_df.loc[:,'time'].to_numpy()
fig = plt.figure(figsize = [15,7])
ax  = fig.subplots(3,1)
ax[0].set_title("Index weak", fontsize=14, fontweight="bold")
ax[0].plot(time, COP_index_weak*1000)
ax[0].set_ylabel("COP [mm]", fontsize=13)

ax[1].plot(time,TFx_index_weak, label="TFv")
ax[1].plot(time,NF_index_weak, label="NFv")
ax[1].legend(fontsize=12)
ax[1].set_ylabel("Forces [N]", fontsize=13)

ax[2].plot(time,np.divide(np.hypot(TFx_index_weak,TFz_index_weak),NF_index_weak))
ax[2].plot(time[slip_indexes_index_weak],mu_index_weak,marker='.',linestyle='')
ax[2].set_ylim([0,3])
ax[2].set_ylabel("TF/NF [-]", fontsize=13)
ax[2].set_xlabel("Time [s]", fontsize=13)

for i in range(0,siz):
    rect0=plt.Rectangle((time[start_search_zones_index_weak[i]],ax[0].get_ylim()[0]),\
                       time[end_search_zones_index_weak[i]-start_search_zones_index_weak[i]],\
                       ax[0].get_ylim()[1]-ax[0].get_ylim()[0],color='k',alpha=0.3)
    rect1=plt.Rectangle((time[start_search_zones_index_weak[i]],ax[1].get_ylim()[0]),\
                       time[end_search_zones_index_weak[i]-start_search_zones_index_weak[i]],\
                       ax[1].get_ylim()[1]-ax[1].get_ylim()[0],color='k',alpha=0.3)
    rect2=plt.Rectangle((time[start_search_zones_index_weak[i]],ax[2].get_ylim()[0]),\
                       time[end_search_zones_index_weak[i]-start_search_zones_index_weak[i]],\
                       ax[2].get_ylim()[1]-ax[2].get_ylim()[0],color='k',alpha=0.3)
    ax[0].add_patch(rect0)
    ax[1].add_patch(rect1)
    ax[2].add_patch(rect2)

#%% Graphe de détection du moment du glissement pour l'essai "thumb strong"
siz= len(start_search_zones_thumb_strong)
time = glm_strong_df.loc[:,'time'].to_numpy()
fig = plt.figure(figsize = [15,7])
ax  = fig.subplots(3,1)
ax[0].set_title("Thumb strong", fontsize=14, fontweight="bold")
ax[0].plot(time, COP_thumb_strong*1000)
ax[0].set_ylabel("COP [mm]", fontsize=13)

ax[1].plot(time,TFx_thumb_strong, label="TFv")
ax[1].plot(time,NF_thumb_strong, label="NFv")
ax[1].legend(fontsize=12)
ax[1].set_ylabel("Forces [N]", fontsize=13)

ax[2].plot(time,np.divide(np.hypot(TFx_thumb_strong,TFz_thumb_strong),NF_thumb_strong))
ax[2].plot(time[slip_indexes_thumb_strong],mu_thumb_strong,marker='.',linestyle='')
ax[2].set_ylim([0,1.5])
ax[2].set_ylabel("TF/NF [-]", fontsize=13)
ax[2].set_xlabel("Time [s]", fontsize=13)

for i in range(0,siz):
    rect0=plt.Rectangle((time[start_search_zones_thumb_strong[i]],ax[0].get_ylim()[0]),\
                       time[end_search_zones_thumb_strong[i]-start_search_zones_thumb_strong[i]],\
                       ax[0].get_ylim()[1]-ax[0].get_ylim()[0],color='k',alpha=0.3)
    rect1=plt.Rectangle((time[start_search_zones_thumb_strong[i]],ax[1].get_ylim()[0]),\
                       time[end_search_zones_thumb_strong[i]-start_search_zones_thumb_strong[i]],\
                       ax[1].get_ylim()[1]-ax[1].get_ylim()[0],color='k',alpha=0.3)
    rect2=plt.Rectangle((time[start_search_zones_thumb_strong[i]],ax[2].get_ylim()[0]),\
                       time[end_search_zones_thumb_strong[i]-start_search_zones_thumb_strong[i]],\
                       ax[2].get_ylim()[1]-ax[2].get_ylim()[0],color='k',alpha=0.3)
    ax[0].add_patch(rect0)
    ax[1].add_patch(rect1)
    ax[2].add_patch(rect2)

#%% Graphe de détection du moment du glissement pour l'essai "thumb medium"
siz= len(start_search_zones_thumb_medium)
time = glm_medium_df.loc[:,'time'].to_numpy()
fig = plt.figure(figsize = [15,7])
ax  = fig.subplots(3,1)
ax[0].set_title("Thumb medium", fontsize=14, fontweight="bold")
ax[0].plot(time, COP_thumb_medium*1000)
ax[0].set_ylabel("COP [mm]", fontsize=13)

ax[1].plot(time,TFx_thumb_medium, label="TFv")
ax[1].plot(time,NF_thumb_medium, label="NFv")
ax[1].legend(fontsize=12)
ax[1].set_ylabel("Forces [N]", fontsize=13)

ax[2].plot(time,np.divide(np.hypot(TFx_thumb_medium,TFz_thumb_medium),NF_thumb_medium))
ax[2].plot(time[slip_indexes_thumb_medium],mu_thumb_medium,marker='.',linestyle='')
ax[2].set_ylim([0,1.5])
ax[2].set_ylabel("TF/NF [-]", fontsize=13)
ax[2].set_xlabel("Time [s]", fontsize=13)

for i in range(0,siz):
    rect0=plt.Rectangle((time[start_search_zones_thumb_medium[i]],ax[0].get_ylim()[0]),\
                       time[end_search_zones_thumb_medium[i]-start_search_zones_thumb_medium[i]],\
                       ax[0].get_ylim()[1]-ax[0].get_ylim()[0],color='k',alpha=0.3)
    rect1=plt.Rectangle((time[start_search_zones_thumb_medium[i]],ax[1].get_ylim()[0]),\
                       time[end_search_zones_thumb_medium[i]-start_search_zones_thumb_medium[i]],\
                       ax[1].get_ylim()[1]-ax[1].get_ylim()[0],color='k',alpha=0.3)
    rect2=plt.Rectangle((time[start_search_zones_thumb_medium[i]],ax[2].get_ylim()[0]),\
                       time[end_search_zones_thumb_medium[i]-start_search_zones_thumb_medium[i]],\
                       ax[2].get_ylim()[1]-ax[2].get_ylim()[0],color='k',alpha=0.3)
    ax[0].add_patch(rect0)
    ax[1].add_patch(rect1)
    ax[2].add_patch(rect2)

#%% Graphe de détection du moment du glissement pour l'essai "thumb weak"
siz= len(start_search_zones_thumb_weak)
time = glm_weak_df.loc[:,'time'].to_numpy()
fig = plt.figure(figsize = [15,7])
ax  = fig.subplots(3,1)
ax[0].set_title("Thumb weak", fontsize=14, fontweight="bold")
ax[0].plot(time, COP_thumb_weak*1000)
ax[0].set_ylabel("COP [mm]", fontsize=13)

ax[1].plot(time,TFx_thumb_weak, label="TFv")
ax[1].plot(time,NF_thumb_weak, label="NFv")
ax[1].legend(fontsize=12)
ax[1].set_ylabel("Forces [N]", fontsize=13)

ax[2].plot(time,np.divide(np.hypot(TFx_thumb_weak,TFz_thumb_weak),NF_thumb_weak))
ax[2].plot(time[slip_indexes_thumb_weak],mu_thumb_weak,marker='.',linestyle='')
ax[2].set_ylim([0,3])
ax[2].set_ylabel("TF/NF [-]", fontsize=13)
ax[2].set_xlabel("Time [s]", fontsize=13)

for i in range(0,siz):
    rect0=plt.Rectangle((time[start_search_zones_thumb_weak[i]],ax[0].get_ylim()[0]),\
                       time[end_search_zones_thumb_weak[i]-start_search_zones_thumb_weak[i]],\
                       ax[0].get_ylim()[1]-ax[0].get_ylim()[0],color='k',alpha=0.3)
    rect1=plt.Rectangle((time[start_search_zones_thumb_weak[i]],ax[1].get_ylim()[0]),\
                       time[end_search_zones_thumb_weak[i]-start_search_zones_thumb_weak[i]],\
                       ax[1].get_ylim()[1]-ax[1].get_ylim()[0],color='k',alpha=0.3)
    rect2=plt.Rectangle((time[start_search_zones_thumb_weak[i]],ax[2].get_ylim()[0]),\
                       time[end_search_zones_thumb_weak[i]-start_search_zones_thumb_weak[i]],\
                       ax[2].get_ylim()[1]-ax[2].get_ylim()[0],color='k',alpha=0.3)
    ax[0].add_patch(rect0)
    ax[1].add_patch(rect1)
    ax[2].add_patch(rect2)


#%% Figure finale pour l'index
x=np.arange(0.2,30,0.02).tolist()
fig = plt.figure(figsize = [15,7])
plt.plot(x,k_index*(x**(n_index-1)))
plt.plot(NF_vector_index,mu_vector_index,linestyle='',marker='.')
plt.ylim([0,3])
plt.xlim([0,30])
plt.title('Coefficient of friction index')
plt.xlabel('Normal Force [N]')
plt.ylabel('Static Friction [-]')

#%% Figure finale pour le pouce
x=np.arange(0.2,30,0.02).tolist()
fig = plt.figure(figsize = [15,7])
plt.plot(x,k_thumb*(x**(n_thumb-1)))
plt.plot(NF_vector_thumb,mu_vector_thumb,linestyle='',marker='.')
plt.ylim([0,3])
plt.xlim([0,30])
plt.title('Coefficient of friction thumb')
plt.xlabel('Normal Force [N]')
plt.ylabel('Static Friction [-]')

#%% Impression de la valeur des variables dans la console
print("Index: the value of k is %f and n is %f" %(k_index,n_index))
print("Thumb: the value of k is %f and n is %f" %(k_thumb,n_thumb))














