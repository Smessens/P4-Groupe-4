#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 16:34:02 2021

@author: sophiealicandro
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal 
import os
import glm_data_processing as glm
import derive as der
#from script_coda  import get_segmentations

maxU = [0,0,0,0,0,0,0,0]
maxD = [0,0,0,0,0,0,0,0]

subjects = ["SophieEBR","SimonEBR","JulEBR","JulianEBR"] #Names of subjects
ntrials = 2 #Number of trials for each subject

# Double for-loop that runs thrgough all subjects and trials
subject_number=0;
for s in subjects:
    for trial in range(1,ntrials+1): 
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)
        glm_df = glm.import_data(glm_path)
        baseline = range(0,400)        
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        time  = glm_df.loc[:,'time'].to_numpy()      
        GF    = glm_df.loc[:,'GF'].to_numpy()        
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        segmentations=np.array([])
        mB0= glm_df['Metronome_b0'].to_numpy()
        mB1= glm_df['Metronome_b1'].to_numpy()
        mB2= glm_df['Metronome_b2'].to_numpy()        
        start=4000
        end=30959
        mB0=mB0[start:end]
        mB1=mB1[start:end]
        mB2=mB2[start:end]
        bip=False
        for i in range( len (mB0)):
            
            if((mB0[i]+ mB1[i]+mB2[i])!=3 and bip==False ):
                segmentations=np.append(segmentations,int(i+start))
                bip=True
            elif(bip and (mB0[i]+mB1[i]+mB2[i])==3 ):
                bip=False
            
        segmentations=segmentations.astype(int)       
        le=int(segmentations[18]-segmentations[17])        
        GFU=np.full(le,0).astype(float)
       
        for i in range (le):    
            GFU[i] = GF[segmentations[17]+i]
            if GFU[i] >= maxU[0]:
                maxU[0] = GFU[i]
        #print(maxU[0])  
        maxU[0]=0
        led=int(segmentations[19]-segmentations[18])
        GFD=np.full(led,0).astype(float)
        
        for i in range(led):
            GFD[i] = GF[segmentations[18]+i]
            if GFD[i] >= maxD[0]:
                maxD[0] = GFD[i]
        #print(maxD[0]) 
        maxD[0]=0


spreadSophieU = [0.53079294, 1.239255547, 2.80378679, 2.946274437]
spreadSophieUR = [1.586643271, 1.36798452, 2.166724393, 2.391039002]
spreadSophieD = [0.700555477, 0.901308539, 3.146430137, 1.41750693]
spreadSophieDR = [1.733071297, 0.977250338, 1.659824311, 0.88180272]

spreadSimonU = [2.156246467, 1.303019104, 2.784038581, 1.430901948]
spreadSimonUR = [0.988143196, 3.953914756, 2.562684066, 2.562684066]
spreadSimonD = [1.466089652, 1.654844374, 1.367962674, 1.148067443]
spreadSimonDR = [2.165007652, 1.782591901, 3.00158859, 1.946509103]

spreadJulienU = [0.839565917, 0.524861828, 1.142932487, 1.141985652]
spreadJulienUR = [0.801779361, 0.398843683, 1.77969223, 1.388055059]
spreadJulienD = [0.941637661, 0.673232949, 1.455420717, 0.881883532]
spreadJulienDR = [1.254305656, 0.979264776, 1.549384961, 1.257301408]

spreadJulianU = [0.597881792, 0.966944216, 1.375373665, 3.084130326]
spreadJulianUR = [3.486713254, 1.28222396, 1.910181688, 2.327042903]
spreadJulianD = [1.001954256, 2.043213794, 2.67675927, 2.261912143]
spreadJulianDR = [2.043213794, 1.075565718, 3.247295673, 1.711930341]


centerSophieU = [7.159153822, 7.017903124, 13.43630521, 15.8920525 ]#SECV, SEF, EH, EB
centerSophieUR =  [9.264205071, 7.848004883, 11.01516796, 9.746792823 ]#SECVR, SEFR, EHR, EBR
centerSophieD = [7.302516614, 7.399769888, 12.51337835, 10.44276115]
centerSophieDR = [10.8719462, 7.13922512, 10.09312117, 9.529586823]

centerSimonU = [10.00257279, 9.274698496, 12.42718742, 9.561645119 ]
centerSimonUR = [12.97598207, 13.28539429, 13.98290533, 16.31885375 ]
centerSimonD = [9.29512317, 10.6822477, 10.13287013, 8.790734636]
centerSimonDR = [15.08755523, 11.21959917, 14.02226815, 15.51076592]

centerJulienU = [9.38408849, 8.338520413, 12.01069969, 9.323454912 ]
centerJulienUR = [10.19930461, 9.103526781, 13.72106279, 13.05104668 ]
centerJulienD = [9.510305298, 9.791070116, 11.56234564, 8.300168739 ]
centerJulienDR = [8.63643242, 10.40932307, 13.24505838, 10.45902418 ]

centerJulianU = [8.84898894, 8.702948364, 20.96475031, 14.76494878 ]
centerJulianUR = [11.50972456, 10.25830437, 15.36982663, 15.93129996 ]
centerJulianD = [9.14097365, 7.183287711, 15.9139846, 9.958285369]
centerJulianDR = [12.67573746, 9.046052123, 13.60225087, 16.4434393]


flier_highSophieU = [7.967650228, 11.02393288, 21.24259111, 20.61160437]
flier_lowSophieU = [6.267279602, 5.722415252, 8.802250623, 10.15794906]
flier_highSophieUR = [14.02830021, 11.9345425, 16.40576611, 13.9608589]
flier_lowSophieUR = [6.848568785, 6.366428126, 8.131020634, 6.668295734]
flier_highSophieD = [8.882098716, 10.45790304, 19.61371061, 13.25934135]
flier_lowSophieD = [6.258077458, 5.745974539, 8.371519732, 7.665358155]
flier_highSophieDR = [15.48226258, 9.712270891, 14.8799095, 11.18099932]
flier_lowSophieDR = [8.61440215, 5.586936934, 8.028695722, 8.100115663]

flier_highSimonU = [14.27380933, 12.24205933, 18.02608845, 12.67150423]
flier_lowSimonU = [7.03081822, 6.850115923, 8.09503243, 7.575295068]
flier_highSimonUR = [14.68679395, 23.46270809, 21.01294535, 21.01294535]
flier_lowSimonUR = [11.50154877, 9.321525709, 10.68975686, 10.68975686]
flier_highSimonD = [13.04140191, 13.43019404, 13.87556592, 11.07301896]
flier_lowSimonD = [7.138172909, 7.848814482, 8.719686897, 7.404327793]
flier_highSimonDR = [20.1748258, 15.63384389, 21.1113327, 19.72632863]
flier_lowSimonDR = [11.61878346, 9.068526968, 10.27014994, 12.0010105]

flier_highJulienU = [10.62198726, 9.318039857, 14.12513365, 11.40087239]
flier_lowJulienU = [8.059736228, 7.385777046, 10.23531983, 8.100215086]
flier_highJulienUR = [11.73030093, 10.06899764, 16.97238526, 16.01207396]
flier_lowJulienUR = [8.397855667, 8.563460991, 10.98098592, 11.06504503]
flier_highJulienD = [11.02331646, 11.01141405, 15.6526829, 10.56651735]
flier_lowJulienD= [7.817665865, 8.56324955, 9.377041673, 6.77000428]
flier_highJulienDR = [11.21060569, 12.01861814, 15.66126716, 14.44679359]
flier_lowJulienDR = [6.814497959, 9.018150049, 10.1194806, 8.968659477]

flier_highJulianU = [10.23439006, 11.19132868, 23.92225808, 19.89803538]
flier_lowJulianU = [8.152478372, 7.021132293, 18.84843661, 9.823433404]
flier_highJulianUR = [17.76364819, 12.5997234, 19.26043431, 19.08611132]
flier_lowJulianUR = [7.750268779, 8.607634759, 12.63633997, 12.1355679]
flier_highJulianD = [12.24889449, 8.69945277, 23.51005309, 14.84207203]
flier_lowJulianD = [7.901882456, 6.001969445, 11.69494008, 6.718504803]
flier_highJulianDR = [16.74561224, 11.39866556, 20.43193977, 19.64340733]
flier_lowJulianDR = [7.457422379, 7.19062089, 8.584751729, 12.61254389]

"""

SECVU=[centerJulianU[0],flier_highJulianU[0],flier_lowJulianU[0]]
SECVRU=[centerJulianUR[0],flier_highJulianUR[0],flier_lowJulianUR[0]]
SEFU=[centerJulianU[1],flier_highJulianU[1],flier_lowJulianU[1]]
SEFRU=[centerJulianUR[1],flier_highJulianUR[1],flier_lowJulianUR[1]]
EHU=[centerJulianU[2],flier_highJulianU[2],flier_lowJulianU[2]]
EHRU=[centerJulianUR[2],flier_highJulianUR[2],flier_lowJulianUR[2]]
EBU =[centerJulianU[3],flier_highJulianU[3],flier_lowJulianU[3]]
EBRU=[centerJulianUR[3],flier_highJulianUR[3],flier_lowJulianUR[3]]

SECVD=[centerJulianD[0],flier_highJulianD[0],flier_lowJulianD[0]]
SECVRD=[centerJulianDR[0],flier_highJulianDR[0],flier_lowJulianDR[0]]
SEFD=[centerJulianD[1],flier_highJulianD[1],flier_lowJulianD[1]]
SEFRD=[centerJulianDR[1],flier_highJulianDR[1],flier_lowJulianDR[1]]
EHD=[centerJulianD[2],flier_highJulianD[2],flier_lowJulianD[2]]
EHRD=[centerJulianDR[2],flier_highJulianDR[2],flier_lowJulianDR[2]]
EBD =[centerJulianD[3],flier_highJulianD[3],flier_lowJulianD[3]]
EBRD=[centerJulianDR[3],flier_highJulianDR[3],flier_lowJulianDR[3]]
"""
"""
SECVU=[centerJulienU[0],flier_highJulienU[0],flier_lowJulienU[0]]
SECVRU=[centerJulienUR[0],flier_highJulienUR[0],flier_lowJulienUR[0]]
SEFU=[centerJulienU[1],flier_highJulienU[1],flier_lowJulienU[1]]
SEFRU=[centerJulienUR[1],flier_highJulienUR[1],flier_lowJulienUR[1]]
EHU=[centerJulienU[2],flier_highJulienU[2],flier_lowJulienU[2]]
EHRU=[centerJulienUR[2],flier_highJulienUR[2],flier_lowJulienUR[2]]
EBU =[centerJulienU[3],flier_highJulienU[3],flier_lowJulienU[3]]
EBRU=[centerJulienUR[3],flier_highJulienUR[3],flier_lowJulienUR[3]]

SECVD=[centerJulienD[0],flier_highJulienD[0],flier_lowJulienD[0]]
SECVRD=[centerJulienDR[0],flier_highJulienDR[0],flier_lowJulienDR[0]]
SEFD=[centerJulienD[1],flier_highJulienD[1],flier_lowJulienD[1]]
SEFRD=[centerJulienDR[1],flier_highJulienDR[1],flier_lowJulienDR[1]]
EHD=[centerJulienD[2],flier_highJulienD[2],flier_lowJulienD[2]]
EHRD=[centerJulienDR[2],flier_highJulienDR[2],flier_lowJulienDR[2]]
EBD =[centerJulienD[3],flier_highJulienD[3],flier_lowJulienD[3]]
EBRD=[centerJulienDR[3],flier_highJulienDR[3],flier_lowJulienDR[3]]
"""
SECVU=[centerJulienU[0],flier_highJulienU[0],flier_lowJulienU[0]]
SECVRU=[centerJulienUR[0],flier_highJulienUR[0],flier_lowJulienUR[0]]
SEFU=[centerJulienU[1],flier_highJulienU[1],flier_lowJulienU[1]]
SEFRU=[centerJulienUR[1],flier_highJulienUR[1],flier_lowJulienUR[1]]
EHU=[centerJulienU[2],flier_highJulienU[2],flier_lowJulienU[2]]
EHRU=[centerJulienUR[2],flier_highJulienUR[2],flier_lowJulienUR[2]]
EBU =[centerJulienU[3],flier_highJulienU[3],flier_lowJulienU[3]]
EBRU=[centerJulienUR[3],flier_highJulienUR[3],flier_lowJulienUR[3]]

SECVD=[centerJulienD[0],flier_highJulienD[0],flier_lowJulienD[0]]
SECVRD=[centerJulienDR[0],flier_highJulienDR[0],flier_lowJulienDR[0]]
SEFD=[centerJulienD[1],flier_highJulienD[1],flier_lowJulienD[1]]
SEFRD=[centerJulienDR[1],flier_highJulienDR[1],flier_lowJulienDR[1]]
EHD=[centerJulienD[2],flier_highJulienD[2],flier_lowJulienD[2]]
EHRD=[centerJulienDR[2],flier_highJulienDR[2],flier_lowJulienDR[2]]
EBD =[centerJulienD[3],flier_highJulienD[3],flier_lowJulienD[3]]
EBRD=[centerJulienDR[3],flier_highJulienDR[3],flier_lowJulienDR[3]]

dataU = [SECVU,SEFU,EHU,EBU]
dataUR = [SECVRU,SEFRU,EHRU,EBRU]
dataD = [SECVD,SEFD,EHD,EBD]
dataDR = [SECVRD,SEFRD,EHRD,EBRD]

labels=['SECV','SEF','EH','EB']
labels2 = ['SECVR','SEFR','EHR','EBR']
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=2, figsize=(9, 9))

bplot1 = ax1[0].boxplot(dataU,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=labels)  # will be used to label x-ticks
ax1[0].set_title("Maxima de la Grip force à l'endroit-UP")
bplot2 = ax1[1].boxplot(dataUR,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=labels2)  # will be used to label x-ticks
ax1[1].set_title("Maxima de la Grip Force à l'envers-UP")
bplot3 = ax2[0].boxplot(dataD,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=labels)  # will be used to label x-ticks
ax2[0].set_title("Maxima de la Grip Force à l'endroit-DOWN")
bplot4 = ax2[1].boxplot(dataDR,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=labels2)  # will be used to label x-ticks
ax2[1].set_title("Maxima de la Grip Force à l'envers-DOWN")


colors = ['pink', 'lightblue', 'lightgreen','yellow']
for bplot in (bplot1, bplot2, bplot3, bplot4):
    for patch, color in zip(bplot['boxes'], colors):
       patch.set_facecolor(color)

