#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 11:06:53 2021

@author: jardaille
"""

# -*- coding: utf-8 -*-
"""
Example script for processing and plotting GLM data
The script uses the data of an oscillation task performed with the 
manipulandum (file TEST_DATA.glm)

Created on Wed Jan 29 11:16:06 2020

@author: opsomerl & fschiltz
"""
#%% Importation des librairies necessaires
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os

import glm_data_processing as glm
import derive as der
from script_coda  import get_segmentations

# Fermeture des figures ouvertes
plt.close('all')
# Double for-loop that runs thrgough all subjects and trials
#subjects = ["JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF"] #Names of subjects
#subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF"] #Names of subjects
#Notes  not existant: EB2, EHR
#subjects = ["SophieEB","SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF"] #Names of subjects
subjects = ["JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
#subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","SophieEB","SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF","JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
ntrials = 2 #Number of trials for each subject
subject_number=0;
for s in subjects:
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        
        #%% CUTTING THE TASK INTO SEGMENTS (your first task)
        pk = signal.find_peaks(accX,  prominence=5, width = (100,1000))
        ipk = pk[0]
        cycle_starts = ipk[:-1]
        cycle_ends = ipk[1:]-1
        
        #%% Compute derivative of LF
        dGF=der.derive(GF,800)
        dGF=glm.filter_signal(dGF,   fs = freqAcq, fc = 10)
        
        name ="%s_00%d.glm" % (s,trial)
        
        #%% segmentations
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
           
           
     
    #  Plot de GF en fonction de LF pour chaque mode avec les données de tout le monde            
       
#        subjectsSECV = ["JulSECV","SimonSECV","SophieSECV","JulianSECV"] #Names of subjects
#        subjectsSECVR = ["julSECVR","SimonSECVR","SophieSECVR","JulianSECVR"] #Names of subjects
#        subjectsSEF = ["julSEF","SimonSEF","SophieSEF","JulianSEF"]
#        subjectsSEFR = ["julSEFR","SimonSEFR","SophieSEFR","JulianSEFR"]
#        subjectsEH = ["JULEH","SimonEH","SophieEH","JulianEH"]
#        subjectsEHR = ["JULEHR","SimonEHR","SophieEHR","JulianEHR"]
#        subjectsEB = ["JULEB","SimonEB","SophieEB","JulianEB"]
#        subjectsEBR = ["JULEBR","SimonEBR","SophieEBR","JulianEBR"]
#        
subjects = ["JulSECV","SimonSECV","SophieSECV","JulianSECV"]#Names of subjects
ntrials = 2 
#kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['b', 'g', 'r', 'c','b', 'g', 'r', 'c']
lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
c = 0
d = 0
k = 0
for s in subjects :
    
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        fig = plt.figure(figsize = [12,6])
        axSECV  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        axSECV[0][0].set_title("SECV - GML data per down movements", fontsize=14, fontweight="bold")
        axSECV[0][1].set_title("SECV - GML data per up movements", fontsize=14, fontweight="bold")  
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)
        
        
        
        
        
            
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                        
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
        
            if(down):
   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]         
            else:
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    

        le=leminb
                   
        axSECV[0][0].plot(arrLB[:le], arrGB[:le],color=colors[k])
     
        le=lemin
                   
        axSECV[0][1].plot(arrL[:le], arrG[:le],color=colors[k])
            
        axSECV[0][0].set_xlabel("LF", fontsize=13)
        axSECV[0][0].set_ylabel("GF", fontsize=13)
           
        axSECV[0][1].set_xlabel("LF", fontsize=13)
        axSECV[0][1].set_ylabel("GF", fontsize=13)
        arrmfB[:le]=0
        arrmf[:le]=0
        
        k = k+1
     
            
        plt.show()
    



#Graphe pour SECVR   
#Graphe pour SECVR   
subjects = ["julSECVR","SimonSECVR","SophieSECVR","JulianSECVR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['b', 'g', 'r', 'c', 'b', 'g', 'r', 'c']
lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
#fig.suptitle("%s_00%d" % (s,trial))
k = 0  

for s in subjects :
    
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        fig = plt.figure(figsize = [12,6])
        axSECVR  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)
        xnum=8000
        x=np.arange(xnum)
        axSECVR[0][0].set_title("SECVR - GML data per down movements", fontsize=14, fontweight="bold")
        axSECVR[0][1].set_title("SECVR - GML data per up movements", fontsize=14, fontweight="bold")
        
        
        
#        fig = plt.figure(figsize = [15,7])
#        axSECV  = fig.subplots(1,2,squeeze=False)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
#            
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
            frac = np.full(xnum,0).astype(float)
            fracB = np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                       
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                         
            if(down):
                #axSECVR[0][0].plot(yLB[:le], yGB[:le],color='g',linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]
                    
#                    
            else:
                    
                #axSECVR[0][1].plot(yL[:le], yG[:le],color='g',linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    
        
    
        le=leminb               
        axSECVR[0][0].plot(arrLB[:le], arrGB[:le],color=colors[k])
     
        le=lemin
        axSECVR[0][1].plot(arrL[:le], arrG[:le],color=colors[k])
            
        axSECVR[0][0].set_xlabel("LF", fontsize=13)
        axSECVR[0][0].set_ylabel("GF", fontsize=13)
           
        axSECVR[0][1].set_xlabel("LF ", fontsize=13)
        axSECVR[0][1].set_ylabel("GF", fontsize=13)
        k = k+1
#        
        plt.show()
        
#Graphe de SEF        
subjects = ["julSEF","SimonSEF","SophieSEF","JulianSEF"]        
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']

lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'r'] 
d = 0            
c = 0
k = 0
for s in subjects :
    
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        fig = plt.figure(figsize = [12,6])
        axSEF  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        axSEF[0][0].set_title("SEF - GML data per down movements", fontsize=14, fontweight="bold")
        axSEF[0][1].set_title("SEF - GML data per up movements", fontsize=14, fontweight="bold") 
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)
                
        
        
        
#        fig = plt.figure(figsize = [12,6])
#        axSECV  = fig.subplots(1,2,squeeze=False)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
#            
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
            frac = np.full(xnum,0).astype(float)
            fracB = np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    fracB[i]=yGB[i]/yLB[i]    
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    frac[i]=yG[i]/yL[i]         
            if(down):
                #axSECV[0][0].plot(x[:le], fracB[:le],color=colors[c],linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]
                    arrfB[i]= arrGB[i]/arrLB[i]
                    
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    arrf[i]= arrG[i]/arrL[i]
                    
            
            
        
        le=leminb
        axSEF[0][0].plot(arrLB[:le], arrGB[:le],color=colors[k])
     
        le=lemin
        axSEF[0][1].plot(arrL[:le], arrG[:le],color=colors[k])
            
        axSEF[0][0].set_xlabel("LF", fontsize=13)
        axSEF[0][0].set_ylabel("GF", fontsize=13)
           
        axSEF[0][1].set_xlabel("LF", fontsize=13)
        axSEF[0][1].set_ylabel("GF", fontsize=13)
        
        k=k+1
        
        plt.show()
        
    
#Graphe pour SEFR -
subjects = ["julSEFR","SimonSEFR","SophieSEFR","JulianSEFR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]

lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
colors = ['b', 'g', 'r', 'c', 'b', 'g', 'r', 'c']  
c = 0
d=0
k=0
for s in subjects :
    
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        
        
        
        
#        fig = plt.figure(figsize = [12,6])
#        axSECV  = fig.subplots(1,2,squeeze=False)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
            
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
        fig = plt.figure(figsize = [12,6])
        axSEFR  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        axSEFR[0][0].set_title("SEFR - GML data per down movements", fontsize=14, fontweight="bold")
        axSEFR[0][1].set_title("SEFR - GML data per up movements", fontsize=14, fontweight="bold")
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
            frac = np.full(xnum,0).astype(float)
            fracB = np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    fracB[i]=yGB[i]/yLB[i]    
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    frac[i]=yG[i]/yL[i]         
            if(down):
                #axSECV[0][0].plot(x[:le], fracB[:le],color=colors[c],linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]
                    arrfB[i]= arrGB[i]/arrLB[i]
                    
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    arrf[i]= arrG[i]/arrL[i]
                    
           
             
        le=leminb
        axSEFR[0][0].plot(arrLB[:le], arrGB[:le],color=colors[k])
     
        le=lemin
        axSEFR[0][1].plot(arrL[:le], arrG[:le],color=colors[k])
            
        axSEFR[0][0].set_xlabel("LF", fontsize=13)
        axSEFR[0][0].set_ylabel("GF", fontsize=13)
           
        axSEFR[0][1].set_xlabel("LF", fontsize=13)
        axSEFR[0][1].set_ylabel("GF", fontsize=13)
        
        
        k = k+1
            
        plt.show()

#Grpahe pour EH - 
subjects = ["JULEH","SimonEH","SophieEH","JulianEH"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'r']
d=0            
c = 0
k = 0
for s in subjects :
    
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        
        
        
        
#        fig = plt.figure(figsize = [15,7])
#        axSECV  = fig.subplots(1,2,squeeze=False)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
        fig = plt.figure(figsize = [12,6])
        axEH  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        axEH[0][0].set_title("EH - GML data per down movements", fontsize=14, fontweight="bold")
        axEH[0][1].set_title("EH - GML data per up movements", fontsize=14, fontweight="bold")
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)    
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
            frac = np.full(xnum,0).astype(float)
            fracB = np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    fracB[i]=yGB[i]/yLB[i]    
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    frac[i]=yG[i]/yL[i]         
            if(down):
                #axSECV[0][0].plot(x[:le], fracB[:le],color=colors[c],linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]
                    arrfB[i]= arrGB[i]/arrLB[i]
                    
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    arrf[i]= arrG[i]/arrL[i]
                    
                   
        le=leminb
        axEH[0][0].plot(arrGB[:le], arrLB[:le],color=colors[k])
     
        le=lemin
        axEH[0][1].plot(arrG[:le], arrL[:le],color=colors[k])
            
        axEH[0][0].set_xlabel("LF", fontsize=13)
        axEH[0][0].set_ylabel("GF", fontsize=13)
           
        axEH[0][1].set_xlabel("LF", fontsize=13)
        axEH[0][1].set_ylabel("GF", fontsize=13)
        arrmfB[:le]=0
        arrmf[:le]=0
        
        k = k+1
            
        plt.show()




#Graphe  pour EHR
subjects = ["JULEHR","SimonEHR","SophieEHR","JulianEHR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'r']
d=0           
c = 0
k=0
for s in subjects :
    
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        
        
        
        
#        fig = plt.figure(figsize = [15,7])
#        axSECV  = fig.subplots(1,2,squeeze=False)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
        fig = plt.figure(figsize = [12,6])
        axEHR  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        axEHR[0][0].set_title("EHR - GML data per down movements", fontsize=14, fontweight="bold")
        axEHR[0][1].set_title("EHR - GML data per up movements", fontsize=14, fontweight="bold")  
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)    
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
            frac = np.full(xnum,0).astype(float)
            fracB = np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    fracB[i]=yGB[i]/yLB[i]    
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    frac[i]=yG[i]/yL[i]         
            if(down):
                #axSECV[0][0].plot(x[:le], fracB[:le],color=colors[c],linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]
                    arrfB[i]= arrGB[i]/arrLB[i]
                    
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    arrf[i]= arrG[i]/arrL[i]
                       
        le=leminb
        axEHR[0][0].plot(arrLB[:le], arrGB[:le],color=colors[k])
     
        le=lemin
        axEHR[0][1].plot(arrL[:le], arrG[:le],color=colors[k])
            
        axEHR[0][0].set_xlabel("LF", fontsize=13)
        axEHR[0][0].set_ylabel("GF", fontsize=13)
           
        axEHR[0][1].set_xlabel("LF", fontsize=13)
        axEHR[0][1].set_ylabel("GF", fontsize=13)
        arrmfB[:le]=0
        arrmf[:le]=0
        k = k+1
        
            
        plt.show()
        
        
#Graphe pour EB - 
subjects = ["JULEB","SimonEB","SophieEB","JulianEB"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'r']
c = 0
d=0
k =0
for s in subjects :
    
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        fig = plt.figure(figsize = [12,6])
        axEB  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        axEB[0][0].set_title("EB - GML data per down movements", fontsize=14, fontweight="bold")
        axEB[0][1].set_title("EB - GML data per up movements", fontsize=14, fontweight="bold")  
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)
        
        
        
#        fig = plt.figure(figsize = [15,7])
#        axSECV  = fig.subplots(1,2,squeeze=False)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
#            
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
            frac = np.full(xnum,0).astype(float)
            fracB = np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    fracB[i]=yGB[i]/yLB[i]    
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    frac[i]=yG[i]/yL[i]         
            if(down):
                #axSECV[0][0].plot(x[:le], fracB[:le],color=colors[c],linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]
                    arrfB[i]= arrGB[i]/arrLB[i]
                    
                    
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    arrf[i]= arrG[i]/arrL[i]
                    
                          
        le=leminb
        axEB[0][0].plot(arrLB[:le], arrGB[:le],color=colors[k])
     
        le=lemin
        axEB[0][1].plot(arrL[:le], arrG[:le],color=colors[k])
            
        axEB[0][0].set_xlabel("LF", fontsize=13)
        axEB[0][0].set_ylabel("GF", fontsize=13)
           
        axEB[0][1].set_xlabel("LF", fontsize=13)
        axEB[0][1].set_ylabel("GF", fontsize=13)
        arrmfB[:le]=0
        arrmf[:le]=0
        k = k+1
        
            
        plt.show()

#Graphe pour EBR - 
subjects = ["JULEBR","SimonEBR","SophieEBR","JulianEBR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
lab = ['Julien','Julien','Simon','Simon','Sophie','Sophie','Julian','Julian']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'r']
d=0            
c = 0
for s in subjects :
    k=0
    for trial in range(1,ntrials+1): 
        # Set data path
        glm_path = "DataGroupe4/%s_00%d.glm" % (s,trial)

        # Import data 
        glm_df = glm.import_data(glm_path)
        
        baseline = range(0,400)        
        # Normal Force exerted by the thumb
        NF_thumb = glm_df.loc[:,'Fygl']-np.nanmean(glm_df.loc[baseline,'Fygl'])
        # Vertical Tangential Force exerted by the thumb
        TFx_thumb  = glm_df.loc[:,'Fxgl']-np.nanmean(glm_df.loc[baseline,'Fxgl'])
        #Horizontal Tangential Force exerted by the thumb
        TFz_thumb  = glm_df.loc[:,'Fzgl']-np.nanmean(glm_df.loc[baseline,'Fzgl'])


        # Normal Force exerted by the index
        NF_index = -(glm_df.loc[:,'Fygr']-np.nanmean(glm_df.loc[baseline,'Fygr']))
        # Vertical Tangential Force exerted by the index
        TFx_index = glm_df.loc[:,'Fxgr']-np.nanmean(glm_df.loc[baseline,'Fxgr'])
        #Horizontal Tangential Force exerted by the index
        TFz_index = glm_df.loc[:,'Fzgr']-np.nanmean(glm_df.loc[baseline,'Fzgr'])
        
        
        
        #%% Get acceleration, LF and GF
        time  = glm_df.loc[:,'time'].to_numpy()
        accX  = glm_df.loc[:,'LowAcc_X'].to_numpy()*(-9.81)
        accX  = accX-np.nanmean(accX[baseline])
        GF    = glm_df.loc[:,'GF'].to_numpy()
        GF    = GF-np.nanmean(GF[baseline])
        LFv   = TFx_thumb+TFx_index
        LFh   = TFz_thumb+TFz_index
        LF    = np.hypot(LFv,LFh)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        accX = glm.filter_signal(accX, fs = freqAcq, fc = freqFiltAcc)
        GF   = glm.filter_signal(GF,   fs = freqAcq, fc = freqFiltForces)
        LF   = glm.filter_signal(LF,   fs = freqAcq, fc = freqFiltForces)
        LFv   = glm.filter_signal(LFv,   fs = freqAcq, fc = freqFiltForces)
        LFh   = glm.filter_signal(LFh,   fs = freqAcq, fc = freqFiltForces)
        fig = plt.figure(figsize = [12,6])
        axEBR  = fig.subplots(1,2,squeeze=False)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        axEBR[0][0].set_title("EBR - GML data per down movements", fontsize=14, fontweight="bold")
        axEBR[0][1].set_title("EBR - GML data per up movements", fontsize=14, fontweight="bold")  
        arrm = np.full(xnum,0).astype(float)
        arrm2 = np.full(xnum,0).astype(float)
        arrmB = np.full(xnum,0).astype(float)
        arrmB2 = np.full(xnum,0).astype(float)
        arrmf = np.full(xnum,0).astype(float)
        arrmfB = np.full(xnum,0).astype(float)
        
        
#        
#        fig = plt.figure(figsize = [15,7])
#        axSECV  = fig.subplots(1,2,squeeze=False)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
#            
        down=True
        arrL=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        arrf = np.full(xnum,0).astype(float)
        arrfB = np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
    
        for e in range (len(segmentations)-1):   
            yL=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)
            frac = np.full(xnum,0).astype(float)
            fracB = np.full(xnum,0).astype(float)
                
            le=int(segmentations[e+1]-segmentations[e]) 
                
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                        
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    fracB[i]=yGB[i]/yLB[i]    
                else:
                        
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    frac[i]=yG[i]/yL[i]         
            if(down):
                #axSECV[0][0].plot(x[:le], fracB[:le],color=colors[c],linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]=yLB[i]
                    arrGB[i]=yGB[i]
                    arrfB[i]= arrGB[i]/arrLB[i]
                    
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]=yL[i]
                    arrG[i]=yG[i]
                    arrf[i]= arrG[i]/arrL[i]
                    
        le=leminb
        axEBR[0][0].plot(arrLB[:le], arrGB[:le],color=colors[k])
     
        le=lemin
        axEBR[0][1].plot(arrL[:le], arrG[:le],color=colors[k])
            
        axEBR[0][0].set_xlabel("LF", fontsize=13)
        axEBR[0][0].set_ylabel("GF", fontsize=13)
           
        axEBR[0][1].set_xlabel("LF", fontsize=13)
        axEBR[0][1].set_ylabel("GF", fontsize=13)
        arrmfB[:le]=0
        arrmf[:le]=0
        k = k+1
        c=c+1
        plt.show()