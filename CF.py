#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 11:43:21 2021

@author: sophiealicandro
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal 
import os
import glm_data_processing as glm
import derive as der
from script_coda  import get_segmentations

###GRAPHE DES CF
nf=np.arange(1,50,0.5)

dataJulianT = 1.843672*pow(nf,0.662438-1)
dataJulianI = 1.809796 * pow(nf, 0.695003-1)

dataJulienT = 1.589808 * pow(nf, 0.822373-1)
dataJulienI = 1.653512 * pow(nf,0.848636-1)

dataSimonT = 1.412046 * pow(nf, 0.634535-1)
dataSimonI = 1.163868 * pow(nf, 0.743537-1)

dataSophieT = 1.945243 * pow(nf,0.647709-1)
dataSophieI = 1.474812 * pow(nf, 0.771625-1)



#plt.plot(nf,dataJulianI, label= 'Julian')
#plt.plot(nf, dataJulienI, label='Julien')
#plt.plot(nf, dataSimonI, label = 'Simon')
#plt.plot(nf, dataSophieI, label ='Sophie')

#ax[0][0].set_title("Coéfficient de friction des différents sujets(pouce)", fontsize=14, fontweight="bold")


#plt.title('Coéfficient de friction des différents sujets(index)')
#plt.xlabel('Force normale [N]')
#plt.ylabel('µ')
#plt.legend()



subjects = ["SophieSECV"] #Names of subjects
ntrials = 1 #Number of trials for each subject

# Double for-loop that runs thrgough all subjects and trials
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
        
        GFSP    = glm_df.loc[:,'GF'].to_numpy()
        LFvSP   = TFx_thumb+TFx_index
        LFhSP   = TFz_thumb+TFz_index
        LFSP    = np.hypot(LFvSP,LFhSP)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        GFSP   = glm.filter_signal(GFSP,   fs = freqAcq, fc = freqFiltForces)
        LFSP   = glm.filter_signal(LFSP,   fs = freqAcq, fc = freqFiltForces)
        LFvSP   = glm.filter_signal(LFvSP,   fs = freqAcq, fc = freqFiltForces)
        LFhSP   = glm.filter_signal(LFhSP,   fs = freqAcq, fc = freqFiltForces)
        
      
        
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
        
        le=int(segmentations[2]-segmentations[1])        
        arrLSP=np.full(le,0).astype(float)      
        GFSophieU=np.full(le,0).astype(float)
        
        for i in range (le):    
            arrLSP[i]=LFSP[segmentations[1]+i]
            GFSophieU = GFSP[segmentations[1]+i]
           
        led=int(segmentations[3]-segmentations[2])
        arrLDSP=np.full(led,0).astype(float)
        GFSophieD=np.full(led,0).astype(float)
        
        for i in range(led):
            arrLDSP[i]=LFSP[segmentations[2]+i]
            GFSophieD = GFSP[segmentations[2]+i]
            
            
            
subjects = ["JulianSECV"] #Names of subjects
ntrials = 1 #Number of trials for each subject

# Double for-loop that runs thrgough all subjects and trials
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
       
        GFJI    = glm_df.loc[:,'GF'].to_numpy()
        LFvJI   = TFx_thumb+TFx_index
        LFhJI   = TFz_thumb+TFz_index
        LFJI    = np.hypot(LFvJI,LFhJI)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces
        
        GFJI   = glm.filter_signal(GFJI,   fs = freqAcq, fc = freqFiltForces)
        LFJI   = glm.filter_signal(LFJI,   fs = freqAcq, fc = freqFiltForces)
        LFvJI = glm.filter_signal(LFvJI,   fs = freqAcq, fc = freqFiltForces)
        LFhJI   = glm.filter_signal(LFhJI,   fs = freqAcq, fc = freqFiltForces)
        
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
        #down=True
        
        le=int(segmentations[2]-segmentations[1])
        
        arrLJI=np.full(le,0).astype(float) 
        GFJulianU=np.full(le,0).astype(float)

        #arrLDJI=np.full(le,0)
        
       # yLJI=np.full(le,0).astype(float)                
       # yLDJI=np.full(le,0).astype(float)
                
        for i in range (le):
            #if(down):
        #        yLDJI[i]=LF[(segmentations[1]+i)]
                #arrLDJI[i]+= yLDJI[i]/10.0
                
         #       down=False
            #else:
                arrLJI[i]=LFJI[segmentations[1]+i]
                GFJulianU = GFJI[segmentations[1]+i]

                #arrLJI[i]+=yLJI[i]/10.0
              #  down=True
        led=int(segmentations[3]-segmentations[2])
        arrLDJI=np.full(led,0).astype(float)
        GFJulianD=np.full(led,0).astype(float)

        for i in range(led):
            arrLDJI[i]=LFJI[segmentations[2]+i]
            GFJulianD = GFJI[segmentations[2]+i]

  
subjects = ["SimonSECV"] #Names of subjects
ntrials = 1 #Number of trials for each subject

# Double for-loop that runs thrgough all subjects and trials
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
       
        GFSM    = glm_df.loc[:,'GF'].to_numpy()
        LFvSM   = TFx_thumb+TFx_index
        LFhSM   = TFz_thumb+TFz_index
        LFSM    = np.hypot(LFvSM,LFhSM)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces
        
        GFSM   = glm.filter_signal(GFSM,   fs = freqAcq, fc = freqFiltForces)
        LFSM   = glm.filter_signal(LFSM,   fs = freqAcq, fc = freqFiltForces)
        LFvSM   = glm.filter_signal(LFvSM,   fs = freqAcq, fc = freqFiltForces)
        LFhSM   = glm.filter_signal(LFhSM,   fs = freqAcq, fc = freqFiltForces)
               
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
        
        le=int(segmentations[2]-segmentations[1])
        
        arrLSM=np.full(le,0).astype(float) 
        
        GFSimonU=np.full(le,0).astype(float)

                
        for i in range (le):
       
                arrLSM[i]=LFSM[segmentations[1]+i]
                GFSimonU = GFSM[segmentations[1]+i]

               
        led=int(segmentations[3]-segmentations[2])
        arrLDSM=np.full(led,0).astype(float)
        GFSimonD=np.full(led,0).astype(float)

        for i in range(led):
            arrLDSM[i]=LFSM[segmentations[2]+i]
            GFSimonD = GFSM[segmentations[2]+i]

       
            
            
            
subjects = ["JulSECV"] #Names of subjects
ntrials = 1 #Number of trials for each subject

# Double for-loop that runs thrgough all subjects and trials
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
       
        GFJE    = glm_df.loc[:,'GF'].to_numpy()
        LFvJE   = TFx_thumb+TFx_index
        LFhJE   = TFz_thumb+TFz_index
        LFJE    = np.hypot(LFvJE,LFhJE)
        
        # %%Filter data
        freqAcq=800 #Frequence d'acquisition des donnees
        freqFiltAcc=20 #Frequence de coupure de l'acceleration
        freqFiltForces=20 #Frequence de coupure des forces

        GFJE   = glm.filter_signal(GFJE,   fs = freqAcq, fc = freqFiltForces)
        LFJE   = glm.filter_signal(LFJE,   fs = freqAcq, fc = freqFiltForces)
        LFvJE   = glm.filter_signal(LFvJE,   fs = freqAcq, fc = freqFiltForces)
        LFhJE   = glm.filter_signal(LFhJE,   fs = freqAcq, fc = freqFiltForces)
        
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
        le=int(segmentations[3]-segmentations[2])
        
        arrLJE=np.full(le,0).astype(float) 
        GFJulienU=np.full(le,0).astype(float)

       
        for i in range (le):
         
                arrLJE[i]=LFJE[segmentations[2]+i]
                GFJulienU = GFJE[segmentations[2]+i]

                
        led=int(segmentations[2]-segmentations[1])
        arrLDJE=np.full(led,0).astype(float)
        GFJulienD=np.full(led,0).astype(float)

        
        for i in range(led):
            arrLDJE[i]=LFJE[segmentations[1]+i]
            GFJulienD = GFJE[segmentations[1]+i]

          
 
monaxe=np.arange(1200)               
###GRAPHE DES SF UP 
LFJulianTU = pow(arrLJI/2*1.843672,1/0.662438)
LFJulienTU = pow(arrLJE/2*1.589808, 1/0.822373)
LFSimonTU = pow(arrLSM/2*1.412046,1/0.634535)
LFSophieTU = pow(arrLSP/(2*1.945243),1/0.647709)

LFJulianIU = pow(arrLJI/2*1.809796, 1/0.695003)
LFJulienIU = pow(arrLJE/2*1.653512, 1/0.848636)
LFSimonIU = pow(arrLSP/2*1.163868,1/0.743537)
LFSophieIU = pow(arrLSP/(2*1.474812),1/0.771625)

###GRAPHE DES SF DOWN
LFJulianTD = pow(arrLDJI/2*1.843672,1/0.662438)
LFJulienTD = pow(arrLDJE/2*1.589808, 1/0.822373)
LFSimonTD = pow(arrLDSM/2*1.412046,1/0.634535)
LFSophieTD = pow(arrLDSP/(2*1.945243),1/0.647709)

LFJulianID = pow(arrLDJI/2*1.809796, 1/0.695003)
LFJulienID = pow(arrLDJE/2*1.653512, 1/0.848636)
LFSimonID = pow(arrLDSM/2*1.163868,1/0.743537)
LFSophieID = pow(arrLDSP/(2*1.474812),1/0.771625)

#Graphes pour UP

fig = plt.figure(figsize = [15,7])
ax = fig.subplots(3)

"""
ax[0].set_title("Slip force des différents sujets - Pouce - UP", fontsize=14, fontweight="bold")
ax[1].set_title("Slip force des différents sujets - Index - UP", fontsize=14, fontweight="bold")
ax[2].set_title("Force load des différents sujets - UP", fontsize=14, fontweight="bold")
        
        
ax[0].plot(np.arange(len(LFJulianTU)), LFJulianTU, label='Julian thumb', color = 'b')
ax[2].plot(np.arange(len(arrLJI)), arrLJI, color ='b',linestyle='dotted', label = 'LF Julian')
ax[0].plot(np.arange(len(LFJulienTU)), LFJulienTU, label ='Julien thumb', color='r')
ax[2].plot(np.arange(len(arrLJE)), arrLJE, color='r',linestyle='dotted', label = 'LF Julien')
ax[0].plot(np.arange(len(LFSimonTU)), LFSimonTU, label='Simon thumb', color='g')
ax[2].plot(np.arange(len(arrLSM)), arrLSM, color='g',linestyle='dotted', label = 'LF Simon')
ax[0].plot(np.arange(len(LFSophieTU)), LFSophieTU, label = 'Sophie thumb', color='m')
ax[2].plot(np.arange(len(arrLSP)), arrLSP, color='m',linestyle='dotted', label = 'LF Sophie')
ax[1].plot(np.arange(len(LFSophieIU)), LFSophieIU, label = 'Sophie Index',  color='m', linestyle='dashed')
ax[1].plot(np.arange(len(LFSimonIU)), LFSimonIU, label ='Simon Index', color='g', linestyle='dashed')
ax[1].plot(np.arange(len(LFJulienIU)), LFJulienIU, label='Julien Index', color='r', linestyle='dashed')
ax[1].plot(np.arange(len(LFJulianIU)), LFJulianIU, label ='Julian Index',color ='b', linestyle='dashed')

plt.legend()
ax[0].set_xlabel("Time [iteration]", fontsize=13)
ax[0].set_ylabel("SF[N]", fontsize=13)
ax[1].set_xlabel("Time [iteration]", fontsize=13)
ax[1].set_ylabel("SF[N]", fontsize=13)
ax[2].set_xlabel("Time [iteration]", fontsize=13)
ax[2].set_ylabel("LF[N]", fontsize=13)
"""

#Graphes pour DOWN

ax[0].set_title("Slip force des différents sujets - Pouce - DOWN", fontsize=14, fontweight="bold")
ax[1].set_title("Slip force des différents sujets - Index - DOWN", fontsize=14, fontweight="bold")
ax[2].set_title("Force load des différents sujets - DOWN", fontsize=14, fontweight="bold")
        
        
ax[0].plot(np.arange(len(LFJulianTD)), LFJulianTD, label='Julian thumb', color = 'b')
ax[2].plot(np.arange(len(arrLDJI)), arrLDJI, color ='b',linestyle='dotted', label = 'LF Julian')
ax[0].plot(np.arange(len(LFJulienTD)), LFJulienTD, label ='Julien thumb', color='r')
ax[2].plot(np.arange(len(arrLDJE)), arrLDJE, color='r',linestyle='dotted', label = 'LF Julien')
ax[0].plot(np.arange(len(LFSimonTD)), LFSimonTD, label='Simon thumb', color='g')
ax[2].plot(np.arange(len(arrLDSM)), arrLDSM, color='g',linestyle='dotted', label = 'LF Simon')
ax[0].plot(np.arange(len(LFSophieTD)), LFSophieTD, label = 'Sophie thumb', color='m')
ax[2].plot(np.arange(len(arrLDSP)), arrLDSP, color='m',linestyle='dotted', label = 'LF Sophie')
ax[1].plot(np.arange(len(LFSophieID)), LFSophieID, label = 'Sophie Index',  color='m', linestyle='dashed')
ax[1].plot(np.arange(len(LFSimonID)), LFSimonID, label ='Simon Index', color='g', linestyle='dashed')
ax[1].plot(np.arange(len(LFJulienID)), LFJulienID, label='Julien Index', color='r', linestyle='dashed')
ax[1].plot(np.arange(len(LFJulianID)), LFJulianID, label ='Julian Index',color ='b', linestyle='dashed')

plt.legend()
ax[0].set_xlabel("Time [iteration]", fontsize=13)
ax[0].set_ylabel("SF[N]", fontsize=13)
ax[1].set_xlabel("Time [iteration]", fontsize=13)
ax[1].set_ylabel("SF[N]", fontsize=13)
ax[2].set_xlabel("Time [iteration]", fontsize=13)
ax[2].set_ylabel("LF[N]", fontsize=13)


#SAEFTY MARGIN UP

SMJITU = LFJulianTU - GFJulianU
SMJIIU = LFJulianIU - GFJulianU

SMJETU = LFJulienTU - GFJulienU
SMJEIU = LFJulienIU - GFJulienU

SMSMTU = LFSimonTU - GFSimonU
SMSMIU = LFSimonIU - GFSimonU

SMSPTU = LFSophieTU - GFSophieU
SMSPIU = LFSophieIU - GFSophieU 

SMJITD = LFJulianTD - GFJulianD
SMJIID = LFJulianID - GFJulianD

SMJETD = LFJulienTD - GFJulienD
SMJEID = LFJulienID - GFJulienD

SMSMTD = LFSimonTD - GFSimonD
SMSMID = LFSimonID - GFSimonD

SMSPTD = LFSophieTD - GFSophieD
SMSPID = LFSophieID - GFSophieD 

#graphes UP
"""
fig = plt.figure(figsize = [15,7])
ax = fig.subplots(2)

ax[0].set_title('Safety margin des différents sujets-Pouce-UP')
ax[0].plot(np.arange(len(SMJITU)),SMJITU, color='b', label ='Julian')
ax[0].plot(np.arange(len(SMJETU)),SMJETU, color='r', label = 'Julien')
ax[0].plot(np.arange(len(SMSMTU)),SMSMTU, color='g', label = 'Simon')
ax[0].plot(np.arange(len(SMSPTU)), SMSPTU, color = 'm', label = 'Sophie')
ax[0].legend()
ax[0].set_xlabel('Time [iteration]')
ax[0].set_ylabel('Safety margin[N]')

ax[1].set_title('Safety margin des différents sujets-Index-UP')
ax[1].plot(np.arange(len(SMJIIU)),SMJIIU, color='b', label ='Julian')
ax[1].plot(np.arange(len(SMJEIU)),SMJEIU, color='r', label = 'Julien')
ax[1].plot(np.arange(len(SMSMIU)),SMSMIU, color='g', label = 'Simon')
ax[1].plot(np.arange(len(SMSPIU)), SMSPIU, color = 'm', label = 'Sophie')
ax[1].legend()
ax[1].set_xlabel('Time [iteration]')
ax[1].set_ylabel('Safety margin[N]')
"""
#graphes DOWn

"""
fig = plt.figure(figsize = [15,7])
ax = fig.subplots(2)

ax[0].set_title('Safety margin des différents sujets-Pouce-DOWN')
ax[0].plot(np.arange(len(SMJITD)),SMJITD, color='b', label ='Julian')
ax[0].plot(np.arange(len(SMJETD)),SMJETD, color='r', label = 'Julien')
ax[0].plot(np.arange(len(SMSMTD)),SMSMTD, color='g', label = 'Simon')
ax[0].plot(np.arange(len(SMSPTD)), SMSPTD, color = 'm', label = 'Sophie')
ax[0].legend()
ax[0].set_xlabel('Time [iteration]')
ax[0].set_ylabel('Safety margin[N]')

ax[1].set_title('Safety margin des différents sujets-Index-DOWN')
ax[1].plot(np.arange(len(SMJIID)),SMJIID, color='b', label ='Julian')
ax[1].plot(np.arange(len(SMJEID)),SMJEID, color='r', label = 'Julien')
ax[1].plot(np.arange(len(SMSMID)),SMSMID, color='g', label = 'Simon')
ax[1].plot(np.arange(len(SMSPID)), SMSPID, color = 'm', label = 'Sophie')
ax[1].legend()
ax[1].set_xlabel('Time [iteration]')
ax[1].set_ylabel('Safety margin[N]')
"""

 