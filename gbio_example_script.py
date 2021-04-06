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

# Fermeture des figures ouvertes
plt.close('all')

#subjects = ["JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF"] #Names of subjects
#subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF"] #Names of subjects
#Notes  not existant: EB2, EHR
subjects =["Julien","Julian","Simon","Sophie"]

kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
expe=["EB","EBR","EH","EHR","SECV","SECVR","SEF","SEFR"]#,"EH","SECV","SEF","EBR","EHR","SECVR","SEFR"]

ntrials = 2 #Number of trials for each subject

masterG={"EB":np.zeros((80,600)),"EBcount":0,"EBD":np.zeros((80,600)),"EBDcount":0,
        "EBR":np.zeros((80,600)),"EBRcount":0,"EBRD":np.zeros((80,600)),"EBRDcount":0,
        "EH":np.zeros((80,600)),"EHcount":0,"EHD":np.zeros((80,600)),"EHDcount":0,
        "EHR":np.zeros((80,600)),"EHRcount":0,"EHRD":np.zeros((80,600)),"EHRDcount":0,
        "SECV":np.zeros((80,600)),"SECVcount":0,"SECVD":np.zeros((80,600)),"SECVDcount":0,
        "SECVR":np.zeros((80,600)),"SECVRcount":0,"SECVRD":np.zeros((80,600)),"SECVRDcount":0,
        "SEF":np.zeros((80,600)),"SEFcount":0,"SEFD":np.zeros((80,600)),"SEFDcount":0,
        "SEFR":np.zeros((80,600)),"SEFRcount":0,"SEFRD":np.zeros((80,600)),"SEFRDcount":0}

masterL={"EB":np.zeros((80,600)),"EBcount":0,"EBD":np.zeros((80,600)),"EBDcount":0,
        "EBR":np.zeros((80,600)),"EBRcount":0,"EBRD":np.zeros((80,600)),"EBRDcount":0,
        "EH":np.zeros((80,600)),"EHcount":0,"EHD":np.zeros((80,600)),"EHDcount":0,
        "EHR":np.zeros((80,600)),"EHRcount":0,"EHRD":np.zeros((80,600)),"EHRDcount":0,
        "SECV":np.zeros((80,600)),"SECVcount":0,"SECVD":np.zeros((80,600)),"SECVDcount":0,
        "SECVR":np.zeros((80,600)),"SECVRcount":0,"SECVRD":np.zeros((80,600)),"SECVRDcount":0,
        "SEF":np.zeros((80,600)),"SEFcount":0,"SEFD":np.zeros((80,600)),"SEFDcount":0,
        "SEFR":np.zeros((80,600)),"SEFRcount":0,"SEFRD":np.zeros((80,600)),"SEFRDcount":0}


ntrials = 2 #Number of trials for each subject




# Double for-loop that runs thrgough all subjects and trials
subject_number=0;
for s in subjects:
    for exp in expe:
        for trial in range(1,ntrials+1): 
            # Set data path
            glm_path = "DataGroupe4/%s%s_00%d.glm" % (s,exp,trial)
    
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

                        


            
            xnum=10000
            x=np.arange(0,xnum/800,1/800)
            
    
            
            down=True
            
            arrx=np.full(xnum,0).astype(float)
            arry=np.full(xnum,0).astype(float)
            arrz=np.full(xnum,0).astype(float)
    
    
            
            lemin =1000
            leminb=10000
            
            count=0

            

            for e in range (len(segmentations)-1):
                yL=np.full(xnum,0).astype(float)
                yLB=np.full(xnum,0).astype(float)
                yG=np.full(xnum,0).astype(float)
                yGB=np.full(xnum,0).astype(float)
    
    
    
                le=int(segmentations[e+1]-segmentations[e]) 
                for i in range (int(segmentations[e+1]-segmentations[e])):
                    if(down):
                        yLB[i]=LF[0][(segmentations[e]+i)]
                        yGB[i]=GF[0][(segmentations[e]+i)]

                    else:
                        yL[i]=LF[0][(segmentations[e]+i)]
                        yG[i]=GF[0][(segmentations[e]+i)]

    
                if(down):
                    masterL[s+"D"][masterL[s+"Dcount"]]=yLB[:600]
                    masterL[s+"Dcount"]+=1
                    masterG[s+"D"][masterG[s+"Dcount"]]=yGB[:600]
                    masterG[s+"Dcount"]+=1
                    
                    down=False

                    if(leminb>le):leminb=le
    
    
    
                        
                else:
                    masterL[s][masterL[s+"count"]]=yL[:600]
                    masterL[s+"count"]+=1
                    masterG[s][masterG[s+"count"]]=yG[:600]
                    masterG[s+"count"]+=1
                    
                    down=True
                    if(lemin>le):lemin=le

            
