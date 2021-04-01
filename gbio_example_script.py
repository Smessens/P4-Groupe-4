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

#subjects = ["JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF"] #Names of subjects
#subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF"] #Names of subjects
#Notes  not existant: EB2, EHR
subjects = ["SophieEB","SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF"] #Names of subjects
#subjects = ["JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
#subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","SophieEB","SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF","JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects

kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
ntrials = 2 #Number of trials for each subject




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
        """

        #%% Basic plot of the data
        fig = plt.figure(figsize = [15,7])
        ax  = fig.subplots(3,1)
        fig.suptitle("%s_00%d" % (s,trial))


        for i in range(0,len(segmentations)-1):
                            ax[0].axvline(x=time[int(segmentations[i])])



        
        
        ax[0].plot(time, accX)
      #  ax[0].plot(time[ipk],accX[ipk], linestyle='', marker='o', 
      #             markerfacecolor='None', markeredgecolor='r')
        ax[0].set_ylabel("Acceleration [m/s^2]", fontsize=13)
        ax[0].set_title("GLM data", fontsize=14, fontweight="bold")
        ax[0].set_xlim([5,40])
        
        # Putting grey patches for cycles

        
        ax[1].plot(time,LF, label="LF")
        ax[1].plot(time,GF, label="GF")
        ax[1].legend(fontsize=12)
        ax[1].set_xlabel("Time [s]", fontsize=13)
        ax[1].set_ylabel("Forces [N]", fontsize=13)
        ax[1].set_xlim([5,40])
        
        ax[2].plot(time,dGF)
        ax[2].set_xlabel("Time [s]", fontsize=13)
        ax[2].set_ylabel("GF derivative [N/s]", fontsize=13)
        ax[2].set_xlim([5,40])
    
        #%% Save the figure as png file. Creates a folder "figures" first if it
        # doesn't exist
        if not os.path.exists('figures'):
            os.makedirs('figures')
        
        for e in kind:
            if e  in  name:
                fig.savefig("figures/" + e +  "/%s_%d_GML_DATA.png" %(s,trial))

        """

                    
        fig = plt.figure(figsize = [15,7])
        ax  = fig.subplots(4,2)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=8000
        x=np.arange(xnum)
        ax[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
        ax[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")
        

        
        down=True
        
        arra=np.full(xnum,0).astype(float)
        arrL=np.full(xnum,0).astype(float)
        arrd=np.full(xnum,0).astype(float)
        arrG=np.full(xnum,0).astype(float)

        arraB=np.full(xnum,0).astype(float)
        arrLB=np.full(xnum,0).astype(float)
        arrdB=np.full(xnum,0).astype(float)
        arrGB=np.full(xnum,0).astype(float)
        lemin=10000
        leminb=10000
        

        for e in range (len(segmentations)-1):
            ya=np.full(xnum,0).astype(float)
            yL=np.full(xnum,0).astype(float)
            yd=np.full(xnum,0).astype(float)
            yaB=np.full(xnum,0).astype(float)
            yLB=np.full(xnum,0).astype(float)
            ydB=np.full(xnum,0).astype(float)
            yG=np.full(xnum,0).astype(float)
            yGB=np.full(xnum,0).astype(float)

            le=int(segmentations[e+1]-segmentations[e]) 
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                    yaB[i]=accX[(segmentations[e]+i)]
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    ydB[i]=dGF[segmentations[e]+i]
                else:
                    ya[i]=accX[segmentations[e]+i]
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    yd[i]=dGF[segmentations[e]+i]       
            if(down):
                ax[0][0].plot(x[:le], yaB[:le],color=(0.8,0.8,0.8))
                ax[1][0].plot(x[:le], yLB[:le],color=(0.8,0.8,0.8))
                ax[2][0].plot(x[:le], yGB[:le],color=(0.8,0.8,0.8))
                ax[3][0].plot(x[:le], ydB[:le],color=(0.8,0.8,0.8))
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                    
                    arraB[i]+=yaB[i]/10
                    arrLB[i]+=yLB[i]/10
                    arrdB[i]+=ydB[i]/10
                    arrGB[i]+=yGB[i]/10
            else:
                ax[0][1].plot(x[:le], ya[:le],color=(0.8,0.8,0.8))
                ax[1][1].plot(x[:le], yL[:le],color=(0.8,0.8,0.8))
                ax[2][1].plot(x[:le], yG[:le],color=(0.8,0.8,0.8))
                ax[3][1].plot(x[:le], yd[:le],color=(0.8,0.8,0.8))
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                    
                    arra[i]+=ya[i]/10
                    arrL[i]+=yL[i]/10
                    arrd[i]+=yd[i]/10
                    arrG[i]+=yG[i]/10
                    
            
                    
        le=leminb
        ax[0][0].plot(x[:le], arraB[:le],color=(0,0,0))
        ax[1][0].plot(x[:le], arrLB[:le],color=(0,0,0))
        ax[2][0].plot(x[:le], arrGB[:le],color=(0,0,0))
        ax[3][0].plot(x[:le], arrdB[:le],color=(0,0,0))     
        
        le=lemin
        ax[0][1].plot(x[:le], arra[:le],color=(0,0,0))
        ax[1][1].plot(x[:le], arrL[:le],color=(0,0,0))
        ax[2][1].plot(x[:le], arrG[:le],color=(0,0,0))
        ax[3][1].plot(x[:le], arrd[:le],color=(0,0,0))    
        
        ax[0][0].set_xlabel("Time [iteration]", fontsize=13)
        ax[0][0].set_ylabel("Acceleration [N/s]", fontsize=13)

        
        ax[1][0].set_xlabel("Time [iteration]", fontsize=13)
        ax[1][0].set_ylabel("LF[N]", fontsize=13)
        
        ax[2][0].set_xlabel("Time [iteration]", fontsize=13)
        ax[2][0].set_ylabel("GF [N]", fontsize=13)
        
        ax[3][0].set_xlabel("Time [iteration]", fontsize=13)
        ax[3][0].set_ylabel("GF derivative [N/s]", fontsize=13)
        plt.show()
        for e in kind:
            if e  in  name:
                fig.savefig("figures/" + e +  "/%s_%dMovements.png" %(s,trial))
    
        
