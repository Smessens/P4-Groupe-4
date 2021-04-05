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
#subjects = ["SophieEB","SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF"] #Names of subjects
subjects = ["JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
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

                    
#        fig = plt.figure(figsize = [15,7])
#        ax  = fig.subplots(4,2)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        ax[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
#        ax[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")
#        
#
#        
#        down=True
#        
#        arra=np.full(xnum,0).astype(float)
#        arrL=np.full(xnum,0).astype(float)
#        arrd=np.full(xnum,0).astype(float)
#        arrG=np.full(xnum,0).astype(float)
#
#        arraB=np.full(xnum,0).astype(float)
#        arrLB=np.full(xnum,0).astype(float)
#        arrdB=np.full(xnum,0).astype(float)
#        arrGB=np.full(xnum,0).astype(float)
#        lemin=10000
#        leminb=10000
#        
#
#        for e in range (len(segmentations)-1):
#            ya=np.full(xnum,0).astype(float)
#            yL=np.full(xnum,0).astype(float)
#            yd=np.full(xnum,0).astype(float)
#            yaB=np.full(xnum,0).astype(float)
#            yLB=np.full(xnum,0).astype(float)
#            ydB=np.full(xnum,0).astype(float)
#            yG=np.full(xnum,0).astype(float)
#            yGB=np.full(xnum,0).astype(float)
#
#            le=int(segmentations[e+1]-segmentations[e]) 
#            for i in range (int(segmentations[e+1]-segmentations[e])):
#                if(down):
#                    yaB[i]=accX[(segmentations[e]+i)]
#                    yLB[i]=LF[(segmentations[e]+i)]
#                    yGB[i]=GF[segmentations[e]+i]
#                    ydB[i]=dGF[segmentations[e]+i]
#                else:
#                    ya[i]=accX[segmentations[e]+i]
#                    yL[i]=LF[segmentations[e]+i]
#                    yG[i]=GF[segmentations[e]+i]
#                    yd[i]=dGF[segmentations[e]+i]       
#            if(down):
#                ax[0][0].plot(x[:le], yaB[:le],color=(0.8,0.8,0.8))
#                ax[1][0].plot(x[:le], yLB[:le],color=(0.8,0.8,0.8))
#                ax[2][0].plot(x[:le], yGB[:le],color=(0.8,0.8,0.8))
#                ax[3][0].plot(x[:le], ydB[:le],color=(0.8,0.8,0.8))
#                down=False
#                if(leminb>le):leminb=le
#                for i in range (le):
#                    
#                    arraB[i]+=yaB[i]/10
#                    arrLB[i]+=yLB[i]/10
#                    arrdB[i]+=ydB[i]/10
#                    arrGB[i]+=yGB[i]/10
#            else:
#                ax[0][1].plot(x[:le], ya[:le],color=(0.8,0.8,0.8))
#                ax[1][1].plot(x[:le], yL[:le],color=(0.8,0.8,0.8))
#                ax[2][1].plot(x[:le], yG[:le],color=(0.8,0.8,0.8))
#                ax[3][1].plot(x[:le], yd[:le],color=(0.8,0.8,0.8))
#                down=True
#                if(lemin>le):lemin=le
#                for i in range (le):
#                    
#                    arra[i]+=ya[i]/10
#                    arrL[i]+=yL[i]/10
#                    arrd[i]+=yd[i]/10
#                    arrG[i]+=yG[i]/10
#                    
#            
#                    
#        le=leminb
#        ax[0][0].plot(x[:le], arraB[:le],color=(0,0,0))
#        ax[1][0].plot(x[:le], arrLB[:le],color=(0,0,0))
#        ax[2][0].plot(x[:le], arrGB[:le],color=(0,0,0))
#        ax[3][0].plot(x[:le], arrdB[:le],color=(0,0,0))     
#        
#        le=lemin
#        ax[0][1].plot(x[:le], arra[:le],color=(0,0,0))
#        ax[1][1].plot(x[:le], arrL[:le],color=(0,0,0))
#        ax[2][1].plot(x[:le], arrG[:le],color=(0,0,0))
#        ax[3][1].plot(x[:le], arrd[:le],color=(0,0,0))    
#        
#        ax[0][0].set_xlabel("Time [iteration]", fontsize=13)
#        ax[0][0].set_ylabel("Acceleration [N/s]", fontsize=13)
#
#        
#        ax[1][0].set_xlabel("Time [iteration]", fontsize=13)
#        ax[1][0].set_ylabel("LF[N]", fontsize=13)
#        
#        ax[2][0].set_xlabel("Time [iteration]", fontsize=13)
#        ax[2][0].set_ylabel("GF [N]", fontsize=13)
#        
#        ax[3][0].set_xlabel("Time [iteration]", fontsize=13)
#        ax[3][0].set_ylabel("GF derivative [N/s]", fontsize=13)
#        plt.show()
#        for e in kind:
#            if e  in  name:
#                fig.savefig("figures/" + e +  "/%s_%dMovements.png" %(s,trial))
#                
#                
                
                
                
                
       #Plot de GF en fonction de LF pour toutes les données de Julian
       
subject_number=0;
fig = plt.figure(figsize = [15,7])
axj = fig.subplots(1,2,squeeze=False)
xnum=8000
x=np.arange(xnum)
axj[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axj[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold") 
colors = ['r','r','b','b','g','g','c','c']
arrm1 = np.full(xnum,0).astype(float)
arrm12 = np.full(xnum,0).astype(float)
arrm2 = np.full(xnum,0).astype(float)
arrm22 = np.full(xnum,0).astype(float)
c = 0
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
               
       
       
        
        #fig.suptitle("%s_00%d" % (s,trial))
         
        
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
            moy1 = np.full(xnum,0).astype(float)
            moy12 = np.full(xnum,0).astype(float)
            moy2 = np.full(xnum,0).astype(float)
            moy22 = np.full(xnum,0).astype(float)
            
            le=int(segmentations[e+1]-segmentations[e]) 
            
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                    
                    yLB[i]=LF[(segmentations[e]+i)]
                    yGB[i]=GF[segmentations[e]+i]
                    fracB[i]=yGB[i]/yLB[i]
                    if(trial == 1):
                        moy1[i]+=fracB[i]/8 
                        
                    else: 
                        moy12[i]+=fracB[i]/8
                    
                else:
                    
                    yL[i]=LF[segmentations[e]+i]
                    yG[i]=GF[segmentations[e]+i]
                    frac[i]=yG[i]/yL[i]
                    if(trial == 1):
                        moy2[i]+=frac[i]/8
                    else: 
                        moy22[i]+=frac[i] 
                                             
            if(down):
                #axj[0][0].plot(x[:le], fracB[:le],color=(0.8,0.8,0.8))
               
                
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                    
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(trial == 1):
                        arrm1[i]+= arrfB[i]/8
                        #print(arrm1[i])
                        
                    else:
                        arrm12[i]+= arrfB[i]/8
                        #print(arrm12[i])
                    
            else:
                
                #axj[0][1].plot(x[:le], frac[:le],color=(0.8,0.8,0.8))
                
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                      
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(trial == 1):
                        arrm2[i]+= arrf[i]/8
                        #print(arrm2[i])
                    else:
                        arrm22[i]+= arrf[i]/8
                        #print(arrm22[i])
                 
le=leminb

        #axj[0][0].plot(x[:le], arrfB[:le],color=colors[c], label = "%s_00%d" % (s,trial))
axj[0][0].plot(x[:le], arrm1[:le],color='r', label = 'essai 1')
axj[0][0].plot(x[:le], arrm12[:le],color='b',label = 'essai 2')
 
le=lemin
        
        #axj[0][1].plot(x[:le], arrf[:le],color = colors[c],label = "%s_00%d" % (s,trial))
axj[0][1].plot(x[:le], arrm2[:le],color = 'r',label = 'essai 1')
axj[0][1].plot(x[:le], arrm22[:le],color = 'b',label = 'essai 2')
        
axj[0][0].set_xlabel("Time [iteration]", fontsize=13)
axj[0][0].set_ylabel("GF/LF", fontsize=13)
        
axj[0][1].set_xlabel("Time [iteration]", fontsize=13)
axj[0][1].set_ylabel("GF/LF", fontsize=13)
c = c+1
        
plt.show()
        #for e in kind:
         #   if e  in  name:
           #     fig.savefig("figures/" + e +  "/%s_%dMovements.png" %(s,trial))
        
        
        
#        
#        fig = plt.figure(figsize = [15,7])
#        ax1  = fig.subplots(len(subjects),2)
#        fig.suptitle("%s_00%d" % (s,trial))
#        xnum=8000
#        x=np.arange(xnum)
#        Le = 0
#        
#        for s in subjects : 
#            if(Le == len(subjects)-1):
#               break     
#            le = leminb
#            ax1[Le][0].plot(arrLB[:le], arrGB[:le],color='g')
#            le = lemin
#            ax1[Le][1].plot(arrLB[:le], arrG[:le],color='m')
#            ax1[Le][0].set_xlabel("LF [N]", fontsize=13)
#            ax1[Le][0].set_ylabel("GF [N]", fontsize=13)
#            Le = Le +1
           
           
           
           
           
           
     
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
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c','b','b','y','c']
fig = plt.figure(figsize = [15,7])
axSECV  = fig.subplots(1,2,squeeze=False)
        #fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axSECV[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axSECV[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
arrm = np.full(xnum,0).astype(float)
arrmB = np.full(xnum,0).astype(float)
flex1 = np.full(xnum,0).astype(float)
flex2 = np.full(xnum,0).astype(float)
c = 0
d = 0
for s in subjects :
    k = 1
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
            fracm = np.full(xnum,0).astype(float)
            fracmB = np.full(xnum,0).astype(float)
                
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
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    arrmB[i]=arrmB[i]+arrfB[i]
                    if(k==2):
                        arrmB[i]= arrmB[i]/2
                        flex1[i]=arrmB[i]
                        arrmB[i]=0
                    
                    
                    
                    
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    arrm[i]+= arrf[i]
                    if(k==2):
                        arrm[i]= arrm[i]/2
                        flex2[i]=arrm[i]
                        arrm[i]=0
                        
                
        if(k%2!=0):        
            k=2
        else:         
            k=1
        le=leminb
        
        axSECV[0][0].plot(x[:le], flex1[:le],color=colors[c])
     
        le=lemin
            
        axSECV[0][1].plot(x[:le], flex2[:le],color=colors[c])
            
        axSECV[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axSECV[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axSECV[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axSECV[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
     
            
plt.show()
    



#Graphe pour SECVR   
subjects = ["julSECVR","SimonSECVR","SophieSECVR","JulianSECVR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']
fig = plt.figure(figsize = [15,7])
axSECVR  = fig.subplots(1,2,squeeze=False)
#fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axSECVR[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axSECVR[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
d=0           
c = 0
for s in subjects :
    k=1
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
                #axSECVR[0][0].plot(x[:le], fracB[:le],color=colors[c],linewidth = 0.2)
                   
                down=False
                if(leminb>le):leminb=le
                for i in range (le):
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(k==2): 
                        arrfB[i]+=arrfB[i]
            else:
                    
                #axSECVR[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(k==2): 
                        arrf[i]+=arrf[i]
            if(k ==2):
                arrfB[i]= arrfB[i]/2
                arrf[i]= arrf[i]/2
            k = k+1
        if(c%2!=0):
            d = c-1 
        else : 
            d = c                 
        le=leminb
        if(c%2==0):
            c = c-1
        axSECVR[0][0].plot(x[:le], arrfB[:le],color=colors[d])
     
        le=lemin
            
        axSECVR[0][1].plot(x[:le], arrf[:le],color=colors[d])
            
        axSECVR[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axSECVR[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axSECVR[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axSECVR[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
            
plt.show()
        
        
#Graphe de SEF        
subjects = ["julSEF","SimonSEF","SophieSEF","JulianSEF"]        
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']
fig = plt.figure(figsize = [15,7])
axSEF  = fig.subplots(1,2,squeeze=False)
       # fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axSEF[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axSEF[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
d = 0            
c = 0
for s in subjects :
    k=1
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
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(k==2): 
                        arrfB[i]+=arrfB[i]
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(k==2): 
                        arrf[i]+=arrf[i]
            if(k ==2):
                arrfB[i]= arrfB[i]/2
                arrf[i]= arrf[i]/2
            k = k+1
        if(c%2!=0):
            d = c-1 
        else : 
            d = c
        le=leminb
    
        axSEF[0][0].plot(x[:le], arrfB[:le],color=colors[d])
     
        le=lemin
            
        axSEF[0][1].plot(x[:le], arrf[:le],color=colors[d])
            
        axSEF[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axSEF[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axSEF[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axSEF[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
            
plt.show()
        
    
#Graphe pour SEFR
subjects = ["julSEFR","SimonSEFR","SophieSEFR","JulianSEFR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']
fig = plt.figure(figsize = [15,7])
axSEFR  = fig.subplots(1,2,squeeze=False)
        #fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axSEFR[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axSEFR[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
c = 0
d=0
for s in subjects :
    k=1
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
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(k==2): 
                        arrfB[i]+=arrfB[i]
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(k==2): 
                        arrf[i]+=arrf[i]
            if(k ==2):
                arrfB[i]= arrfB[i]/2
                arrf[i]= arrf[i]/2
            k = k+1
        if(c%2!=0):
            d = c-1 
        else : 
            d = c            
        le=leminb
    
        axSEFR[0][0].plot(x[:le], arrfB[:le],color=colors[d])
     
        le=lemin
            
        axSEFR[0][1].plot(x[:le], arrf[:le],color=colors[d])
            
        axSEFR[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axSEFR[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axSEFR[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axSEFR[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
            
plt.show()

#Grpahe pour EH
subjects = ["JULEH","SimonEH","SophieEH","JulianEH"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']
fig = plt.figure(figsize = [15,7])
axEH  = fig.subplots(1,2,squeeze=False)
#fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axEH[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axEH[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
d=0            
c = 0
for s in subjects :
    k=1
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
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(k==2): 
                        arrfB[i]+=arrfB[i]
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(k==2): 
                        arrf[i]+=arrf[i]
            if(k ==2):
                arrfB[i]= arrfB[i]/2
                arrf[i]= arrf[i]/2
            k = k+1
        if(c%2!=0):
            d = c-1 
        else : 
            d = c            
        le=leminb
    
        axEH[0][0].plot(x[:le], arrfB[:le],color=colors[d])
     
        le=lemin
            
        axEH[0][1].plot(x[:le], arrf[:le],color=colors[d])
            
        axEH[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axEH[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axEH[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axEH[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
            
plt.show()




#Graphe  pour EHR
subjects = ["JULEHR","SimonEHR","SophieEHR","JulianEHR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']
fig = plt.figure(figsize = [15,7])
axEHR  = fig.subplots(1,2,squeeze=False)
#fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axEHR[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axEHR[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
d=0           
c = 0
for s in subjects :
    k=1
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
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(k==2): 
                        arrfB[i]+=arrfB[i]
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(k==2): 
                        arrf[i]+=arrf[i]
            if(k ==2):
                arrfB[i]= arrfB[i]/2
                arrf[i]= arrf[i]/2
            k = k+1
        if(c%2!=0):
            d = c-1 
        else : 
            d = c             
        le=leminb
    
        axEHR[0][0].plot(x[:le], arrfB[:le],color=colors[d])
     
        le=lemin
            
        axEHR[0][1].plot(x[:le], arrf[:le],color=colors[d])
            
        axEHR[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axEHR[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axEHR[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axEHR[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
            
plt.show()
        
        
#Graphe pour EB
subjects = ["JULEB","SimonEB","SophieEB","JulianEB"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']
fig = plt.figure(figsize = [15,7])
axEB  = fig.subplots(1,2,squeeze=False)
        #fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axEB[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axEB[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
c = 0
d=0
for s in subjects :
    k=1
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
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(k==2): 
                        arrfB[i]+=arrfB[i]
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(k==2): 
                        arrf[i]+=arrf[i]
            if(k ==2):
                arrfB[i]= arrfB[i]/2
                arrf[i]= arrf[i]/2
            k = k+1
        if(c%2!=0):
            d = c-1 
        else : 
            d = c            
        le=leminb
    
        axEB[0][0].plot(x[:le], arrfB[:le],color=colors[d])
     
        le=lemin
            
        axEB[0][1].plot(x[:le], arrf[:le],color=colors[d])
            
        axEB[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axEB[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axEB[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axEB[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
            
plt.show()

#Graphe pour EBR
subjects = ["JULEBR","SimonEBR","SophieEBR","JulianEBR"]
ntrials = 2 
kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
colors = ['r','b','g','c']
fig = plt.figure(figsize = [15,7])
axEBR  = fig.subplots(1,2,squeeze=False)
fig.suptitle("%s_00%d" % (s,trial))
xnum=8000
x=np.arange(xnum)
axEBR[0][0].set_title("GML data per down movements", fontsize=14, fontweight="bold")
axEBR[0][1].set_title("GML data per up movements", fontsize=14, fontweight="bold")  
d=0            
c = 0
for s in subjects :
    k=1
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
                        
                    arrLB[i]+=yLB[i]/10
                    arrGB[i]+=yGB[i]/10
                    arrfB[i]= arrGB[i]/arrLB[i]
                    if(k==2): 
                        arrfB[i]+=arrfB[i]
            else:
                    
                #axSECV[0][1].plot(x[:le], frac[:le],color=colors[c],linewidth = 0.2)
                down=True
                if(lemin>le):lemin=le
                for i in range (le):
                          
                    arrL[i]+=yL[i]/10
                    arrG[i]+=yG[i]/10
                    arrf[i]= arrG[i]/arrL[i]
                    if(k==2): 
                        arrf[i]+=arrf[i]
            if(k ==2):
                arrfB[i]= arrfB[i]/2
                arrf[i]= arrf[i]/2
            k = k+1
        if(c%2!=0):
            d = c-1 
        else : 
            d = c            
        le=leminb
    
        axEBR[0][0].plot(x[:le], arrfB[:le],color=colors[d])
     
        le=lemin
            
        axEBR[0][1].plot(x[:le], arrf[:le],color=colors[d])
            
        axEBR[0][0].set_xlabel("Time [iteration]", fontsize=13)
        axEBR[0][0].set_ylabel("GF/LF", fontsize=13)
           
        axEBR[0][1].set_xlabel("Time [iteration]", fontsize=13)
        axEBR[0][1].set_ylabel("GF/LF", fontsize=13)
        c=c+1
plt.show()