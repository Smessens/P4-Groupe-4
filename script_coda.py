# -*- coding: utf-8 -*-
"""
Example script for processing and plotting CODA data

Created on Wed Mar 7  2021

@author: opsomerl
"""
import matplotlib.pyplot as plt

import coda_tools as coda
import processing_tools as tool
import numpy  as np

from scipy import signal

import glm_data_processing as glm
import derive as der

# Fermeture des figures ouvertes
plt.close('all')

#subjects = ["JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF"] #Names of subjects
#subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF"] #Names of subjects
#Notes  not existant: EB2, EHR
subjects = ["SophieEB"]#,"SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF"] #Names of subjects
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
        coda_path = "Groupe_4_codas/%s_000%d.txt" % (s,trial)
        name ="%s_00%d.glm" % (s,trial)

        #%% Import data from the files
        coda_df = coda.import_data(coda_path)
        glm_df = glm.import_data(glm_path)

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
                segmentations=np.append(segmentations,int((i+start)//4))
                bip=True
            elif(bip and (mB0[i]+mB1[i]+mB2[i])==3 ):
                bip=False
            
        segmentations=segmentations.astype(int)
    

        
        
        #%% Compute the coordinates of the center of the manipulandum 
        
        # Markers of the manipulandum (order matters! See doc of manipulandum_center)
        markers_id = [6,7,8,9] 
        
        # Center position
        pos = coda.manipulandum_center(coda_df,markers_id)
        pos = pos/1000;
        
        # Filter position signal
        #pos = tool.filter_signal(pos,axis=1,fs=200,fc=10,N=4)
        
        # Derive position to get velocity
        vel = tool.derive(pos,200,axis=1)
    
        
        #%% Basic plot of the data
        fig = plt.figure(figsize = [15,7])
        ax  = fig.subplots(3,2)
        fig.suptitle("%s_00%d" % (s,trial))
        xnum=10000
        x=np.arange(0,xnum/800,1/800)

        ax[0][0].set_title("Coda data per down movements", fontsize=14, fontweight="bold")
        ax[0][1].set_title("Coda data per up movements", fontsize=14, fontweight="bold")
        

        
        down=True
        
        arrx=np.full(xnum,0).astype(float)
        arry=np.full(xnum,0).astype(float)
        arrz=np.full(xnum,0).astype(float)

        arrxB=np.full(xnum,0).astype(float)
        arryB=np.full(xnum,0).astype(float)
        arrzB=np.full(xnum,0).astype(float)
        
        lemin =10000
        leminb=10000
        
        count=0

        for e in range (len(segmentations)-1):
            yx=np.full(xnum,0).astype(float)
            yy=np.full(xnum,0).astype(float)
            yz=np.full(xnum,0).astype(float)
            yxB=np.full(xnum,0).astype(float)
            yyB=np.full(xnum,0).astype(float)
            yzB=np.full(xnum,0).astype(float)

            le=int(segmentations[e+1]-segmentations[e]) 
            for i in range (int(segmentations[e+1]-segmentations[e])):
                if(down):
                    yxB[i]=pos[0][(segmentations[e]+i)]
                    yyB[i]=pos[1][(segmentations[e]+i)]
                    yzB[i]=pos[2][segmentations[e]+i]
                else:
                    yx[i]=pos[0][segmentations[e]+i]
                    yy[i]=pos[1][segmentations[e]+i]
                    yz[i]=pos[2][segmentations[e]+i]
            if(down):
                ax[0][0].plot(x[:le], yxB[:le],color=(0.8,0.8,0.8))
                ax[1][0].plot(x[:le], yyB[:le],color=(0.8,0.8,0.8))
                ax[2][0].plot(x[:le], yzB[:le],color=(0.8,0.8,0.8))
                down=False
                if(leminb>le):
                    leminb=le

-                for i in range (le):
                    arrxB[i]+=(yxB[i]/10.0)
                    arryB[i]+=(yyB[i]/10.0)
                    arrzB[i]+=(yzB[i]/10.0)
                    
            else:
                ax[0][1].plot(x[:le], yx[:le],color=(0.8,0.8,0.8))
                ax[1][1].plot(x[:le], yy[:le],color=(0.8,0.8,0.8))
                ax[2][1].plot(x[:le], yz[:le],color=(0.8,0.8,0.8))
                down=True
                if(lemin>le):
                    lemin=le
                    
                for i in range (le):
                    arrx[i]+=(yx[i]/10.0)
                    arry[i]+=yy[i]/10.0
                    arrz[i]+=yz[i]/10.0
                    
        
            
                    
        le=leminb
        ax[0][0].plot(x[:le], arrxB[:le],color=(0,0,0))
        ax[1][0].plot(x[:le], arryB[:le],color=(0,0,0))
        ax[2][0].plot(x[:le], arrzB[:le],color=(0,0,0))     
        
        le=lemin
        ax[0][1].plot(x[:le], arrx[:le],color=(0,0,0))
        ax[1][1].plot(x[:le], arry[:le],color=(0,0,0))
        ax[2][1].plot(x[:le], arrz[:le],color=(0,0,0))    
        
        ax[0][0].set_xlabel("Time [s]", fontsize=13)
        ax[0][0].set_ylabel("Pos X [m]", fontsize=13)

        
        ax[1][0].set_xlabel("Time [s]", fontsize=13)
        ax[1][0].set_ylabel("Pos y[m]", fontsize=13)
        
        ax[2][0].set_xlabel("Time [iteration]", fontsize=13)
        ax[2][0].set_ylabel("Pos [m]", fontsize=13)
        

        plt.show()
        for e in kind:
            if e  in  name:
                fig.savefig("figures/" + e +  "/%s_%dCoda_Movements.png" %(s,trial))
    

    