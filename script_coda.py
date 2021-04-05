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
#subjects = ["SophieEB"]#,"SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF"] #Names of subjects
#subjects = ["JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
#subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","SophieEB","SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF","JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
subjects =["Sophie"]#"Julien","Julian","Simon","Sophie"]

kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
expe=["SECV","SEF","EH","EB","SECVR","SEFR","EHR","EBR"]#,"EH","SECV","SEF","EBR","EHR","SECVR","SEFR"]

ntrials = 2 #Number of trials for each subject

master={"EB":np.zeros((80,600)),"EBcount":0,"EBD":np.zeros((80,600)),"EBDcount":0,
        "EBR":np.zeros((80,600)),"EBRcount":0,"EBRD":np.zeros((80,600)),"EBRDcount":0,
        "EH":np.zeros((80,600)),"EHcount":0,"EHD":np.zeros((80,600)),"EHDcount":0,
        "EHR":np.zeros((80,600)),"EHRcount":0,"EHRD":np.zeros((80,600)),"EHRDcount":0,
        "SECV":np.zeros((80,600)),"SECVcount":0,"SECVD":np.zeros((80,600)),"SECVDcount":0,
        "SECVR":np.zeros((80,600)),"SECVRcount":0,"SECVRD":np.zeros((80,600)),"SECVRDcount":0,
        "SEF":np.zeros((80,600)),"SEFcount":0,"SEFD":np.zeros((80,600)),"SEFDcount":0,
        "SEFR":np.zeros((80,600)),"SEFRcount":0,"SEFRD":np.zeros((80,600)),"SEFRDcount":0}


# Double for-loop that runs thrgough all subjects and trials
subject_number=0;
for s in subjects:
    for exp in expe:
        for trial in range(1,ntrials+1): 
            if("%s%s_00%d" % (s,exp,trial)!="SimonSECV_002" 
               and "%s%s_00%d" % (s,exp,trial)!="JulienSECV_001"
               and "%s%s_00%d" % (s,exp,trial)!="JulianSEF_001"
               and "%s%s_00%d" % (s,exp,trial)!="SimonEB_002" ):#coumés
                # Set data path
                glm_path = "DataGroupe4/%s%s_00%d.glm" % (s,exp,trial)
                coda_path = "Groupe_4_codas/%s%s_000%d.txt" % (s,exp,trial)
                name ="%s%s_00%d.glm" % (s,exp,trial)
        
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
                    yx=np.full(xnum,0).astype(float)
        
                    yxB=np.full(xnum,0).astype(float)
        
        
                    le=int(segmentations[e+1]-segmentations[e]) 
                    for i in range (int(segmentations[e+1]-segmentations[e])):
                        if(down):
                            yxB[i]=pos[0][(segmentations[e]+i)]
        
                        else:
                            yx[i]=pos[0][segmentations[e]+i]
        
                    if(down):
                        master[exp+"D"][master[exp+"Dcount"]]=yxB[:600]
                        master[exp+"Dcount"]+=1
        
                        down=False
    
                        if(leminb>le):leminb=le
        
        
        
                            
                    else:
                        master[exp][master[exp+"count"]]=yx[:600]
                        master[exp+"count"]+=1
                        down=True
                        if(lemin>le):lemin=le

                            


                    
"""                    
for exp in expe:
    name = exp+"_"

    fig = plt.figure(figsize = [15,7])
    ax  = fig.subplots(2)
    fig.suptitle(exp)
    xnum=10000
    x=np.arange(0,xnum/800,1/800)
    le=200                    
    ax[0].set_title("", fontsize=14, fontweight="bold")
    ax[1].set_title("", fontsize=14, fontweight="bold")
    ax[0].set_ylim([-1,0])
    ax[1].set_ylim([-1,0])

    
        
    arrxB=np.full(xnum,0).astype(float)
    arrxstdB=np.full(xnum,0).astype(float)
    for  i in range (600):   
        y=[row[i] for row in  master[exp+"D"]]
        y=y[:master[exp+"Dcount"]]
        arrxB[i]=np.mean(y)
        #arrxstdB[i]=np.std(y)
    
    for col in master[exp+"D"]:
        if(np.sum(col)!=0):
            ax[0].plot(x[:le], col[:le],color=(0.9,0.9,0.9))
     
        
    arrx=np.full(xnum,0).astype(float)
    arrxstd=np.full(xnum,0).astype(float)
    for  i in range (600):   
        y=[row[i] for row in  master[exp]]
        y=y[:master[exp+"count"]]
        arrx[i]=np.mean(y)
       # arrxstd[i]=np.std(y)
        
        
    for col in  master[exp]:
        if(np.sum(col)!=0):
            ax[1].plot(x[:le], col[:le],color=(0.9,0.9,0.9))
            
            
            
    ax[0].plot(x[:le], arrxB[:le],color=(0,0,0))

    #ax[0][0].plot(x[:le], arrxstdB[:le],color=(1,0,0))
    
    ax[1].plot(x[:le], arrx[:le],color=(0,0,0))
    #ax[0][1].plot(x[:le], arrxstd[:le],color=(1,0,0))
    
    
    ax[0].set_xlabel("Time [s]", fontsize=13)
    ax[0].set_ylabel("Pos X [m]", fontsize=13)
    
    
    
    plt.show()

    fig.savefig("figures/" + exp +"_/%s_mean_Coda_Movements.png" %(exp))
        
  
"""    
fig = plt.figure(figsize = [20,9])
ax  = fig.subplots(2,8)
fig.suptitle("All Movements Sophie")
ax[0][0].set_ylabel("Pos X [m]", fontsize=13)

for j in range (len(expe)):
    name = expe[j]+"_"
    
    xnum=10000
    x=np.arange(0,xnum/800,1/800)
    le=200 
                   
    ax[0][j].set_title(expe[j], fontsize=14, fontweight="bold")
    ax[0][j].set_ylim([-1,0])
    ax[1][j].set_ylim([-1,0])
    
    
    
    arrxB=np.full(xnum,0).astype(float)
    arrxstdB=np.full(xnum,0).astype(float)
    for  i in range (600):   
        y=[row[i] for row in  master[expe[j]+"D"]]
        y=y[:master[expe[j]+"Dcount"]]
        arrxB[i]=np.mean(y)
        #arrxstdB[i]=np.std(y)
    
    for col in master[expe[j]+"D"]:
        if(np.sum(col)!=0):
            ax[0][j].plot(x[:le], col[:le],color=(0.9,0.9,0.9))
     
        
    arrx=np.full(xnum,0).astype(float)
    arrxstd=np.full(xnum,0).astype(float)
    for  i in range (600):   
        y=[row[i] for row in  master[expe[j]]]
        y=y[:master[expe[j]+"count"]]
        arrx[i]=np.mean(y)
       # arrxstd[i]=np.std(y)
        
        
    for col in  master[expe[j]]:
        if(np.sum(col)!=0):
            ax[1][j].plot(x[:le], col[:le],color=(0.9,0.9,0.9))
            
            
            
    ax[0][j].plot(x[:le], arrxB[:le],color=(0,0,0))
    
    #ax[0][0].plot(x[:le], arrxstdB[:le],color=(1,0,0))
    
    ax[1][j].plot(x[:le], arrx[:le],color=(0,0,0))
    #ax[0][1].plot(x[:le], arrxstd[:le],color=(1,0,0))
    
    
    ax[0][j].set_xlabel("Time [s]", fontsize=13)

    
plt.show()

fig.savefig("figures/MasterMovements")
    
            