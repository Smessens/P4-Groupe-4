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

mast=["SECV","SEF","EH","EB","SECVR","SEFR","EHR","EBR"]#,"EH","SECV","SEF","EBR","EHR","SECVR","SEFR"]

fig = plt.figure(figsize = [20,9])
ax  = fig.subplots(2)
fig.suptitle("Moyenne position en x Julien")

for m in range (len( mast)):
     
    #subjects = ["JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF"] #Names of subjects
    #subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF"] #Names of subjects
    #Notes  not existant: EB2, EHR
    #subjects = ["SophieEB"]#,"SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF"] #Names of subjects
    #subjects = ["JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
    #subjects = ["SimonEH","SimonSECV","SimonSECVR","SimonSEFR","SimonSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","JULEB","JULEBR","JULEH","JULEHR","JulSECV","julSECVR","julSEFR","julSEF","SophieEB","SophieEBR","SophieEH","SophieEHR","SophieSECV","SophieSECVR","SophieSEFR","SophieSEF","JulianEB","JulianEBR","JulianEH","JulianEHR","JulianSECV","JulianSECVR","JulianSEFR","JulianSEF"] #Names of subjects
    subjects =["Julien"]#"Julien","Julian","Simon","Sophie"]
    
    kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
    expe=[mast[m]]#"SECV","SEF","EH","EB","SECVR","SEFR","EHR","EBR"]#,"EH","SECV","SEF","EBR","EHR","SECVR","SEFR"]
    
    ntrials = 2 #Number of trials for each subject
    
    master={"Sophie":np.zeros((160,600)),"Sophiecount":0,"SophieD":np.zeros((160,600)),"SophieDcount":0,
            "Simon":np.zeros((160,600)),"Simoncount":0,"SimonD":np.zeros((160,600)),"SimonDcount":0,
            "Julian":np.zeros((160,600)),"Juliancount":0,"JulianD":np.zeros((160,600)),"JulianDcount":0,
            "Julien":np.zeros((160,600)),"Juliencount":0,"JulienD":np.zeros((160,600)),"JulienDcount":0}
    
    
    
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
                            master[s+"D"][master[s+"Dcount"]]=yxB[:600]
                            master[s+"Dcount"]+=1
            
                            down=False
        
                            if(leminb>le):leminb=le
            
            
            
                                
                        else:
                            master[s][master[s+"count"]]=yx[:600]
                            master[s+"count"]+=1
                            down=True
                            if(lemin>le):lemin=le
    


    
    ax[0].set_ylabel("Pos X [m]", fontsize=13)
    
    
    #color={"Simon":"green","Julien":"red","Julian":"blue","Sophie":"purple"}
    color={"SECV":"green","SEF":"red","EH":"blue","EB":"purple","SECVR":"black","SEFR":"purple","EHR":"yellow","EBR":"pink"}
    for s in (subjects):
        
        xnum=10000
        x=np.arange(0,xnum/800,1/800)
        le=200 
                       
       # ax[0].set_title(, fontsize=14, fontweight="bold")
        ax[0].set_ylim([-1,0])
        ax[1].set_ylim([-1,0])
        
        
        
        arrxB=np.full(xnum,0).astype(float)
        arrxstdB=np.full(xnum,0).astype(float)
        for  i in range (600):   
            y=[row[i] for row in  master[s+"D"]]
            y=y[:master[s+"Dcount"]]
            arrxB[i]=np.mean(y)
            #arrxstdB[i]=np.std(y)
        
    
        arrx=np.full(xnum,0).astype(float)
        arrxstd=np.full(xnum,0).astype(float)
        for  i in range (600):   
            y=[row[i] for row in  master[s]]
            y=y[:master[s+"count"]]
            arrx[i]=np.mean(y)
           # arrxstd[i]=np.std(y)
           
       
                
        ax[0].plot(x[:le], arrxB[:le],color=color[mast[m]],label=mast[m])

        
        ax[1].plot(x[:le], arrx[:le],color=color[mast[m]],label=mast[m])
        
        ax[0].set_xlabel("Time [s]", fontsize=13)
    
    ax[0].set_title("Mouvement vers le bas", fontsize=11, fontweight="bold")
    ax[0].legend()
    
    ax[1].set_title("Mouvement vers le haut", fontsize=11, fontweight="bold")
    ax[1].legend()


plt.show()

fig.savefig("figures/Moyenne_Tout_sujets")
        
                