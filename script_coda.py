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

subjects =["Julian"]#"Julien","Julian","Simon","Sophie"]

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

                            


                    
   
fig = plt.figure(figsize = [20,9])
ax  = fig.subplots(2)
fig.suptitle("Boxplots")
data=np.zeros((len(expe),3))
dataD=np.zeros((len(expe),3))

for j in range (len(expe)):
    name = expe[j]+"_"
    
    xnum=10000
    x=np.arange(0,xnum/800,1/800)
    le=200 
                  
    arr=np.full(master[expe[j]+"Dcount"],0).astype(float)
  
    

    #amplitude
    for  i in range (master[expe[j]+"Dcount"]): 
    
        arr[i]=abs(master[expe[j]+"D"][i][0]-master[expe[j]+"D"][i][le])

    dataD[j]=[np.mean(arr),np.max(arr),np.min(arr)]



    arr=np.full(master[expe[j]+"count"],0).astype(float)
        #amplitude
    for  i in range (master[expe[j]+"count"]): 
    
        arr[i]=abs(master[expe[j]][i][0]-master[expe[j]][i][le])

    data[j]=[np.mean(arr),np.max(arr),np.min(arr)]

    
    
ax[0].set_ylabel("Amplitude [cm]", fontsize=13)
ax[1].set_ylabel("Amplitude [cm]", fontsize=13)

ax[0].set_title(" Mouvement bas")
ax[1].set_title(" Mouvement Haut")

ax[0].boxplot(np.transpose(dataD),labels=expe)
ax[1].boxplot(np.transpose(data), labels=expe)

    
plt.show()

    
            