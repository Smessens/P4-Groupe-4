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

subjects =["Julian","Sophie"]

kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]
expe=["EB","EH","SECV","SEF","EBR","EHR","SECVR","SEFR"]

ntrials = 2 #Number of trials for each subject

master={"EB":np.zeros((160,600)),"EBcount":0,"EBD":np.zeros((160,600)),"EBDcount":0,
        "EBR":np.zeros((160,600)),"EBRcount":0,"EBRD":np.zeros((160,600)),"EBRDcount":0,
        "EH":np.zeros((160,600)),"EHcount":0,"EHD":np.zeros((160,600)),"EHDcount":0,
        "EHR":np.zeros((160,600)),"EHRcount":0,"EHRD":np.zeros((160,600)),"EHRDcount":0,
        "SECV":np.zeros((160,600)),"SECVcount":0,"SECVD":np.zeros((160,600)),"SECVDcount":0,
        "SECVR":np.zeros((160,600)),"SECVRcount":0,"SECVRD":np.zeros((160,600)),"SECVRDcount":0,
        "SEF":np.zeros((160,600)),"SEFcount":0,"SEFD":np.zeros((160,600)),"SEFDcount":0,
        "SEFR":np.zeros((160,600)),"SEFRcount":0,"SEFRD":np.zeros((160,600)),"SEFRDcount":0}


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

                            


                    
   
fig = plt.figure(figsize = [9,9])
ax  = fig.subplots(1)

data=np.zeros((2,3))

#exp1=["EH","EHR"]
#exp2=["EB","EBR"]

#exp1=["EB","EHR"]
#exp2=["EH","EBR"]

exp1=["EH","EB","SECV","SEF"]
exp2=["EHR","EBR","SECVR","SEFR"]



e1=np.array([])

for j in range (len(exp1)):
    
    name = exp1[j]+"_"
    
    xnum=10000
    x=np.arange(0,xnum/800,1/800)
    le=200 
                  
    arr=np.full(master[exp1[j]+"Dcount"],0).astype(float)
  
    

    #amplitude
    for  i in range (master[exp1[j]+"Dcount"]): 
    
        arr[i]=abs(master[exp1[j]+"D"][i][0]-master[exp1[j]+"D"][i][le])
        
    e1=np.append(e1,arr)
 
data[0]=[np.mean(e1),np.max(e1),np.min(e1)]


e2=np.array([])

for j in range (len(exp2)):
    
    name = exp2[j]+"_"
    
    xnum=10000
    x=np.arange(0,xnum/800,1/800)
    le=200 
                  
    arr=np.full(master[exp2[j]+"Dcount"],0).astype(float)
  
    

    #amplitude
    for  i in range (master[exp2[j]+"Dcount"]): 
    
        arr[i]=abs(master[exp2[j]+"D"][i][0]-master[exp2[j]+"D"][i][le])
        
    e2=np.append(e2,arr)
    
    

print(np.mean(e1),np.mean(e2))
data[1]=[np.mean(e2),np.max(e2),np.min(e2)]

    
   
    
    
ax.set_ylabel("Amplitude [m]", fontsize=13)
ax.set_ylim(0,0.75)


b0 = ax.boxplot(np.transpose(data),
                   vert=True,  # vertical box alignment
                   patch_artist=True,
                   labels=["Sujet droit","Sujet retourné"])

    

# fill with colors
color={"SECV":"green","SEF":"red","EH":"blue","EB":"purple","SECVR":"black","SEFR":"orange","EHR":"yellow","EBR":"pink"}

colors = [color[i] for i in expe]
ax.axhline(y=0.45, color='b', linestyle=':')

for bplot in ([b0]):
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

# adding horizontal grid lines
for ax in [ax]:
    ax.yaxis.grid(True)
    ax.set_xlabel('condition')
    ax.set_ylabel('Amplitudes [m]')

[plt.axvline(x, color = 'r', linestyle='--') for x in [1.5]]

plt.show()
    

    
            