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

subjects =["Julian"]

kind=["EB_","EH_","SECV_","SEF_","EBR_","EHR_","SECVR_","SEFR_"]#"EH","EB","EHR","EBR"]#
expe=["SECV","SEF","EH","EB","SECVR","SEFR","EHR","EBR"]#,"EH","SECV","SEF","EBR","EHR","SECVR","SEFR"]

ntrials = 2 #Number of trials for each subject

master={"EB":np.zeros((400,600)),"EBcount":0,"EBD":np.zeros((400,600)),"EBDcount":0,
        "EBR":np.zeros((400,600)),"EBRcount":0,"EBRD":np.zeros((400,600)),"EBRDcount":0,
        "EH":np.zeros((400,600)),"EHcount":0,"EHD":np.zeros((400,600)),"EHDcount":0,
        "EHR":np.zeros((400,600)),"EHRcount":0,"EHRD":np.zeros((400,600)),"EHRDcount":0,
        "SECV":np.zeros((400,600)),"SECVcount":0,"SECVD":np.zeros((400,600)),"SECVDcount":0,
        "SECVR":np.zeros((400,600)),"SECVRcount":0,"SECVRD":np.zeros((400,600)),"SECVRDcount":0,
        "SEF":np.zeros((400,600)),"SEFcount":0,"SEFD":np.zeros((400,600)),"SEFDcount":0,
        "SEFR":np.zeros((400,600)),"SEFRcount":0,"SEFRD":np.zeros((400,600)),"SEFRDcount":0}


# Double for-loop that runs thrgough all subjects and trials
subject_number=0;
for s in subjects:
    for exp in expe:
        for trial in range(1,ntrials+1): 

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

                            


                    
   
fig = plt.figure(figsize = [8,8])
ax  = fig.subplots(1,1)
fig.suptitle("")
data=np.zeros((18*len(subjects),len(expe)))



for j in range (len(expe)):
    name = expe[j]+"_"
    
    xnum=10000
    x=np.arange(0,xnum/800,1/800)
    le=200 
                  
    arr=np.full(master[expe[j]+"count"],0).astype(float)
        #amplitude
    for  i in range (min(master[expe[j]+"count"],18*len(subjects))):     
        data[i][j]=abs(master[expe[j]][i][0]-master[expe[j]][i][le])
        
data[:][:]-=np.min(data[:][:])
print(len(data[:][:]))   






ax.set_ylabel("Amplitude [cm]", fontsize=13)
ax.set_ylim(0,0.68)
ax.set_title("Répartition de l'amplitude pour les conditions élastique")
ax.axhline(y=0.45, color='b', linestyle=':')
ax.grid( linestyle=':', linewidth=2)


ax.violinplot(data)


# fill with colors
color={"SECV":"grey","SEF":"lightpink","EH":"red","EB":"gold","SECVR":"lightgreen","SEFR":"turquoise","EHR":"mediumpurple","EBR":"deeppink"}
#color={"SECV":"green","SEF":"red","EH":"blue","EB":"purple","SECVR":"black","SEFR":"orange","EHR":"yellow","EBR":"pink"}
#color={"SECV":"green","SEF":"green","EH":"green","EB":"green","SECVR":"blue","SEFR":"blue","EHR":"blue","EBR":"blue"}
#color={"SECV":"green","SEF":"green","EH":"green","EB":"blue","SECVR":"blue","SEFR":"blue","EHR":"green","EBR":"blue"}
"""

b0 = ax.boxplot(np.transpose(dataD),
                   vert=True,  # vertical box alignment
                   patch_artist=True,  # fill with color
                   labels=expe)

colors = [color[i] for i in expe]

for bplot in ([b0]):
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

# adding horizontal grid lines
for ax in [ax]:
    ax.yaxis.grid(True)
    ax.set_xlabel('Conditions')
    ax.set_ylabel('Amplitudes [cm]')
"""
plt.show()


    

    
            
