# -*- coding: utf-8 -*-
"""
Example script for processing and plotting CODA data

Created on Wed Mar 7  2021

@author: opsomerl
"""
import matplotlib.pyplot as plt

import coda_tools as coda
import processing_tools as tool


#%% Import data from the txt file
file_path = "DataGroupe4/JulEB_001.txt"
coda_df = coda.import_data(file_path)


#%% Compute the coordinates of the center of the manipulandum 

# Markers of the manipulandum (order matters! See doc of manipulandum_center)
markers_id = [6,7,8,9] 

# Center position
pos = coda.manipulandum_center(coda_df,markers_id)
pos = pos/1000;

# Filter position signal
pos = tool.filter_signal(pos,axis=1,fs=200,fc=10,N=4)

# Derive position to get velocity
vel = tool.derive(pos,200,axis=1)


#%% Basic plot of the data
fig = plt.figure(figsize = [15,7])
ax  = fig.subplots(3,2)

time = coda_df.time.to_numpy()
for i in range(3):
    ax[i,0].plot(time,pos[i])
    ax[i,1].plot(time,vel[i])
    
    

pos[0]-=min(pos[0])
ma=max(pos[0])
up=True
for i in range(len(pos[1])):
    percent=(pos[0][i]/ma)

    if(up):
        if(percent<=0.4):
            ax[0,0].axvline(x=time[i])
            ax[0,1].axvline(x=time[i])
            up=False
    elif(percent>=0.6):
        ax[0,0].axvline(x=time[i])
        ax[0,1].axvline(x=time[i])
        up=True
    

    
    
    
ax[0,0].set_ylabel("X pos [m]", fontsize=12)
ax[1,0].set_ylabel("Y pos [m]", fontsize=12)
ax[2,0].set_ylabel("Z pos [m]", fontsize=12)

ax[0,1].set_ylabel("X vel [m/s]", fontsize=12)
ax[1,1].set_ylabel("Y vel [m/s]", fontsize=12)
ax[2,1].set_ylabel("Z vel [m/s]", fontsize=12)
    
    





