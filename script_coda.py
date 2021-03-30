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





def filter_signal(y, axis=1, fs=200, fc=10, N=4, type='low'):
    """Filters signal y by using a Butterworth filter of order N and a cut-off 
    frequency fc."""
    
    # Converts the cut-off frequency to [pi rad/s]
    Wn = fc/(fs/2)
    
    
    # Create butterworth digital filter
    b,a = signal.butter(N,Wn,btype=type,analog=False)
    
    # Filter y with a zero-phase forward and reverse digital IIR
    ys = signal.filtfilt(b,a,y,axis=axis)
    
    return ys



def get_segmentations(s,v):
    file_path = "Groupe_4_codas/%s_000%d.txt" % (s,v)
    
    coda_df = coda.import_data(file_path)


    seg=np.array([])
    #%% Compute the coordinates of the center of the manipulandum 
    
    # Markers of the manipulandum (order matters! See doc of manipulandum_center)
    markers_id = [6,7,8,9] 

    # Center position
    pos = coda.manipulandum_center(coda_df,markers_id)
    pos = pos/1000;

    # Filter position signal
    #pos = filter_signal(pos,axis=1,fs=200,fc=10,N=4)
    

    
    pos[0]-=min(pos[0])
    ma=max(pos[0])
    up=True
    for i in range(len(pos[1])):
        percent=(pos[0][i]/ma)
    
        if(up):
            if(percent<=0.4):
                seg=np.append(seg,int(i*4))

                up=False
        elif(percent>=0.6):
            seg=np.append(seg,int(i*4))
            up=True
            
    return seg
    

if (__name__ == "__main__"):
     #%% Import data from the txt file
    file_path = "DataGroupe4/JULEB_001.txt"
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
        
        
    
    arr=pos[0]-min(pos[0])
    ma=max(arr)
    up=True
    for i in range(len(arr)):
        percent=(arr[i]/ma)
    
        if(up):
            if(percent<=0.5):
                ax[0,0].axvline(x=time[i])
                ax[0,1].axvline(x=time[i])
                up=False
        elif(percent>=0.6):
            ax[0,0].axvline(x=time[i])
            ax[0,1].axvline(x=time[i])
            up=True
        
    
    seg=np.array([])


    pos[0]-=min(pos[0])
    ma=max(pos[0])
    up=True
    for i in range(len(pos[1])):
        percent=(pos[0][i]/ma)
    
        if(up):
            if(percent<=0.4):
                seg=np.append(seg,i)

                up=False
        elif(percent>=0.6):
            seg=np.append(seg,i)
            up=True
            
            


    

        
        
        
    ax[0,0].set_ylabel("X pos [m]", fontsize=12)
    ax[1,0].set_ylabel("Y pos [m]", fontsize=12)
    ax[2,0].set_ylabel("Z pos [m]", fontsize=12)
    
    ax[0,1].set_ylabel("X vel [m/s]", fontsize=12)
    ax[1,1].set_ylabel("Y vel [m/s]", fontsize=12)
    ax[2,1].set_ylabel("Z vel [m/s]", fontsize=12)
        


 #%% Import data from the txt file
    file_path = "DataGroupe4/JULEB_002.txt"
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

    
 




