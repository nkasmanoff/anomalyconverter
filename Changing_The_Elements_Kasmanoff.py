

import numpy as np
import math as mth


# In[24]:

def length(v):
    return mth.sqrt(np.dot(v, v))

def angle(v1, v2):

    return mth.acos(np.dot(v1, v2) / (length(v1) * length(v2)))

def normalize(x):
    xnorm = []
    mag = 0
    for i in range(len(x)):
        mag += x[i]**2    
    for j in range(len(x)):
        xnorm.append(x[j]/np.sqrt(mag))
        
    return xnorm


# In[75]:

def cart2OE(x,y,z,vx,vy,vz):
    #first the things that will be used later... 
    GM=1
    r=np.array([x,y,z])
    rmag = length(r)
    print(rmag)

    v=np.array([vx,vy,vz])
    vmag = length(v)
    h=np.cross(r,v)
    I = np.array([1,0,0])  #x axis vector 
    K = np.array([0,0,1])  #z axis of this setup
    n = np.cross(K,h) #node vector, should be zero when uninclined
    if length(n)==0:
        n=I
    nmag = length(n)  #magnitude of node vector, use this to decide if equatorial orbit or not. 
    #a
    a = 1/((2/rmag)-((vmag**2)/GM))
    #e
    evec=np.cross(v,h)/GM-normalize(r) 
    e=length(evec)
    #i
    hmag = length(h)
    hnorm = normalize(h)
    if z==0 and vz==0:
        i=0
    else:
        i=np.arccos(hnorm[2])  #normalized w z axis, distance btw... seems to pass every test so far!
    #O
    if i>0 and i<np.pi:
        if n[1]>=0:
            O = np.arccos(np.dot(n,I)/length(n))
        elif n[1]<0:
            O = 2*np.pi - np.arccos(np.dot(n,I)/length(n))
    else:
        O = 0
    #w
    if e==0:
        w=0
    else:
        if i > 0 and i < np.pi: 
            if evec[2]>=0:
                w = np.arccos(np.dot(n,evec)/(length(evec)*length(n)))
            elif evec[2]<0:
                w = 2*np.pi - np.arccos(np.dot(n,evec)/(length(evec)*length(n)))
        #elif i==np.pi/2:
        
        else:
            if h[2]>=0:
                w = 2*np.pi + np.arctan2(evec[1],evec[0])
                if w >= 2*np.pi:
                    w = w - 2*np.pi
            else: 
                w = 2*np.pi - np.arctan2(evec[1],evec[0])
                if w >=2*np.pi:
                    w = w - 2*np.pi
    #v
    if np.dot(r,v)>0:
        true = np.arccos(np.dot(evec,r)/(length(evec)*length(r)))
    elif np.dot(r,v)<0:
        true = 2*np.pi - np.arccos(np.dot(evec,r)/(length(evec)*length(r)))
    elif np.dot(r,v)==0:
        true = np.arccos(np.dot(I,r)/length(r))
  #  print("h = " + str(h))
    print("a,e,i,O,w,v")

    return a , e , np.rad2deg(i) , np.rad2deg(O), np.rad2deg(w) , np.rad2deg(true) 


# In[76]:

def OE2cart(a,e,i,O,w,v):
    #first to put it in terms of radians...
    i = np.deg2rad(i)    
    O = np.deg2rad(O)
    w = np.deg2rad(w)            
    v = np.deg2rad(v)
    if a==0: 
        print("Error! a can't be zero.")
        return
    if e>=1 and a>0: 
        print("Error! For e>1, a must be negative. ")
        return 
    if e<1 and a<0: 
        print("Error! For 0<=e<1, a must be positive. ")
        return
    if e<0: 
        print("Error! e cannot be negative.")
        return 
    if i>np.pi or i<0: 
        print("Error! The inclination must be bewtween 0 and 180 degrees")
        return 
    if O<-np.pi or O>2*np.pi: 
        print("Error! The Longitude of the Ascending Node must be between -180 and 360 degrees. ")
        return
    if w<-np.pi or w>2*np.pi: 
        print("Eror! The Argument of Pericenter must be between -180 and 360 degrees. ")
        return
    
    if e<1 and v<-np.pi or v>2*np.pi:
        print("Error! For Elliptical Orbits, the True Anomaly must be between -180 and 360 degrees. ")
    if e>1:
       #define restricted range of here based on e... then make sure it is less than that! 
        if abs(v)>(180-np.rad2deg(np.arccos(1/e))):
            print("Error! For your chosen Eccentricity (e=" + str(e) + "), the True Anomaly is restricted to the range: " + str(-(180-np.rad2deg(np.arccos(1/e)))) + " < v < " + str((180-np.rad2deg(np.arccos(1/e)))) + ".")  
            return

    GM=1
    p=a*(1-e**2)
#start with work in a perifocal coordinate frame... 
    r=p/(1+e*np.cos(v))
    #magnitude of r right here...
    #now magnitude of v...
#v=sqrt(GM/p)[-sinvX+(e+cosv)Y]

    vperifocal=np.sqrt(GM*((2/r)-(1/a))) 
    
    rP=r*np.cos(v)
    rQ=r*np.sin(v)
    rperifocal=[[rP],[rQ],[0]]
    vP=np.sqrt(GM/p)*(-np.sin(v))
    vQ=np.sqrt(GM/p)*(e+np.cos(v))
    vperifocal = [[vP],[vQ],[0]]
    #obtained from Fundamental of Astrodynamics, using a matrix transformation from perifocal to Cartesion coord systems. 
    R11=(np.cos(O)*np.cos(w)-np.sin(O)*np.sin(w)*np.sin(i))
    R12 = -np.cos(O)*np.sin(w)-np.sin(O)*np.cos(w)*np.cos(i)
    R13 = np.sin(O)*np.sin(i)
    R21 = np.sin(O)*np.cos(w) + np.cos(O)*np.sin(w)*np.cos(i)
    R22 = -np.sin(O)*np.sin(w) + np.cos(O)*np.cos(w)*np.cos(i)
    R23 = -np.cos(O)*np.sin(i)
    R31 = np.sin(w)*np.sin(i)
    R32 = np.cos(w)*np.sin(i)
    R33 = np.cos(i)
    R=[[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]]
     
    rcart=[[0],[0],[0]]
    vcart=[[0],[0],[0]]
    for i in range(len(R)):
         for j in range(len(rperifocal[0])):
            for k in range(len(rperifocal)):
                rcart[i][j] += R[i][k] * rperifocal[k][j]
                vcart[i][j] +=R[i][k] *vperifocal[k][j]
    print("x,y,z,vx,vy,vz")
    return rcart[0][0],rcart[1][0],rcart[2][0],vcart[0][0],vcart[1][0],vcart[2][0]


# This is my testing ground... easy cut and paste of one set of variable into the other, and see if the desired result is recovered after going through both of the functions. 

