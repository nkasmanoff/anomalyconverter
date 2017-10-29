
# coding: utf-8

# In[ ]:

import matplotlib.pyplot as plt
import math
import numpy as np

#including the case of e->1, and letting mean anomaly vary from -180 to 180
#like its done in the online sources 

def E2T(E,e): 
    if E>360 or E<0:
        print("Has to be from 0 to 360 degrees!")
    else:
        if e<0: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!")
        elif e==0:
            print("This is a circular orbit")
            E=np.deg2rad(E)
            v=E #true for circular orbits...
        elif e>1: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!") 
        elif e==1:
            e=0.999999 #approximating 1 to be very close to 1. 
            E=np.deg2rad(E)
            v=2*math.atan(np.sqrt((1+e)/(1-e))*math.tan(E/2)) 
        else:
            E=np.deg2rad(E)
            v=2*math.atan(np.sqrt((1+e)/(1-e))*math.tan(E/2))
        return np.rad2deg(v)

def T2E(v,e):
#given v and e, convert to E
    if v>180 or v<-180:
        print("Has to be from -180 to 180 degrees!")
    else:
        if e<0: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!")
        elif e==0:
            print("This is a circular orbit")
            v=np.deg2rad(v)
            E=v
            if E<0:
                E=E+2*np.pi
            return np.rad2deg(E)
        elif e>1: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!") 
        elif e==1: 
            e=0.999999
            v=np.deg2rad(v)
            E=2*math.atan(np.sqrt((1-e)/(1+e))*math.tan(v/2))
            if E<0:
                E=E+2*np.pi
        else:
            v=np.deg2rad(v)
            E=2*math.atan(np.sqrt((1-e)/(1+e))*math.tan(v/2))
            if E<0:
                E=E+2*np.pi
        return np.rad2deg(E)

def E2M(E,e):
#given E and e, convert to mean anomaly M. 
    if E>360 or E<0:
        print("Has to be from 0 to 360 degrees!")
    else:
        if e<0: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!")
        elif e==0:
            print("This is a circular orbit.")
            E=np.deg2rad(E)
            v=E
        elif e>1: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!") 
        elif e==1: 
            e=0.999999
            E=np.deg2rad(E)
            M=E-e*np.sin(E)
        else:
            E=np.deg2rad(E)
            M=E-e*np.sin(E)
        return np.rad2deg(M)

def M2E(M,e):
#given M and e, convert to E
    if M>360 or M<0:
        print("Has to be from 0 to 360 degrees!")
    else:
        if e<0: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!")
        elif e==0:
            print("This is a circular orbit.")
            M=np.deg2rad(M)
            E=M #true for circular orbits...
            if E<0:
                E=E+2*np.pi
            return np.rad2deg(v)
        elif e>1: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!") 
        elif e==1:
            e=0.999999
            M=np.deg2rad(M)
    #other conditions about range of M and e, will come in soon 
    #now to apply newton Method
            Eo=M
            for i in range(0,1000):
                E1=Eo-(Eo-e*np.sin(Eo)-M)/(1-e*np.cos(Eo))
                Eo=E1
                
        else:   
            M=np.deg2rad(M)
    #other conditions about range of M and e, will come in soon 
    #now to apply newton Method
            Eo=M
            for i in range(0,1000):
                E1=Eo-(Eo-e*np.sin(Eo)-M)/(1-e*np.cos(Eo))
                Eo=E1
                
        return np.rad2deg(Eo)

def T2M(v,e):
    if v>180 or v<-180:
        print("Has to be from -180 to 180 degrees!")
    else:
        if e<0: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!")
        elif e==0:
            print("This is a circular orbit.")
            v=np.deg2rad(v)
            E=v
            if E<0:
                E=E+2*np.pi
            return np.rad2deg(E)
        elif e>1: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!") 
        elif e==1:
            e=0.999999
            E=T2E(v,e)

        else:
            E=T2E(v,e)
        return E2M(E,e)

def M2T(M,e):
    if M>360 or M<0:
        print("Has to be from 0 to 360 degrees!")
    else:
        if e<0: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!")
        elif e==0:
            print("cool, its a circle!")
            M=np.deg2rad(M)
            E=M #true for circular orbits...
            if E<0:
                E=E+2*np.pi
            return np.rad2deg(v)
        elif e>1: 
            print("Sorry, the eccentricity of an orbit is between 0 and 1!") 
        elif e==1:
            e=0.99999
            E=M2E(M,e)
        else:  
            E=M2E(M,e)
        return E2T(E,e)

