# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:31:49 2021

@author: drewj
"""

import electronclass as ec
import allfunctions as af
import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from random import uniform
from timeit import default_timer as timer
from tqdm import tqdm



#########################################
########    CONSTANTS           #########
#########################################

eps= 8.854*10**-12   #C^2/(N*m^2)
k= 1/(4*math.pi*eps)   #kg⋅m³⋅s⁻²⋅C^-2
mass= 9.109*10**-31     #kg
elecCharge=-1.602*10**-19  #C




#########################################
########    ODE Model           #########
#########################################
def model1(z,t):
        return[z[1],electron.getFor()/mass]
    
    
########################################
########   INITIAL SETTINGS    #########
########################################

#sets the timepoints
n=50000
t=np.linspace(0,20,n)

#number of electrons in system
num_e = 100

#list for initial positions
initpos=[]  

#lists of electrons
electrons,staticelectrons = [],[]

#creates random starting positions for each electrons b/w a certain range

for i in range (num_e):
    r= uniform(-50,50)
    if r not in initpos:
        initpos.append(r)
        
initpos=sorted(initpos)

for i in range (num_e):
    electrons.append(ec.electron(initpos[i],0,0,n))
    
    

# use specified initial positions for testing 
# electrons.append(ec.electron(-16.28813296526114,0,0,n))
# electrons.append(ec.electron(-11.782512519758592,0,0,n))
# electrons.append(ec.electron(-6.981923432847982,0,0,n))
# electrons.append(ec.electron(-2.0803642640444395,0,0,n))
# electrons.append(ec.electron(2.807485724143389,0,0,n))
# electrons.append(ec.electron(7.5977441171776965,0,0,n))
# electrons.append(ec.electron(12.203863752813945,0,0,n))
# electrons.append(ec.electron(16.482007581328062,0,0,n))

leftbound,rightbound=-7,7
staticelectrons.append(ec.electron(rightbound,0,0,n))   #right barrier
staticelectrons.append(ec.electron(leftbound,0,0,n))  #left barrier
    
    
#creates arrays for energy values
ke, pe, totale = np.zeros(n),np.zeros(n),np.zeros(n)


#starts timer to check how long the main loop takes



########################################
########      MAIN LOOP        #########
########################################

for x in tqdm(range (0,n)):
    start=timer()
    ke_now,pe_now,counter = 0,0,0
    pos_now = []
    
    for electron in electrons:                   #makes an array of current electron positions and calculates KE
        pos_now.append(electron.getPos())
        ke_now +=af.Kinetic(electron.getVel())  
    
    for electron in electrons:      #solves for PE
        currpos=electron.getPos()
        for i in range (counter+1,num_e):
            pe_now +=af.PotE(abs(af.getDistance(currpos, pos_now[i])))
        counter +=1  
        pe_now+=(af.PotE(abs(af.getDistance(currpos,rightbound)))
                  +af.PotE(abs(af.getDistance(currpos,leftbound))))
        fnet = 0
        #solves for total force on each electron
        for i in range (num_e):    
            distance = af.getDistance(currpos,pos_now[i])
            if (distance==0):
                fnet += 0
            elif (distance>0):
                fnet += k*(elecCharge**2)/(distance**2) #force will push electron right
            else:
                fnet -= k*(elecCharge**2)/(distance**2) #force will push electron left
                
        #pe_now+=(af.PotE(abs(af.getDistance(currpos,rightbound)))+af.PotE(abs(af.getDistance(currpos,leftbound))))
        #force from the 2 boundary electrons
        fnet+= (-k*(elecCharge**2)/(af.getDistance(currpos,rightbound)**2)
                +k*(elecCharge**2)/(af.getDistance(currpos,leftbound)**2))
        
        #attrition force
        if(x>0):
            fnet-=electron.getVel()*(500*10**-30)
        
        electron.ChangeFor(fnet)
        z = odeint(model1,electron.getPosVel(),t)   #gets the new pos and vel
        electron.ChangePosVel(z[1][0],z[1][1],x)    #updates pos and vel
    
    # for electron in staticelectrons:
    #     electron.ChangePos(electron.getPos(),x)

   # pe_now=pe_now/2
    ke[x],pe[x] = ke_now,pe_now         #adds the new KE and PE to arrays
    totale[x]=pe[x] + ke[x]             #addds total energy to array
    end=timer()
    #print(end-start)
    
for electron in electrons:
    print(electron.getPos())
    
avgdeltax=0
for i in range(num_e):
    if(i==0):
        avgdeltax+=abs(electrons[i].getPos()-leftbound)
    elif(i==num_e-1):
        avgdeltax+=abs(rightbound-electrons[i].getPos())
    else:
        avgdeltax+=abs(electrons[i+1].getPos()-electrons[i].getPos())

print(avgdeltax/num_e)


#########################################
########       PLOTTING         #########
#########################################


#position vs time
for electron in electrons:
    plt.plot(t,electron.getPosList())
    plt.xlabel('time')
    plt.ylabel('position')
# for electron in staticelectrons:
#     plt.plot(t,electron.getPosList(),'--')
#     plt.xlabel('time')
#     plt.ylabel('position')
# plt.axis([14,25,-55,55])
plt.show()

# plt.plot(t,electrons[0].getPosList())
# plt.xlabel('time')
# plt.ylabel('position')
# plt.show()

    
# #velocity vs time
for electron in electrons:
    plt.plot(t,electron.getVelList())
    plt.xlabel('time')
    plt.ylabel('velocity')
plt.show()


# #energy vs time
plt.plot(t,totale,'g-')
plt.plot(t,pe,'r-')
plt.plot(t,ke,'b-')
plt.xlabel('time')
plt.ylabel('energy')
plt.legend(['total','PE','KE'])
# #plt.axis([-0.5,60,0,(1.5)*(10**-28)])
plt.show()


#phase space
for electron in electrons:
    plt.plot(electron.getPosList(),electron.getVelList())
    plt.xlabel("position")
    plt.ylabel("velocity")
plt.show()
