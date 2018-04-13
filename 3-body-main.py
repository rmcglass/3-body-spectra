import numpy as np
##### initialize variables
maxTime = 50
G=1
# make masses
m_1 = 10
m_2 = 10
m_3 = 10
#initial velocity vectors
vi_1 = [50,0,0]
vi_2 = [0,50,0]
vi_3 = [0,0,50]
#intitial positions
posi_1 = [100,0,0]
posi_2 = [0,100,0]
posi_3 = [0,0,100]
#initial positions
def accel(posa,posb,mb):
    #returns the acceleration of particle a caused by its gravitational interaction with particle b
    return -G*mb/np.abs(posa-posb)

def velHalfStep(vprev,timestep,posa,posb,posc,mb,mc):
    #returns the updated velocity of particle a based on the gravitational interaction with particles b and c.
    accelbona=accel(posa,posb,mb)
    accelcona=accel(posa,posc,mc)
    totalaccel=accelbona+accelcona
    return vprev+totalaccel*timestep

def posNew(pos, dt, Vn1half):
    return pos + dt*Vn1half
def Vn3half():