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


def accel(posa,posb,ma):
    #returns the acceleration of particle b caused by its gravitational interaction with particle a
    return -G*ma/np.abs(posa-posb)
def posNew(pos, dt, Vn1half):
    return pos + dt*Vn1half
def Vn3half():