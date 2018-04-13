import numpy as np
import matplotlib.pyplot as plt
##### initialize variables
maxTime = 50
G=1
# make masses
m_1 = 10
m_2 = 10
m_3 = 10
#initial velocity vectors
vi_1 = np.array([50,0,0])
vi_2 = np.array([0,50,0])
vi_3 = np.array([0,0,50])
#intitial positions
pos1 = np.array([10,0,0])
pos2 = np.array([0,10,0])
pos3 = np.array([0,0,10])
#initial positions

def accel(posa,posb,mb):
    #returns the acceleration of particle a caused by its gravitational interaction with particle b
    return G*mb*(posa-posb)/np.linalg.norm(posa-posb)**3

#initial half-step velocity
v05_1 = vi_1+0.05*(accel(pos1,pos2,m_2)+accel(pos1,pos3,m_3))
v05_2 = vi_2+0.05*(accel(pos2,pos1,m_2)+accel(pos2,pos3,m_3))
v05_3 = vi_3+0.05*(accel(pos3,pos1,m_2)+accel(pos3,pos2,m_3))

def velHalfStep(vprev,timestep,posa,posb,posc,mb,mc):
    #returns the updated velocity of particle a based on the gravitational interaction with particles b and c.
    accelbona=accel(posa,posb,mb)
    accelcona=accel(posa,posc,mc)
    totalaccel=accelbona+accelcona
    velhalf = vprev+totalaccel*timestep
    return velhalf

def posNew(pos,dt,velupdated):
    return pos + dt*velupdated

velupdated1 = 0.0
velupdated2 = 0.0
velupdated3 = 0.0


for t in np.arange(0,100,0.1):
    dt = 0.1
    if t == 0:
        velupdated1 = velHalfStep(v05_1,dt,pos1,pos2,pos3,m_2,m_3)
        velupdated2 = velHalfStep(v05_2,dt,pos2,pos3,pos1,m_3,m_1)
        velupdated3 = velHalfStep(v05_3,dt,pos3,pos1,pos2,m_1,m_2)
    else:
        velupdated1 = velHalfStep(velupdated1,dt,pos1,pos2,pos3,m_2,m_3)
        velupdated2 = velHalfStep(velupdated2,dt,pos2,pos3,pos1,m_3,m_1)
        velupdated3 = velHalfStep(velupdated3,dt,pos3,pos1,pos2,m_1,m_2)
    x1updated = posNew(pos1,dt,velupdated1)
    x2updated = posNew(pos2,dt,velupdated2)
    x3updated = posNew(pos3,dt,velupdated3)
    print x1updated

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax=fig.add_subplot(111,projection='3d')