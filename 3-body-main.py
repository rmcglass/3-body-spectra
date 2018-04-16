import numpy as np
import matplotlib.pyplot as plt
##### initialize variables
maxTime = 15
G=1
# make masses
m_1 = 1
m_2 = 1
m_3 = 1
#initial velocity vectors
p1 = 0.51394
p2 = 0.30474
vi_1 = np.array([p1,p2])
vi_2 = np.array([p1,p2])
vi_3 = np.array([-2*p1,-2*p2])
#intitial positions
pos1 = np.array([-1.0,0.0])
pos2 = np.array([1.0,0.0])
pos3 = np.array([0.0,0.0])
#position arrays
posArray1 = []
posArray2=[]
posArray3=[]

def accel(posa,posb,mb):
    #returns the acceleration of particle a caused by its gravitational interaction with particle b
    normvec = (posb - posa)/ np.linalg.norm(posa-posb)
    accel = G*mb*normvec/(np.linalg.norm(posa-posb)**2.0)
    return accel

def computeEnergy(vel1,vel2,vel3,pos1,pos2,pos3):
    #we can use this function to make sure that our numerical solution conserves energy.
    kineticterm = 0.5*(m_1*np.linalg.norm(vel1)**2 + m_2*np.linalg.norm(vel2)**2 + m_3*np.linalg.norm(vel3)**2)
    potentialterm = -2*G*(m_1*m_2/np.linalg.norm(pos1-pos2)+m_1*m_3/np.linalg.norm(pos1-pos3)+m_2*m_3/np.linalg.norm(pos2-pos3))
    energy = potentialterm + kineticterm
    return energy

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
x1updated=pos1
x2updated=pos2
x3updated=pos3

energyvals=[]

for t in np.arange(0,maxTime,0.1):
    dt = 0.1
    if t == 0:
        velupdated1 = velHalfStep(v05_1,dt,pos1,pos2,pos3,m_2,m_3)
        velupdated2 = velHalfStep(v05_2,dt,pos2,pos3,pos1,m_3,m_1)
        velupdated3 = velHalfStep(v05_3,dt,pos3,pos1,pos2,m_1,m_2)
    else:
        velupdated1 = velHalfStep(velupdated1,dt,x1updated,x2updated,x3updated,m_2,m_3)
        velupdated2 = velHalfStep(velupdated2,dt,x2updated,x3updated,x1updated,m_3,m_1)
        velupdated3 = velHalfStep(velupdated3,dt,x3updated,x1updated,x2updated,m_1,m_2)
    updatedenergy = computeEnergy(velupdated1, velupdated2, velupdated3, x1updated, x2updated, x3updated)
    energyvals.append(updatedenergy)
    x1updated = posNew(x1updated,dt,velupdated1)
    x2updated = posNew(x2updated,dt,velupdated2)
    x3updated = posNew(x3updated,dt,velupdated3)
    posArray1.append(x1updated)
    posArray2.append(x2updated)
    posArray3.append(x3updated)


from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax=fig.add_subplot(111,projection='3d')
#----the 3d plot methods require that we use an array for x values, an array for y values, and an array for z values-
#this is the transpose of the arrays that we make in the loop above, so i just switch them around here.
arrayforplots1= np.transpose(posArray1)
arrayforplots2=np.transpose(posArray2)
arrayforplots3=np.transpose(posArray3)

plt.plot(arrayforplots1[0],arrayforplots1[1])
plt.plot(arrayforplots2[0],arrayforplots2[1])
plt.plot(arrayforplots3[0],arrayforplots3[1])

#for the time being, we'll try and look at the 2d motion of our particles. this should be fine, because if they
#have no initial z offset or velocity, we shouldnt have any motion or acceleration in that dimension

#ax.plot(arrayforplots1[0],arrayforplots1[1],arrayforplots1[2])
#ax.plot(arrayforplots2[0],arrayforplots2[1],arrayforplots2[2])
#ax.plot(arrayforplots3[0],arrayforplots3[1],arrayforplots3[2])

plt.show()

plt.plot(np.arange(0,maxTime,0.1),energyvals)
plt.show()