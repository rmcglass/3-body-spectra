import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
##### initialize variables
maxTime = 100000
G=1.0
dt = 10
# make asses
m_1 = 1
m_2 = 1
m_3 = 1
#initial velocity vectors
p1 = 0.51394
p2 = 0.30474
p3 = 0
vi_1 = np.array([p1,p2,p3])
vi_2 = np.array([p1,p2,p3])
vi_3 = np.array([-2*p1,-2*p2,-2*p3])
#intitial positions
projPosArray1=[]
projPosArray2=[]
projPosArray3=[]

pos1 = np.array([-1.0,0.0,0.0])
pos2 = np.array([1.0,0.0,0.0])
pos3 = np.array([0.0,0.0,0.0])
# position arrays
posArray1 = []
posArray2 = []
posArray3 = []
#line of sight velocity arrays
losVelArray1=[]
losVelArray2=[]
losVelArray3=[]

velup1 = 0.0
velup2 = 0.0
velup3 = 0.0
x1updated=pos1
x2updated=pos2
x3updated=pos3

def accel(posa,posb,mb):
    #returns the acceleration of particle a caused by its gravitational interaction with particle b
    vecmag = np.linalg.norm(posa-posb)
    normvec = (posb - posa)/ vecmag
    if vecmag < .05:
        accel = 0
    else:
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
v05_2 = vi_2+0.05*(accel(pos2,pos1,m_1,)+accel(pos2,pos3,m_3))
v05_3 = vi_3+0.05*(accel(pos3,pos1,m_1)+accel(pos3,pos2,m_2))

def velHalfStep(vprev,timestep,posa,posb,posc,mb,mc):
    #returns the updated velocity of particle a based on the gravitational interaction with particles b and c.
    accelbona=accel(posa,posb,mb)
    accelcona=accel(posa,posc,mc)

    totalaccel=accelbona+accelcona
    velhalf = vprev+totalaccel*timestep
    return velhalf

def velTwoBody(vprev,timestep,posa,posb,mb):
    #returns the updated velocity of a particle based on interaction with a SINGLE body
    accelbona=accel(posa,posb,mb)
    velhalf=vprev+accelbona*timestep
    return velhalf

def posNew(pos,dt,velupdated):
    return pos + dt*velupdated

def projectToLOS(vel,los):
    return np.dot(vel,los)

# projects 2D coords onto 3D sphere w/ fixed radius
def projPos():
    for i in range(0,len(posArray1)-1):
        oldX = posArray1[i][0]
        oldY = posArray1[i][1]
        newPos1 = np.array([2*oldX/(1+oldX**2+oldY**2),2*oldY/(1+oldX**2+oldY**2),(-1+oldX**2+oldY**2)/(1+oldX**2+oldY**2)])
        projPosArray1.append(newPos1)
    for i in range(0,len(posArray2)-1):
        oldX = posArray2[i][0]
        oldY = posArray2[i][1]
        newPos2 = np.array([2 * oldX / (1 + oldX ** 2 + oldY ** 2), 2 * oldY / (1 + oldX ** 2 + oldY ** 2),(-1 + oldX ** 2 + oldY ** 2) / (1 + oldX ** 2 + oldY ** 2)])
        projPosArray2.append(newPos2)
    for i in range(0,len(posArray3)-1):
        oldX = posArray3[i][0]
        oldY = posArray3[i][1]
        newPos3 = np.array([2 * oldX / (1 + oldX ** 2 + oldY ** 2), 2 * oldY / (1 + oldX ** 2 + oldY ** 2),
                            (-1 + oldX ** 2 + oldY ** 2) / (1 + oldX ** 2 + oldY ** 2)])
        projPosArray3.append(newPos3)

velupdated1 = 0.0
velupdated2 = 0.0
velupdated3 = 0.0
x1updated=pos1
x2updated=pos2
x3updated=pos3

energyvals=[]
def leapfrog(x1updated,x2updated,x3updated,losvec):
    # i wanted this as a function, not as a little few lines of code so i can call it a whole bunch of times
    for t in np.arange(0,maxTime,dt):
        if t == 0:
            velupdated1 = velHalfStep(v05_1,dt,x1updated,x2updated,x3updated,m_2,m_3)
            velupdated2 = velHalfStep(v05_2,dt,x2updated,x3updated,x1updated,m_3,m_1)
            velupdated3 = velHalfStep(v05_3,dt,x3updated,x1updated,x2updated,m_1,m_2)
        else:
            velupdated1 = velHalfStep(velupdated1,dt,x1updated,x2updated,x3updated,m_2,m_3)
            velupdated2 = velHalfStep(velupdated2,dt,x2updated,x3updated,x1updated,m_3,m_1)
            velupdated3 = velHalfStep(velupdated3,dt,x3updated,x1updated,x2updated,m_1,m_2)
        updatedenergy = computeEnergy(velupdated1, velupdated2, velupdated3, x1updated, x2updated, x3updated)
        energyvals.append(updatedenergy)
        x1updated = posNew(x1updated,dt,velupdated1)
        x2updated = posNew(x2updated,dt,velupdated2)
        x3updated = posNew(x3updated,dt,velupdated3)

        # uh.... okay. the animation isn't working for some reason, which is pretty sad i guess. I'll get the line of sight stuff
        # to do stuff rn instead.
        #i think the best way to do this will be to use the velocity at each position and just do a dot product
        #we can save the dotted velocities as an array as well

        losvel1= projectToLOS(velupdated1,losvec)
        losvel2= projectToLOS(velupdated2,losvec)
        losvel3= projectToLOS(velupdated3,losvec)

        losVelArray1.append(losvel1)
        losVelArray2.append(losvel2)
        losVelArray3.append(losvel3)

        posArray1.append(x1updated)
        posArray2.append(x2updated)
        posArray3.append(x3updated)

def leapfrogtwobody(x1updated,x2updated):
    for t in np.arange(0,maxTime,dt):
        if t == 0:
            velupdated1 = velTwoBody(v05_1,dt,x1updated,x2updated,m_2)
            velupdated2 = velTwoBody(v05_2,dt,x2updated,x1updated,m_1)
        else:
            velupdated1 = velTwoBody(velupdated1,dt,x1updated,x2updated,m_2)
            velupdated2 = velTwoBody(velupdated2,dt,x2updated,x1updated,m_1)
        #updatedenergy = computeEnergy(velupdated1, velupdated2, velupdated3, x1updated, x2updated, x3updated)
        #energyvals.append(updatedenergy)
        x1updated = posNew(x1updated,dt,velupdated1)
        x2updated = posNew(x2updated,dt,velupdated2)
        posArray1.append(x1updated)
        posArray2.append(x2updated)

#leapfrog(x1updated,x2updated,x3updated)
#--------------------RESTRICTED 3-BODY PROBLEM----------------------
#now that we have leapfrog defined, we can just change the values of our variables between runs of our simulation.
#we'll want to think about the restricted 3-body problem. in this case, we have two masses that are very, very massive - they will have stable orbits.
#first, i'm going to test this - well try to do a two body version of our code!
m_1=1000.0
m_2=1.0
m_3=10**-6
#m_3 will have no initial velocity
#we also have to reinitialize our position arrays and our initial velocities
posArray1=[]
posArray2=[]
posArray3=[]

losVelArray1=[]
losVelArray2=[]
losVelArray3=[]

vi_1 = np.array([0,-0.001,0])
vi_2 = np.array([0,0.999,-0.001])
vi_3 = np.array([0,0.0,0.0])

v05_1 = vi_1+0.05*(accel(pos1,pos2,m_2))
v05_2 = vi_2+0.05*(accel(pos2,pos1,m_2))

pos1 = np.array([0.0,0.0,0.0])
pos2 = np.array([1000.0,0.0,0.0])
pos3 = np.array([500.0,1.0,0.0])

los=[1,0,0]
leapfrog(pos1,pos2,pos3,los)
# -----------------------------------------------------------------
# -------------- PLOTS AND ANIMATIONS BELOW------------------------
# -----------------------------------------------------------------

plt.plot(np.arange(0,maxTime,dt),losVelArray1,label='massive obj')
plt.plot(np.arange(0,maxTime,dt),losVelArray2,label='less massive obj')
plt.plot(np.arange(0,maxTime,dt),losVelArray3,label='tiny obj')
plt.legend()

plt.show()
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax=fig.add_subplot(111, projection='3d')
line1, = ax.plot([], [], [], 'o-', lw=2)
line2, = ax.plot([], [], [], 'o-', lw=2)
line3, = ax.plot([], [], [], 'o-', lw=2)
#line1Data=np.empty()

#----the 3d plot methods require that we use an array for x values, an array for y values, and an array for z values-
#this is the transpose of the arrays that we make in the loop above, so i just switch them around here.

arrayforplots1= np.transpose(posArray1)
arrayforplots2=np.transpose(posArray2)
arrayforplots3=np.transpose(posArray3)

plt.plot(arrayforplots1[0],arrayforplots1[1], arrayforplots1[2], label='sun size mass') #,arrayforplots1[2])
plt.plot(arrayforplots2[0],arrayforplots2[1], arrayforplots2[2], label='jup size mass') #,arrayforplots2[2])
plt.plot(arrayforplots3[0],arrayforplots3[1], arrayforplots3[2], label='smaller size mass')
#plt.legend()
## make arrays suitable for projection
projPos()
arrayforprojplots1 = np.transpose(projPosArray1)
arrayforprojplots2 = np.transpose(projPosArray2)
arrayforprojplots3 = np.transpose(projPosArray3)

def plotInit():
    '''set up animation with initial conditions'''
    line1.set_data(arrayforplots1[0][0],arrayforplots1[1][0])
    line1.set_3d_properties(arrayforplots1[2][0])
    line2.set_data(arrayforplots2[0][0],arrayforplots2[1][0])
    line2.set_3d_properties(arrayforplots2[2][0])
    line3.set_data(arrayforplots3[0][0],arrayforplots3[1][0])
    line3.set_3d_properties(arrayforplots3[2][0])
    lines = [line1,line2,line3]
    return lines

def animate(i):
    '''step animation'''
    lines = [line1,line2,line3]
    arrays = [arrayforplots1,arrayforplots2,arrayforplots3]
    for line, array in zip(lines,arrays):
        line.set_data(array[0][:i],array[1][:i])
        line.set_3d_properties(array[2][:i])
    fig.canvas.draw()
    return lines

ani = animation.FuncAnimation(fig, animate, frames=len(arrayforplots3[0]), blit=True, init_func=plotInit)

#
# #for the time being, we'll try and look at the 2d motion of our particles. this should be fine, because if they
# #have no initial z offset or velocity, we shouldnt have any motion or acceleration in that dimension
#
plt.xlim(-1500.0,1500.0)
plt.ylim(-1500.0,1500.0)
plt.show()

#plt.plot(np.arange(0,maxTime,dt),energyvals)
#plt.show()

figp = plt.figure()
axp=figp.add_subplot(111, projection='3d')
axp.plot(arrayforprojplots1[0],arrayforprojplots1[1],arrayforprojplots1[2])
axp.plot(arrayforprojplots2[0],arrayforprojplots2[1],arrayforprojplots2[2])
axp.plot(arrayforprojplots3[0],arrayforprojplots3[1],arrayforprojplots3[2])
plt.show()

#-------------------------------------------------------------------------------------------------------------
#--------------------------------------SPECTRAL LINE BS-------------------------------------------------------

#first - define a method that returns a gaussian at the velocity of each of the parts of the system.
#for now: this line will have a size proportional to the mass of the object.
#(we'll plot this in log space, i guess)

def gaussianprofile(velpos,centvel,maxflux,width):
    #creates a gaussianprofile for each
    #Set the normalisation
    norm = maxflux/(width*np.sqrt(2.*np.pi))
    return norm*np.exp(-((velpos-centvel)**2.) / (2.*width**2.))

def createSumGauss(velpos, cvels,fluxes,widths):
    #this method creates the triple gaussian profile of all three lines
    centvel1,centvel2,centvel3 = cvels
    flux1,flux2,flux3 = fluxes
    width1,width2,width3 = widths
    retvals=[gaussianprofile(v,centvel1,flux1,width1)+gaussianprofile(v,centvel2,flux2,width2)+gaussianprofile(v,centvel3,flux3,width3) for v in velpos]
    return retvals
# plt.plot(np.arange(0,maxTime,0.1),energyvals)
# plt.show()