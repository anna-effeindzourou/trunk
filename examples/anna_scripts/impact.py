# -*- coding: utf-8 -*-
from yade import ymport, plot
from math import *

import numpy
#import getpass
#user=getpass.getuser()
#sys.path.append('/home/'+user+'/YADE-git/install/bin')
#from Functions import *

try:
  os.mkdir('data')
except:
  pass
try:
  os.mkdir('paraview')
except:
  pass
      
#################################
### DEFAULT PARAMETERS        ###
#################################

####-------------------------------------
#### set False when running in batch mode
defaultTable = True
#defaultTable = False
####-------------------------------------
####-------------------------------------


utils.readParamsFromTable(
unloading = 8,
noTableOk = True
)

from yade.params.table import *

if defaultTable:
  O.tags['description'] = 'unloading_'+str(unloading)

print O.tags['description'], ' unloading =',unloading
#TAG

#### identify if script is running in batch mode
isBatch = runningInBatch()

###############################
###   LOAD SPHERE PACKING   ###
###############################
### load net, final position for gravity
O.load(O.tags['description']+'.xml.bz2')
O.tags['description'] = 'impact_'+O.tags['description']
gravelMat=O.materials['gravel'].id

######################
###   PARAMETERS   ###
######################
#BOX DIMENSION
a=1
h=1
d=0.5

#INITIAL ROCK POSITION
startxyz=Vector3(-0.75,0,0.8)

Kn=170e6 #normal stiffness Plassiard thesis
Ks=35e6 #tangential stiffness Plassiard thesis

#ENGINES PARAMETERS
g=9.81
damping_coeff = 0.0
dt=1e-5

#gravelPhys
krot=1.8
eta=1.8
alpha=1.0
c=20e3

#HIGH OF THE GRANULAR MATERIAL IN THE BOX
d_s=hMaxSpheres(d)
############################
####   CLUMP CREATION     ###
#############################
 
for i in O.bodies:
  if(type(i.shape)==Sphere):
    for j in range(2):
      i.state.angVel[j]=0
      i.state.vel[j]=0
O.step()
for i in O.interactions:  
  if((O.bodies[i.id1].isClumpMember==False) and (O.bodies[i.id2].isClumpMember==False) and (type(O.bodies[i.id1].shape)==Sphere) and (type(O.bodies[i.id2].shape)==Sphere)):
    O.bodies[i.id1].shape.color=O.bodies[i.id2].shape.color
    O.interactions.erase(i.id1,i.id2)
    O.bodies.clump([i.id1,i.id2])
############################
### DEFINING MATERIAL    ###
############################
n=1960
blocMat = O.materials.append(FrictMat(young=Kn ,poisson=Ks/Kn,frictionAngle=radians(30),density=44.5/(4./3.*pi*0.02**3*n),label='block'))

###############################
### LOADING OF THE BLOCK    ###
###############################

bloc=ymport.text('block.spheres',scale=1,shift=startxyz,mask=2,wire=False,highlight=False,color=[1,0,0],material=blocMat)
### rotate block

p0=startxyz
rot=Quaternion((0,0,1),radians(40))
for b in bloc: b.state.pos=p0+rot*(b.state.pos-p0)


clump,clumpspheres=O.bodies.appendClumped(bloc)

#############################
###   IMPACT PARAMETERS   ###
#############################

###test 1
#x=0
#y=-15
#z=-10
#w=2.56

#test 2
#x=0
#y=0
#z=-15
#w=0

###test 3
x=1.96
y=-1.72
z=-11.55
w=1.95*2*pi

###test 4
#x=0.21
#y=-0.7
#z=-14.5
#w=-2.56

O.bodies[clump].state.vel=Vector3(x,y,z) 
O.bodies[clump].state.angMom=Vector3(O.bodies[clump].state.inertia[0]*w,0.0,0)



### time step definition for simulation
O.dt = dt
############################
###   DEFINING ENGINES   ###
############################

O.engines=[
  ForceResetter(),
  InsertionSortCollider([Bo1_Sphere_Aabb(label='aabb'),Bo1_Facet_Aabb()]), 
  InteractionLoop(
  [Ig2_Sphere_Sphere_ScGeom(label='Ig2ssGeom'),Ig2_Facet_Sphere_ScGeom(label='Ig2fsGeom')],
  [Ip2_FrictMat_FrictMat_GravelPhys(alpha=MatchMaker(matches=((gravelMat,gravelMat,alpha),(blocMat,gravelMat,0.0))),krot=MatchMaker(matches=((gravelMat,gravelMat,krot),(blocMat,gravelMat,0.0))),eta=MatchMaker(matches=((gravelMat,gravelMat,krot),(blocMat,gravelMat,0.0))), psi_unloading=MatchMaker(matches=((gravelMat,gravelMat,unloading),(blocMat,gravelMat,1.0))), c=MatchMaker(matches=((gravelMat,gravelMat,c),(blocMat,gravelMat,0.0))),  label='block_particles')],
  [ Law2_ScGeom_GravelPhys_Plassiard(includeMoment=False, label='Law')]),
  NewtonIntegrator(damping=damping_coeff,gravity=(0,0,-g)),
  PyRunner(initRun=True,iterPeriod=100,command='addPlotData()'),
  VTKRecorder(iterPeriod=500,initRun=True,fileName='paraview/'+O.tags['description']+'-',recorders=['spheres','velocity']),
  ]

#################################
### PLOT + SIMULATION STORAGE ###
#################################

######################################
plot.plots={'t':('v',None,('z','r--'))}


def addPlotData():
	t=O.time
	plot.addData( t=O.time, x=O.bodies[clump].state.pos[0], y=O.bodies[clump].state.pos[1], z=O.bodies[clump].state.pos[2], vx=O.bodies[clump].state.vel[0], vy=O.bodies[clump].state.vel[1], vz=O.bodies[clump].state.vel[2], v=numpy.linalg.norm(O.bodies[clump].state.vel), wx=O.bodies[clump].state.angVel[0], wy=O.bodies[clump].state.angVel[1], wz=O.bodies[clump].state.angVel[2], w=numpy.linalg.norm(O.bodies[clump].state.angVel) )
plot.plot(noShow=False, subPlots=False)

def saveData():
	plot.saveDataTxt('data/'+O.tags['description']+'.dat',vars=('t','x','y','z','vx','vy','vz','v','wx','wy','wz','w'))

	

from yade import timing
O.timingEnabled=1
	
##########################
###  YADE-BATCH MODE   ###
##########################	

  
if not isBatch:
#  VISUALIZATION  
  from yade import qt
  qt.Controller()
  qtv = qt.View()
  qtr = qt.Renderer()
  plot.plot(noShow=False, subPlots=False)
else:
  O.run(100000,True)
  O.wait()
  saveData()