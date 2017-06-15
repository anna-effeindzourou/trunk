# -*- coding: utf-8
from yade import ymport, utils,pack,export,qt
import gts,os
from yade import geom
#import matplotlib
from yade import plot
#from pylab import *
#import os.path, locale


def hMax(n):
	idHMax=0
	hMax=-1000000.0
	for i in O.bodies:
		h=i.state.pos[n]
		if (h>hMax) and (type(i.shape)==Sphere):
			hMax=h
			idHMax=i.id
	return (hMax)	
    
    
def hMin(n):
    idHMin=0
    hMin=100000.0
    for i in O.bodies:
	h=i.state.pos[n]
	if (h<hMin):
	   hMin=h
	   idHMin=i.id
    return (hMin)      
#Function in order to calculate rmin (minimum radius) and rmax (maximum radius)
def MinMax():
    rmax=0
    rmin=10
    r=0
    for i in O.bodies:
      if(type(i.shape)==Sphere):
	r=i.shape.radius
	if(r>rmax):
	  rmax=r
	if(r<rmin):
	  rmin=r
    l=[rmin,rmax]
    return (l)


def sup(d):
	for i in O.bodies:
		if (type(i.shape)==Sphere) and (i.state.pos[2]>d):
			O.bodies.erase(i.id)    
  
def writeFile():
	yade.export.text('spheres_1e-02.txt')




####################
###  MATERIAL    ###
####################
poisson=0.28
E=2*7.9e10*(1+poisson) ##1e11
density=7.8e8
frictionAngle=0.096
frictionAngleW=0.228

O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngleW,label='Wallmat'))
O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngle,label='Smat'))



##########################
###  SPHERE PACKING    ###
##########################
#### Parameters the cylinder ###
L=0.203 #high [m]
l=0.0505	# radius [m]
color=[155./255.,155./255.,100./255.]
radius=1e-02


kwBoxes={'color':[1,0,0],'wire':True,'dynamic':False,'material':1}
O.bodies.append(utils.geom.facetCylinder(center=Vector3(0,0,L/2.), radius=l, height=L, orientation=Quaternion((1, 0, 0), 0),**kwBoxes))

###erase the top and bottom facet of the cylinder
for i in range(0,40,4):
	O.bodies.erase(i)
for i in range(1,38,4):
	O.bodies.erase(i)


predicate=inCylinder(centerBottom=Vector3(0,0,0), centerTop=Vector3(0,0,L+L/2.), radius=l-0.005)	
sp=SpherePack()
sp=pack.randomDensePack(predicate, radius=radius, material='Smat',  cropLayers=10, rRelFuzz=0.0, spheresInCell=100,returnSpherePack=True)
sp.toSimulation()


########################
#### WALL GENERATION  ##
########################
O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngleW,label='Wmat'))
topPlate=utils.wall(position=hMax(2)+radius*10,sense=0, axis=2,color=Vector3(1,0,0),material='Wmat')
O.bodies.append(topPlate)
bottomPlate=utils.wall(position=0,sense=0, axis=2,color=Vector3(1,0,0),material='Wmat')
O.bodies.append(bottomPlate)

######################
#### MOVE TOP WALL  ##
######################
v=1.7e-03
def movewall(v):
	topPlate.state.pos=Vector3(0,0,hMax(2)+radius)
	topPlate.state.vel=Vector3(0,0,-v)
	

def dataCollector():
	S=pi*l**2
	Fnt=O.forces.f(topPlate.id)[2]
	Fnb=O.forces.f(bottomPlate.id)[2]
	#sigma=Fnb/S
	plot.addData(t1=O.time,t2=O.time,Fnb=Fnb,Fnt=Fnt)
    
   
plot.plots={'t1':('Fnb'),'t2':('Fnt')}
plot.plot(noShow=False, subPlots=True)

#########################
### ENGINE DEFINITION  ##
#########################  
O.dt=0.5*PWaveTimeStep()
O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_Sphere_Aabb(),
		Bo1_Wall_Aabb(),
		Bo1_Facet_Aabb(),
	]),
	InteractionLoop([
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_Facet_Sphere_ScGeom(),
		Ig2_Wall_Sphere_ScGeom()
		],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True),
		Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		Law2_ScGeom_FrictPhys_CundallStrack(),
		Law2_ScGridCoGeom_FrictPhys_CundallStrack(),
		Law2_GridCoGridCoGeom_FrictPhys_CundallStrack()
		]
	),
	DomainLimiter(lo=(-l,-l,0),hi=(l,l,1),iterPeriod=200),
	NewtonIntegrator(damping=0.7,gravity=(0,0,-9.81),label='Newton'),
	PyRunner(initRun=True,iterPeriod=1,command='dataCollector()'),
	]	