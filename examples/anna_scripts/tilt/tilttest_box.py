 # -*- coding: utf-8
from yade import ymport,utils,pack,export,qt,bodiesHandling
import gts,os

# for plotting
from math import *
from yade import plot


############################
### DEFINING PARAMETERS  ###
############################

#GEOMETRIC :dimension of the rectangular box
a=.2	 # side dimension 
h=.1	 # height 


#MATERIAL PROPERTIES
frictionAngle=radians(35)
density=2700
Kn=25e7 #normal stiffness Plassiard thesis
Ks=5e7 #tangential stiffness Plassiard thesis


#PARTICLE SIZE
Rs=0.015# mean particle radius
Rf=0.5 # dispersion (Rs±Rf*Rs)
nSpheres=2000# number of particles

#ENGINES PARAMETERS
g=9.81
damping_coeff = 0.5

#MATERIAL'S PARAMETERS
#box
young_box=30e5
poison_box=0.25
friction_box=radians(10)
density_box=1576.

#gravel
G=40e5
poisson_g=0.25
young_g=30e5
m=1.878
v=(4*pi*pow(Rs,3))/3

#TIME STEP
dt=1.e-7
#TAG

####################
###  ENGINES     ###
####################
O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_Sphere_Aabb(),
		Bo1_GridConnection_Aabb(),
		Bo1_PFacet_Aabb()
	]),
	InteractionLoop([
		Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		Ig2_GridConnection_PFacet_ScGeom(),
		Ig2_Sphere_GridConnection_ScGridCoGeom(),
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_PFacet_PFacet_ScGeom()
		],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=False),
		Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		Law2_ScGeom_FrictPhys_CundallStrack(),
		Law2_ScGridCoGeom_FrictPhys_CundallStrack(),Law2_GridCoGridCoGeom_FrictPhys_CundallStrack()
		]
	),
]



############################
### DEFINING MATERIAL    ###
############################
O.materials.append(CohFrictMat(young=young_box,poisson=poison_box,density=density,frictionAngle=friction_box,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=True,alphaKr=5,label='NodeMat'))

boxMat = O.materials.append(FrictMat(young=young_box ,poisson=poison_box,frictionAngle=friction_box,density=density))
sphereMat = O.materials.append(FrictMat(young=young_box ,poisson=poison_box,frictionAngle=friction_box,density=density))
rWall=0.05

x=a/2.+rWall
y=0
z=0
fac=2.+rWall
clump=0
color=Vector3(1,0,0)
fixed=False
rWall=0.008

#### define box

nodesIds=[]

nodesIds.append( O.bodies.append(gridNode([-x,-y,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,-y,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,-y,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,-y,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,2*x,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,2*x,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,2*x,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,2*x,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

cylIds=[]
pfIds=[]

pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[1],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[2],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))

#pfIds.append(pfacetCreator3(nodesIds[1],nodesIds[5],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
#pfIds.append(pfacetCreator3(nodesIds[5],nodesIds[6],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[5],nodesIds[6],nodesIds[1],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[6],nodesIds[2],nodesIds[1],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))

pfIds.append(pfacetCreator3(nodesIds[7],nodesIds[6],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[6],nodesIds[5],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))

#pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[4],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
#pfIds.append(pfacetCreator3(nodesIds[4],nodesIds[7],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[7],nodesIds[0],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[7],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=True,mask=-1))


nbodies=[]

#O.step()
for i in nodesIds:
		O.bodies[i].state.blockedDOFs='x'
		#nbodies.append(O.bodies[i])

#clump,clumpIds=O.bodies.appendClumped(nbodies)

#O.bodies[clump].state.blockedDOFs='xyzXYZ'

fixed=True
nodesIds=[]
nodesIds.append( O.bodies.append(gridNode([-x,-y,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,-y,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,-y,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,-y,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,2*x,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,2*x,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,2*x,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,2*x,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

cylIds=[]
pfIds=[]

pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[1],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[2],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))

pfIds.append(pfacetCreator3(nodesIds[1],nodesIds[5],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[5],nodesIds[6],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))

pfIds.append(pfacetCreator3(nodesIds[7],nodesIds[6],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[6],nodesIds[5],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))

pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[4],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[4],nodesIds[7],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))

pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[1],nodesIds[5],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))
pfIds.append(pfacetCreator3(nodesIds[4],nodesIds[0],nodesIds[5],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=None,color=color,fixed=fixed,mask=-1))



############################
###   DEFINING ENGINES   ###
############################

init=0
para=1000
coll=100
def rotate():
	global init
	init=O.iter
	Newton_Integrator.damping=0.2
	#for i in nodesIds1:
		#O.bodies[i].state.blockedDOFs='x'
	#O.bodies[clump].state.blockedDOFs='x'	
	O.engines+=[
		RotationEngine(ids=idsToRotate,angularVelocity=angularVelocity,rotateAroundZero=True,zeroPoint=(-rWall,0,0),rotationAxis=(1,0,0),label='rotationEngine'),
		#VTKRecorder(iterPeriod=para,dead=True,initRun=True,fileName='paraview/'+O.tags['description']+'_',recorders=['spheres','velocity','intr','materialId'],label='parav'),
		#PyRunner(initRun=True,iterPeriod=coll,command='dataCollector()')
	]




angularVelocity=0.1
idsToRotate=nodesIds

O.dt=0.0001
O.engines+=[NewtonIntegrator(gravity=(0,0,-9.81),damping=0.2,label='Newton_Integrator')]
############################
###   VISUALIZATION      ###
############################
renderer = qt.Renderer()
qt.View()
O.saveTmp()
#O.run()


def getAngle():
	alpha=angularVelocity*O.iter*O.dt/pi*180.
	print 'alpha = ', alpha
