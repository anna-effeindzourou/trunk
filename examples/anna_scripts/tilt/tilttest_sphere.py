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
young_box=30e7
poison_box=0.25
friction_box=27
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
		Bo1_Wall_Aabb(),
		Bo1_PFacet_Aabb(),
		Bo1_Facet_Aabb(),
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
		Law2_ScGridCoGeom_FrictPhys_CundallStrack(),
		]
	),
]

############################
### DEFINING MATERIAL    ###
############################
O.materials.append(CohFrictMat(young=0.01*young_box,poisson=poison_box,density=density,frictionAngle=friction_box,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=True,alphaKr=5,label='NodeMat'))

boxMat = O.materials.append(FrictMat(young=0.01*young_box ,poisson=poison_box,frictionAngle=friction_box,density=density))
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
### DEFINING OBJECTS     ###
############################

#O.bodies.append(geom.facetBox(center=(0.,0.,0),extents=(a,a,h),orientation=Quaternion((0.,0,0),0),wallMask=31,wire=True,dynamic=False,material=boxMat,color=(0,1,1)))
# Create spheres
#Using dense packing
#pred=pack.inAlignedBox((-a,-h,-d),(a,h,d))
#spheres=pack.randomDensePack(pred,spheresInCell=nSpheres,radius=Rs,rRelFuzz=Rf)
#O.bodies.append(spheres)


O.bodies.append(sphere(center=(-x+5*0.05,-y,h+0.05),radius=.01,fixed=False))

O.bodies.append(sphere(center=(-x,-y+5*0.05,h+0.05),radius=.01,fixed=False))

O.bodies.append(sphere(center=(-x+5*0.05,2*x,h+0.05),radius=.01,fixed=False))



O.bodies.append(sphere(center=(x,2*x-0.05,h+0.05),radius=.01,fixed=False))


O.bodies.append(sphere(center=(-x+2*0.05,-x+6*0.05,h+0.05),radius=.01,fixed=False))


############################
###   DEFINING ENGINES   ###
############################
angularVelocity=0.1
idsToRotate=nodesIds
#idsToRotate.append(bottomBox)
O.dt=0.00001
O.engines+=[

	NewtonIntegrator(gravity=(0,0,-9.81),damping=0.2),
	#NewtonIntegrator(gravity=[0,0,-g],damping=0.,label='Newton_Integrator'),
]
def rotate():
	O.engines+=[RotationEngine(ids=idsToRotate,angularVelocity=angularVelocity,rotateAroundZero=True,zeroPoint=(-rWall,0,0),rotationAxis=(1,0,0),label='rotationEngine'),]
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
