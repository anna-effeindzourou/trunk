# -*- coding: utf-8
from yade import ymport, utils,pack,export,qt
import gts,os
from yade import geom
#import matplotlib
from yade import plot
#from pylab import *
#import os.path, locale
  
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
		Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
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
]

Ey=12.5e8
poisson=1
density=2.60e4
Kr=6250
frictionAngle=radians(00)
O.materials.append(CohFrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=frictionAngle,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=True,alphaKr=0,label='mat'))
matStif=O.materials.append(CohFrictMat(young=2*Ey,poisson=poisson,density=density,frictionAngle=frictionAngle,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=True,alphaKr=0,label='matStiff'))


O.materials.append(FrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=frictionAngle,label='Pmat'))

oriBody = Quaternion(Vector3(0,0,1),(pi/2))
nodesIds=[]
nodesIds1=[]
cylIds=[]
pfIds=[]
radius=.001

croisillons = 1


#### Parameters of a rectangular grid ###
L=0.203 #length [m]
l=0.0505	#half width	(radius) [m]
nbL=20 #number of nodes for the length	[#] doit etre paire
nbl=16	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!
#nbl=16
r=radius
color=[155./255.,155./255.,100./255.]
nodesIds=[]


predicate=inCylinder(centerBottom=Vector3(0,0,0), centerTop=Vector3(0,0,L+L/3.), radius=l-0.005)	
kwBoxes={'color':[1,0,0],'wire':True,'dynamic':False,'material':2}
O.bodies.append(utils.geom.facetCylinder(center=Vector3(0,0,+L/2.), radius=l, height=L, orientation=Quaternion((1, 0, 0), 0),**kwBoxes))
sp=SpherePack()
sp=pack.randomDensePack(predicate, radius=2.5e-03, material='Pmat',  cropLayers=10, rRelFuzz=0.0, spheresInCell=100,returnSpherePack=True)
sp.toSimulation()
print 'spheres=', len(O.bodies)

O.dt=1.140175425099138e-05

for i in range(0,40,4):
	O.bodies.erase(i)

#########################
### ENGINE DEFINITION  ##
#########################  			            
O.engines=O.engines+[
	BoxFactory(maxParticles=100000, extents=(l,l,0.1),center=(0,0,L-0.1),vMin=0,vMax=10.0,PSDsizes=(0.002, 0.0025, 0.003), PSDcum=(0.33, 0.66, 1.0), PSDcalculateMass=True, exactDiam=True,vAngle=0.,massFlowRate=30000.,normal=(0.,0.,-1.),label='factory',materialId=2),
	DomainLimiter(lo=(-l,-l,0),hi=(l,l,1),iterPeriod=200),
	#PyRunner(iterPeriod=1,command='T=max(T+2,700); print "T=",T',label="Tupdate"),
	#PyRunner(iterPeriod=1,dead=False,command='Apply_confiningpressure(T)'),
	NewtonIntegrator(damping=0.7,gravity=(0,0,-9.81),label='Newton'),
	#PyRunner(iterPeriod=1,command='normalstresses()'),
	
	]	