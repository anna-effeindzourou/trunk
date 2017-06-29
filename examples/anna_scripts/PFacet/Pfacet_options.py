# encoding: utf-8
from yade import utils,pack,geom,qt
import matplotlib
from pylab import *
import math


O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_Sphere_Aabb(),
		Bo1_Wall_Aabb(),
		Bo1_PFacet_Aabb(),
	],sortThenCollide=True),
		InteractionLoop(
		[Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		Ig2_Wall_PFacet_ScGeom(),Ig2_Wall_Sphere_ScGeom()
		],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True),
		Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		Law2_ScGeom_FrictPhys_CundallStrack(),
		Law2_ScGridCoGeom_FrictPhys_CundallStrack(),
		Law2_GridCoGridCoGeom_FrictPhys_CundallStrack()
		]
	),
	NewtonIntegrator(gravity=(0,-9.81,0),damping=0.1,label='newton')
]


poisson=1
angle=30
young=1e7
O.materials.append(CohFrictMat(young=young,poisson=poisson,density=1e2,frictionAngle=radians(angle),normalCohesion=3e7,shearCohesion=3e7,momentRotationLaw=True,label='gridNodeMat'))
O.materials.append(FrictMat(young=young,poisson=poisson,density=1e2,frictionAngle=radians(angle),label='gridConnectionMat'))

nodesIds=[]
cylIds=[]

color=[255./255.,102./255.,0./255.]
#position of the node in the middle
a=0.00
r=0.03

fixed = False
###########################
#####   PFacet        #####
###########################
## Option 1: pfacet(id1,id2,id3) -> based on 3 nodes already connected via 3 cylinders

nodesIds.append( O.bodies.append(gridNode([0,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
nodesIds.append( O.bodies.append(gridNode([1,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
nodesIds.append( O.bodies.append(gridNode([0.5,1,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )

cylIds.append(O.bodies.append( gridConnection(0,1,r,color=color) )	)
cylIds.append(O.bodies.append( gridConnection(1,2,r,color=color) ))
cylIds.append(O.bodies.append( gridConnection(2,0,r,color=color) ))
	
O.bodies.append( pfacet(nodesIds[0],nodesIds[1],nodesIds[2],wire=False,color=color,highlight=False,material=O.materials[1])
)


##################################
######   pfacetCreator1        ###
##################################
## Option 2: pfacetCreator1(vertices) -> based on 3 vertices

v1=Vector3(2,0,0)
v2=Vector3(3,0,0)
v3=Vector3(2.5,1,0)
vertices=[]
vertices.append(v1)
vertices.append(v2)
vertices.append(v3)
pfacetCreator1(vertices,r,nodesIds=[],cylIds=[],pfIds=[],wire=False,color=color,fixed=fixed,materialNodes='gridNodeMat',material='gridConnectionMat')

#################################
#####   pfacetCreator2        ###
#################################
## Option 3: pfacetCreator2(id1,id2,vertex) -> based on 2 nodes connected via a cylinder and a vertex

nodesIds.append( O.bodies.append(gridNode([4,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
nodesIds.append( O.bodies.append(gridNode([5,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
vertex=Vector3(4.5,1,0)
cylIds.append(O.bodies.append( gridConnection(nodesIds[3],nodesIds[4],r,color=color) )	)


pfacetCreator2(nodesIds[3],nodesIds[4],vertex,r,nodesIds=nodesIds,wire=True,materialNodes='gridNodeMat',material='gridConnectionMat',color=color,fixed=fixed)

##################################
######   pfacetCreator3        ###
##################################
## Option 4: pfacetCreator3(id1,id2,id3) -> based on 3 nodes

a = O.bodies.append(gridNode([6,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 
b = O.bodies.append(gridNode([7,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 
c = O.bodies.append(gridNode([6.5,1,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 
 
pfacetCreator3(a,b,c,cylIds=[],pfIds=[],wire=False,material=-1,Ecyl=None,color=color)

#####################
#####   Wall      ###
#####################
##
Plate=utils.wall(position=-1,sense=0, axis=1,color=Vector3(1,0,0),material='gridConnectionMat')
O.bodies.append(Plate)


O.dt = 0.01*PWaveTimeStep()
O.saveTmp()
