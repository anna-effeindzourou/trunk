# encoding: utf-8
from yade import utils,pack,geom,qt
import matplotlib
from pylab import *
import math


O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		
	#Bo1_GridConnection_Aabb(),
	Bo1_PFacet_Aabb(),
	
	],sortThenCollide=True),
		InteractionLoop(
		[Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_Sphere_GridConnection_ScGridCoGeom(),
		Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		Ig2_GridConnection_PFacet_ScGeom(),
		Ig2_PFacet_PFacet_ScGeom()
		
		],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True),
		Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		Law2_ScGridCoGeom_FrictPhys_CundallStrack(),
		Law2_GridCoGridCoGeom_FrictPhys_CundallStrack()
		]
	),
	NewtonIntegrator(gravity=(0,0,-9.81),damping=0.1,label='newton')
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


###########################
#####   Cylinders      #####
###########################
factor=0.5
fixed = True
test=3

if(test==1):
	nodesIds.append( O.bodies.append(gridNode([-0.5,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0.5,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0,0,-0.5],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )

	cylIds.append(O.bodies.append( gridConnection(0,1,r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(1,2,r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(0,2,r,color=color) )	)

	pfacetCreator3(0,2,1,cylIds=[],pfIds=[],wire=False,material=-1,Ecyl=None,color=color)

	fixed = False
	nodesIds.append( O.bodies.append(gridNode([0,0,3*r+0.5],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([-0.5*factor,0,3*r],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0.5*factor,0,3*r],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )


	cylIds.append(O.bodies.append( gridConnection(nodesIds[3],nodesIds[4],r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(nodesIds[4],nodesIds[5],r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(nodesIds[3],nodesIds[5],r,color=color) )	)

	pfacetCreator3(nodesIds[3],nodesIds[4],nodesIds[5],cylIds=[],pfIds=[],wire=False,material=-1,Ecyl=None,color=color)
elif(test==2):
	nodesIds.append( O.bodies.append(gridNode([-0.5,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0.5,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0,-0.5,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )

	cylIds.append(O.bodies.append( gridConnection(0,1,r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(1,2,r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(0,2,r,color=color) )	)

	pfacetCreator3(0,2,1,cylIds=[],pfIds=[],wire=False,material=-1,Ecyl=None,color=color)

	fixed = False
	nodesIds.append( O.bodies.append(gridNode([0,-0.2,3*r+0.5],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([-0.5*factor,-0.2,3*r],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0.5*factor,-0.2,3*r],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )


	cylIds.append(O.bodies.append( gridConnection(nodesIds[3],nodesIds[4],r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(nodesIds[4],nodesIds[5],r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(nodesIds[3],nodesIds[5],r,color=color) )	)

	pfacetCreator3(nodesIds[3],nodesIds[4],nodesIds[5],cylIds=[],pfIds=[],wire=False,material=-1,Ecyl=None,color=color)
else:
	nodesIds.append( O.bodies.append(gridNode([-0.5,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0.5,0,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([0,-0.5,0],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )

	cylIds.append(O.bodies.append( gridConnection(0,1,r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(1,2,r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(0,2,r,color=color) )	)

	pfacetCreator3(0,2,1,cylIds=[],pfIds=[],wire=False,material=-1,Ecyl=None,color=color)

	fixed = False
	nodesIds.append( O.bodies.append(gridNode([0,-0.4,0.2],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([-1*factor,-0.2,0.2],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )
	nodesIds.append( O.bodies.append(gridNode([1*factor,-0.2,0.2],r,wire=False,fixed=fixed,material='gridNodeMat',color=color)) )


	cylIds.append(O.bodies.append( gridConnection(nodesIds[3],nodesIds[4],r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(nodesIds[4],nodesIds[5],r,color=color) )	)
	cylIds.append(O.bodies.append( gridConnection(nodesIds[3],nodesIds[5],r,color=color) )	)

	pfacetCreator3(nodesIds[3],nodesIds[4],nodesIds[5],cylIds=[],pfIds=[],wire=False,material=-1,Ecyl=None,color=color)


def printInteractions():
	for i in O.interactions:
		print "i.id1 = ",i.id1,", i.id2= ",i.id2
		F=i.phys.normalForce+i.phys.shearForce
		print F.norm()
O.dt = 0.01*PWaveTimeStep()
O.saveTmp()
