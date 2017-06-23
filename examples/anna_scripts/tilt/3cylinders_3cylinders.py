## encoding: utf-8
from yade import utils

O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_GridConnection_Aabb()]),
		InteractionLoop(
		[Ig2_GridNode_GridNode_GridNodeGeom6D(),Ig2_GridConnection_GridConnection_GridCoGridCoGeom()],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True),
		Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		Law2_GridCoGridCoGeom_FrictPhys_CundallStrack()
		]
	),
	NewtonIntegrator(gravity=(0,0,-9.81),damping=0.1,label='newton')
]

O.materials.append(CohFrictMat(young=1e7,poisson=1,density=1e3,frictionAngle=radians(30),normalCohesion=3e7,shearCohesion=3e7,momentRotationLaw=True,label='gridNodeMat'))
O.materials.append(FrictMat(young=1e7,poisson=1,density=1e3,frictionAngle=radians(30),label='gridConnectionMat'))

color=[255./255.,102./255.,0./255.]
###########################
#####   Cylinders      #####
###########################
factor=1
fixed = True
O.bodies.append(gridNode([-0.5,0,0],0.03,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 
O.bodies.append(gridNode([0.5,0,0],0.03,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 
O.bodies.append(gridNode([0,0,-0.5],0.03,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 

O.bodies.append( gridConnection(0,1,0.03,color=color) )	
O.bodies.append( gridConnection(1,2,0.03,color=color) )	
O.bodies.append( gridConnection(0,2,0.03,color=color) )	


fixed = False
O.bodies.append(gridNode([0,0,2*0.03+0.5],0.03,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 
O.bodies.append(gridNode([-0.5*factor,0,2*0.03],0.03,wire=False,fixed=fixed,material='gridNodeMat',color=color))
O.bodies.append(gridNode([0.5*factor,0,2*0.03],0.03,wire=False,fixed=fixed,material='gridNodeMat',color=color)) 


O.bodies.append( gridConnection(6,7,0.03,color=color) )
O.bodies.append( gridConnection(7,8,0.03,color=color) )
O.bodies.append( gridConnection(6,8,0.03,color=color) )

O.dt = 0.01*PWaveTimeStep()


def printInteractions():
	print 'Number of Real interactions =', O.interactions.countReal()
	for i in O.interactions:
		print "i.id1 = ",i.id1,", i.id2= ",i.id2

O.dt = 0.01*PWaveTimeStep()
O.step()
printInteractions()
O.saveTmp()

		