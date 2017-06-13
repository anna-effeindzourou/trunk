
from yade import pack,ymport,qt
import gts, os.path, locale

locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')   #gts is locale-dependend.  If, for example, german locale is used, gts.read()-function does not import floats normally

'''
if you get "Error: unsupported locale setting"
-> type as root: "dpkg-reconfigure locales"
-> choose "en_US.UTF-8" (press space to choose)
'''

################
### ENGINES  ###
################
g=0

O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_PFacet_Aabb(), 
		Bo1_Sphere_Aabb(), 
		Bo1_GridConnection_Aabb(), 
	]),
	
	InteractionLoop([
		Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		Ig2_Sphere_GridConnection_ScGridCoGeom(),
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		#Ig2_PFacet_PFacet_ScGeom(),
		],
		
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=False),
		Ip2_FrictMat_FrictMat_FrictPhys()
		],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		Law2_ScGeom_FrictPhys_CundallStrack(),
		Law2_ScGridCoGeom_FrictPhys_CundallStrack(),
		Law2_GridCoGridCoGeom_FrictPhys_CundallStrack()
		]
	),
	NewtonIntegrator(gravity=(0,0,g),damping=0.,label='newton')
]


nodesIds=[]
cylIds=[]
pfIds=[]
wire=False
fixed=False
color=[1,0,0]
young=1e-02
angle=20
radius=0.00001
#3_5_8_15
O.materials.append(CohFrictMat(young=young*1e6,poisson=0.3,density=2650,frictionAngle=radians(angle),normalCohesion=1e100,shearCohesion=1e100,momentRotationLaw=True,label='spheremat2'))
O.materials.append(FrictMat(young=young*1e6,poisson=0.3,density=2650,frictionAngle=radians(angle),label='boxmat2'))





########################################
### GENERATE THE WALL AND THE SPHERE ###
########################################
nodesIds2= []
r=0.0005
#WALL

#a=0.25 #horse
#z=-0.08#horse


a=1
z=-0.57
#nodesIds2.append( O.bodies.append(gridNode([-a,-a,z],r,wire=False,fixed=True,material='spheremat2',color=color)) )
#nodesIds2.append( O.bodies.append(gridNode([a,-a,z],r,wire=False,fixed=True,material='spheremat2',color=color)) )
#nodesIds2.append( O.bodies.append(gridNode([-a,a,z],r,wire=False,fixed=True,material='spheremat2',color=color)) )
#nodesIds2.append( O.bodies.append(gridNode([a,a,z],r,wire=False,fixed=True,material='spheremat2',color=color)) )
#O.bodies.append( gridConnection(0, 1,r,color=color,material='spheremat2') )
#O.bodies.append( gridConnection(2,3,r,color=color,material='spheremat2') )
#O.bodies.append( gridConnection(2,1,r,color=color,material='spheremat2') )
#O.bodies.append( gridConnection(2,0,r,color=color,material='spheremat2') )
#O.bodies.append( gridConnection(3,1,r,color=color,material='spheremat2') )
#O.bodies.append(pfacet(2,1,0,wire=False,material='boxmat2',color=color))
#O.bodies.append(pfacet(2,3,1,wire=False,material='boxmat2',color=color))






color=[0,0,1]

#O.bodies.append([ sphere(center=(0,0,0.2),radius=.05,fixed=False), sphere((0,0,0.3),.05,fixed=False)])
#O.bodies[0].state.vel=Vector3(0,0,-1)

#ymport.gtsPFacet('head.gts',shift=0,scale=1.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=wire,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)


#ymport.gtsPFacet('bunny.gts',shift=0,scale=1.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=wire,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)

#ymport.gtsPFacet('tie.gts',shift=0,scale=1.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=wire,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)

#ymport.gtsPFacet('seashell.gts',shift=0,scale=5.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=wire,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)




#O.bodies.append([ sphere(center=(0,0,0.2),radius=.05,fixed=False), sphere((0,0,0.3),.05,fixed=False)])
#O.bodies[0].state.vel=Vector3(0,0,-1)
#O.bodies[1].state.vel=Vector3(0,0,-1)
#fixed=False
#wire=True
#ymport.gtsPFacet('goblet.gts',shift=0,scale=1.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=wire,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)
#ymport.gtsPFacet('goblet.gts',shift=1,scale=1.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=True,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)

ymport.gtsPFacet('donut.gts',shift=0,scale=1.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=wire,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)
#ymport.gtsPFacet('seashell.gts',shift=1,scale=1.0,radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=True,fixed=fixed,materialNodes='spheremat2',material='boxmat2',color=color)


#ymport.gmshPFacet(meshfile="donut.mesh",shift=Vector3.Zero,scale=1.0,orientation=Quaternion.Identity,nodesIds=[],cylIds=[],pfIds=[],radius=1.0,wire=False,fixed=True,materialNodes=-1,material=-1,color=None)

O.tags['description']='sphere_n1_'+str(young)




O.dt=PWaveTimeStep()
##########
## VIEW ##
##########

qt.Controller()
qtv = qt.View()
qtr = qt.Renderer()
qtr.light2=True        # back
qtr.lightPos=Vector3(1200,1500,500)
qtr.bgColor=[1,1,1]
qtv.ortho=True

