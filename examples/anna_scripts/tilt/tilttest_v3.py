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
a=400e-3	 # side dimension 
h=150e-3	 # height 


#MATERIAL PROPERTIES
frictionAngle=radians(35)
density=1500
Kn=25e7 #normal stiffness Plassiard thesis
Ks=5e7 #tangential stiffness Plassiard thesis


#PARTICLE SIZE
Rs=0.0322*1# mean particle radius
Rf=0.18 # dispersion (Rs±Rf*Rs)
nSpheres=2000# number of particles

#ENGINES PARAMETERS
g=9.81
damping_coeff = 0.5

#MATERIAL'S PARAMETERS
#box
young_box=1e6
poison_box=0.25
friction_box=radians(5)
density_box=6e4

#gravel
G=40e5
poisson_g=0.25
young_g=30e5
m=1.878
v=(4*pi*pow(Rs,3))/3

#TIME STEP
dt=1.e-7
#TAG
O.tags['description']='tilttest_rmean'+str(Rs)
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

clumpbox=True

############################
### DEFINING MATERIAL    ###
############################
O.materials.append(CohFrictMat(young=young_box,poisson=poison_box,density=density_box,frictionAngle=friction_box,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=False,alphaKr=5,label='NodeMat'))

boxMat = O.materials.append(FrictMat(young=young_box ,poisson=poison_box,frictionAngle=friction_box,density=density_box))
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

mesh=2
nodesIds1=[]

nodesIds1.append( O.bodies.append(gridNode([-x,-y,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds1.append( O.bodies.append(gridNode([x,-y,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )
if(mesh==2):
	nodesIds1.append( O.bodies.append(gridNode([0,-y,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds1.append( O.bodies.append(gridNode([x,-y,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds1.append( O.bodies.append(gridNode([-x,-y,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )
if(mesh==2):
	nodesIds1.append( O.bodies.append(gridNode([-x,x,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds1.append( O.bodies.append(gridNode([-x,2*x,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds1.append( O.bodies.append(gridNode([x,2*x,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

if(mesh==2):
	nodesIds1.append( O.bodies.append(gridNode([0,2*x,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds1.append( O.bodies.append(gridNode([x,2*x,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds1.append( O.bodies.append(gridNode([-x,2*x,2*h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

if(mesh==2):
	nodesIds1.append( O.bodies.append(gridNode([x,x,h+fac*rWall],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

cylIds=[]
pfIds=[]

Ecyl=100*young_box

if(mesh==1):
	pfIds.append(pfacetCreator3(nodesIds1[0],nodesIds1[1],nodesIds1[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[0],nodesIds1[2],nodesIds1[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[5],nodesIds1[6],nodesIds1[1],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[6],nodesIds1[2],nodesIds1[1],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))

	pfIds.append(pfacetCreator3(nodesIds1[7],nodesIds1[5],nodesIds1[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[7],nodesIds1[6],nodesIds1[5],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))



	pfIds.append(pfacetCreator3(nodesIds1[0],nodesIds1[4],nodesIds1[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[4],nodesIds1[7],nodesIds1[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))

if(mesh==2):
	pfIds.append(pfacetCreator3(nodesIds1[0],nodesIds1[4],nodesIds1[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[2],nodesIds1[4],nodesIds1[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[2],nodesIds1[3],nodesIds1[1],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))

	pfIds.append(pfacetCreator3(nodesIds1[0],nodesIds1[5],nodesIds1[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[5],nodesIds1[10],nodesIds1[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[5],nodesIds1[6],nodesIds1[10],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	
	
	
	pfIds.append(pfacetCreator3(nodesIds1[7],nodesIds1[8],nodesIds1[9],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[9],nodesIds1[8],nodesIds1[10],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[8],nodesIds1[6],nodesIds1[10],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	
	
	
	pfIds.append(pfacetCreator3(nodesIds1[1],nodesIds1[11],nodesIds1[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[3],nodesIds1[11],nodesIds1[9],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds1[11],nodesIds1[7],nodesIds1[9],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	




nbodies=[]


if(clumpbox==True):
	for i in nodesIds1:
		nbodies.append(O.bodies[i])

	clump,clumpIds=O.bodies.appendClumped(nbodies)
	O.bodies[clump].state.blockedDOFs='xyzXYZ'
else:
	for i in nodesIds1:
		O.bodies[i].state.blockedDOFs='xyzXYZ'


fixed=True
nodesIds=[]
nodesIds.append( O.bodies.append(gridNode([-x,-y,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,-y,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )
if(mesh==2):
	nodesIds.append( O.bodies.append(gridNode([0,-y,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,-y,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,-y,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )
if(mesh==2):
	nodesIds.append( O.bodies.append(gridNode([-x,x,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,2*x,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,2*x,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

if(mesh==2):
	nodesIds.append( O.bodies.append(gridNode([0,2*x,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

nodesIds.append( O.bodies.append(gridNode([x,2*x,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )


nodesIds.append( O.bodies.append(gridNode([-x,2*x,h],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )

if(mesh==2):
	nodesIds.append( O.bodies.append(gridNode([x,x,0],rWall,wire=False,fixed=fixed,material='NodeMat',color=color)) )
cylIds=[]
pfIds=[]

if(mesh==1):
	pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[1],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[2],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))

	pfIds.append(pfacetCreator3(nodesIds[1],nodesIds[5],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[5],nodesIds[6],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))

	pfIds.append(pfacetCreator3(nodesIds[7],nodesIds[6],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[6],nodesIds[5],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))

	pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[4],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[4],nodesIds[7],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))

	pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[1],nodesIds[5],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[4],nodesIds[0],nodesIds[5],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=fixed,mask=-1))

if(mesh==2):
	pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[4],nodesIds[2],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[2],nodesIds[4],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[2],nodesIds[3],nodesIds[1],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))

	pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[5],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[5],nodesIds[10],nodesIds[4],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[5],nodesIds[6],nodesIds[10],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	
	
	
	pfIds.append(pfacetCreator3(nodesIds[7],nodesIds[8],nodesIds[9],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[9],nodesIds[8],nodesIds[10],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[8],nodesIds[6],nodesIds[10],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	
	
	
	pfIds.append(pfacetCreator3(nodesIds[1],nodesIds[11],nodesIds[3],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[3],nodesIds[11],nodesIds[9],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[11],nodesIds[7],nodesIds[9],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
		#bottom
	pfIds.append(pfacetCreator3(nodesIds[0],nodesIds[2],nodesIds[5],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[2],nodesIds[5],nodesIds[8],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[5],nodesIds[6],nodesIds[8],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	
	
	
	pfIds.append(pfacetCreator3(nodesIds[8],nodesIds[11],nodesIds[7],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[2],nodesIds[11],nodesIds[8],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	pfIds.append(pfacetCreator3(nodesIds[2],nodesIds[1],nodesIds[11],cylIds=cylIds,pfIds=pfIds,wire=False,material=boxMat,Ecyl=Ecyl,color=color,fixed=True,mask=-1))
	

################
##   SAMPLE  ###
################
#k1={'mask':3}
#sp=pack.SpherePack()
#sp.makeCloud((-x+2*rWall,-y+2*rWall,2*rWall),(x-2*rWall,2*x-2*rWall,6.*h),rMean=Rs,rRelFuzz=Rf)
#sp.toSimulation(**k1)



############################
###   DEFINING ENGINES   ###
############################
angularVelocity=0.05
idsToRotate=nodesIds
#idsToRotate.append(bottomBox)
O.dt=.5*utils.PWaveTimeStep()
O.engines+=[

	NewtonIntegrator(gravity=(0,0,-9.81),damping=0.9,label='Newton_Integrator'),
]

init=0
para=1000
coll=100


def test():
	global init
	O.run(100000,True)
	Newton_Integrator.damping=0.2
	if(clumpbox==True):
		O.bodies[clump].state.blockedDOFs='xyz'
	else:
		for i in nodesIds1:
			O.bodies[i].state.blockedDOFs='xyz'	
	O.run(10000,True)
	
	
def rotate():	
	init=O.iter
	#Newton_Integrator.damping=0.2
	if(clumpbox==True):
		O.bodies[clump].state.blockedDOFs='x'
	else:
		for i in nodesIds1:
			O.bodies[i].state.blockedDOFs='x'	
			
			
	
	O.engines+=[
		RotationEngine(ids=idsToRotate,angularVelocity=angularVelocity,rotateAroundZero=True,zeroPoint=(-rWall,0,0),rotationAxis=(1,0,0),label='rotationEngine'),
		VTKRecorder(iterPeriod=para,dead=False,initRun=True,fileName='paraview/'+O.tags['description']+'_',recorders=['spheres','velocity','intr','materialId'],label='parav'),
		PyRunner(initRun=True,iterPeriod=coll,command='dataCollector()')
	]
	
def stop():
	while(unbalancedForce()<0.3):
		rotationEngine.angularVelocity=angularVelocity
		O.run(500,True)
		rotationEngine.angularVelocity=0
		O.run(500,True)
	final=O.iter
	alpha=angularVelocity*(final-init)*O.dt
	alpha=degrees(alpha)
	print 'alpha = ', alpha
	
#####################
#####   SAVE DATA   ##
###################### 	
def dataCollector():
	unbF=unbalancedForce()
	#print unbF
	plot.addData(t1=O.time,t2=O.iter,unbF=unbF,y1=O.bodies[0].state.pos[1],y2=O.bodies[1].state.pos[1],y3=O.bodies[6].state.pos[1],y4=O.bodies[7].state.pos[1],vy=O.bodies[nodesIds1[0]].state.vel[1])
    
plot.plots={'t1':('unbF'),'t2':('y1','y2',None,'y3','y4')}
plot.plot(subPlots=True)
    
def saveData():
	plot.saveDataTxt(O.tags['description']+'.dat.bz2',vars=('t1','t2','unbF','y1','y2','y3','y4','vy'))


Gl1_PFacet.wire=True
nodesIdsbot=[0,1,4,5]
def forces():
	for i in nodesIdsbot:
		print O.forces.f(i)[1],i
		print O.forces.f(i),i
		#print O.forces.f(i).norm()/O.forces.f(0).norm()
############################
###   VISUALIZATION      ###
############################
renderer = qt.Renderer()
qt.View()
O.saveTmp()
#O.run()


def getAngle():
	alpha=angularVelocity*(final-init)*O.dt
	alpha=degrees(alpha)
	print 'alpha = ', alpha
