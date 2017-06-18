# -*- coding: utf-8
from yade import ymport, utils,pack,export, geom, plot
import gts,os

#### set False when running in batch mode
defaultTable = True
####-------------------------------------
####-------------------------------------

utils.readParamsFromTable(
rm = 0.33,
noTableOk = True
)
from yade.params.table import *
print 'rm=',rm

#################################
#####     FUNCTIONS          ####
#################################
def limitfinder():
	for b in O.bodies:
		if(b.state.pos[2]>=L-2*radius):
			if isinstance(b.shape,GridNode):
				top_boundary.append(b.id)
				b.shape.color=(0,0,1)
				b.state.blockedDOFs='z'
		if(b.state.pos[2]<0.1*radius ):
			if isinstance(b.shape,GridNode):
				bottom_boundary.append(b.id)
				b.state.blockedDOFs='z'
				b.shape.color=(0,0,1)	
				
########################
##  GENERATE FOLDERS  ##
########################
try:
	os.mkdir('data')
except:
	pass
try:
	os.mkdir('paraview')
except:
	pass
isBatch = runningInBatch()

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

######################
###  PROPERTIES    ###
######################
#radius=0.0025*rm
radius=0.0075*rm
sigma=-3e6

#### Parameters of a rectangular grid ###
L=0.205 #length [m]
l=0.101/2.	#half width	(radius) [m]
nbL=20#number of nodes for the length	[#] doit etre paire
nbl=24	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!

#nbL=36#number of nodes for the length	[#] doit etre paire
#nbl=44	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!

r=radius
color=[1,0,0]
nodesIds=[]
pfIds=[]


top_boundary=[]
bottom_boundary=[]

####################
###  MATERIAL    ###
####################
poisson=0.28
E=2*7.9e10*(1+poisson) ##1e11
density=7.8e8
density_mem=7.8e3
Et=E/1e3
Emem=E/1e1
frictionAngle=0.096
frictionAngleW=0.228

O.tags['description']='triaxial_rm_'+str(rm)+'_Emem_'+str(Emem)

O.materials.append(CohFrictMat(young=Emem,poisson=poisson,density=density_mem,frictionAngle=0,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=False,alphaKr=0,label='NodeMat'))
O.materials.append(FrictMat(young=Emem,poisson=poisson,density=density_mem,frictionAngle=0,label='memMat'))
O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngleW,label='Wallmat'))
O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngle,label='Smat'))
O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngleW,label='Wallmat'))


##############################
###  SAMPLE GENERATION     ###
##############################
kw={'color':[0.8,0.8,0.8],'wire':False,'dynamic':True,'material':3}
pile=ymport.text('spheres_7.5e-03.txt',**kw)
pile2=O.bodies.append(pile)

L=hMax(2)
#################################
####  MEMBRANE GENERATION     ###
#################################
#Create all nodes first :
for i in range(0,nbL+1):
	for j in range(0,nbl):
		z=i*L/float(nbL)
		y=(l+r)*sin(2*pi*j/float(nbl))
		x=(l+r)*cos(2*pi*j/float(nbl))
		nodesIds.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )

##Create connection between the nodes
for i in range(0,nbL+1):
	for j in range(0,nbl-1):
		O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[i*nbl+j+1],r,color=color,mask=5,material='memMat',Et=Et) )
for i in range(0,nbL,1):	
		for j in range(0,nbl):
			O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],r,color=color,mask=5,material='memMat',Et=Et) )
for i in range(-1,nbL):
	j=nbl
	O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color,mask=5,material='memMat',Et=Et) )
	
for i in range(0,nbL):
	for j in range(0,nbl-1):
			if (j%2==0):
				O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j+1],r,color=color,mask=5,material='memMat',Et=Et) )
			else:
				O.bodies.append( gridConnection(nodesIds[(i+1)*nbl+j],nodesIds[i*nbl+j+1],r,color=color,mask=5,material='memMat',Et=Et) )
for i in range(0,nbL):
	j=nbl
	O.bodies.append( gridConnection(nodesIds[(i-1)*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color,mask=5,material='memMat',Et=Et) )		


###Create PFacets
for i in range(0,nbL):
	for j in range(0,nbl-1):
			if (j%2==0):
				pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],color=color,mask=5,material='memMat')))
				pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j+1],nodesIds[(i)*nbl+j+1],color=color,mask=5,material='memMat')))
			else:
				pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i)*nbl+j+1],color=color,mask=5,material='memMat')))
				pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j+1],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],color=color,mask=5,material='memMat')))
	
for i in range(0,nbL,1):
	j=nbl
	pfIds.append(O.bodies.append(pfacet( nodesIds[i*nbl+j],nodesIds[(i-1)*nbl+j],nodesIds[(i+1)*nbl+j-1],color=color,material='memMat' )))
	pfIds.append(O.bodies.append(pfacet( nodesIds[(i)*nbl+j-1],nodesIds[(i+1)*nbl+j-1],nodesIds[(i-1)*nbl+j],color=color,material='memMat' )))


#########################
##### WALL GENERATION  ##
#########################
topPlate=utils.wall(position=hMax(2)+radius,sense=0, axis=2,color=Vector3(1,0,0),material='Wallmat')
O.bodies.append(topPlate)

bottomPlate=utils.wall(position=0,sense=0, axis=2,color=Vector3(1,0,0),material='Wallmat')
O.bodies.append(bottomPlate)

###################
#### APPLY LOAD  ##
###################
normalVEL=0
S0=pi*l**2
normalSTRESS=sigma
shearing=False
sigmaN=0


#### APPLY CONFINING PRESSURE  
def Apply_load():
	global sigmaN, Fn
	Fn=abs(O.forces.f(topPlate.id)[2])
	sigmaN=Fn/S0
	if abs((normalSTRESS-sigmaN)/normalSTRESS)<0.001:
		topPlate.state.vel[2]=0

					
def Apply_confiningpressure():
	for i in pfIds:
		e0 =O.bodies[i].shape.node3.state.pos - O.bodies[i].shape.node1.state.pos
		e1 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node1.state.pos
		e2 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node3.state.pos
		P=(O.bodies[i].shape.node1.state.pos+O.bodies[i].shape.node2.state.pos+O.bodies[i].shape.node3.state.pos)/3
		v0 = e0
		v1 = e1
		v2 = P - O.bodies[i].shape.node1.state.pos
		
		##// Compute dot products
		dot00 = scalar(v0,v0)
		dot01 = scalar(v0,v1)
		dot02 = scalar(v0,v2)
		dot11 = scalar(v1,v1)
		dot12 = scalar(v1,v2)
		
		##// Compute the barycentric coordinates of the projection P
		invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
		p1 = (dot11 * dot02 - dot01 * dot12) * invDenom
		p2 = (dot00 * dot12 - dot01 * dot02) * invDenom
		p3 = 1-p1-p2

		a = sqrt(scalar(e0,e0))
		b = sqrt(scalar(e1,e1))
		c = sqrt(scalar(e2,e2))

		s=0.5*(a+b+c)
		area= sqrt(s*(s-a)*(s-b)*(s-c))
		
		Fapplied=area*sigma
		
		normal = cross(e0,e1)	
		
		normal=normal/normal.norm()
		
		F=Fapplied
		p1normal=F*p1*normal
		p2normal=F*p2*normal
		p3normal=F*p3*normal
			
		O.forces.addF(O.bodies[i].shape.node1.id,p1normal,permanent=False)
		O.forces.addF(O.bodies[i].shape.node2.id,p2normal,permanent=False)
		O.forces.addF(O.bodies[i].shape.node3.id,p3normal,permanent=False)


sigma3=0
def check_confiningpressure():
	global sigma3
	sigma3=0
	for i in pfIds:
		e0 =O.bodies[i].shape.node3.state.pos - O.bodies[i].shape.node1.state.pos
		e1 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node1.state.pos
		e2 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node3.state.pos

		a = sqrt(scalar(e0,e0))
		b = sqrt(scalar(e1,e1))
		c = sqrt(scalar(e2,e2))

		s=0.5*(a+b+c)
		area= sqrt(s*(s-a)*(s-b)*(s-c))
		F=(O.forces.f(O.bodies[i].shape.node1.id) + O.forces.f(O.bodies[i].shape.node2.id)+O.forces.f(O.bodies[i].shape.node3.id)).norm()
		sigma3=sigma3+F/area
	return sigma3


pos=topPlate.state.pos[2]

def dataCollector():
	global pos,p,q,sigma1
	S=pi*l**2
	Fnt=O.forces.f(topPlate.id)[2]
	Fnb=O.forces.f(bottomPlate.id)[2]
	sigma1=Fnt/S
	sigma3=check_confiningpressure()
	pos=topPlate.state.pos[2]
	q=(sigma1-3e6)
	p=(sigma1+2*3e6)/2
	
	plot.addData(t=O.time,t2=O.time,pos=pos,Fnt=Fnt,Fnb=Fnb,sigma1=sigma1,sigma3=sigma3,unbF=unbalancedForce(),p=p,q=q)
    
    
def saveData():
	plot.saveDataTxt('data/'+O.tags['description']+'.dat',vars=('t','pos','Fnt','Fnb','sigma1','sigma3','unbF','p','q'))

plot.plots={'t':('Fnb'),'t2':('Fnt')}

	
def saveDonnees():
	saveData()
limitfinder()		
#### MOVE TOP AND BOTTOM WALL 
v=1.7e-03
#v=1.7e-05
def moveWall(v):
	topPlate.state.vel=(0,0,-v)


###########################
##### ENGINE DEFINITION  ##
########################### 
#g=-9.81
g=0

iter_simu=int(24/(0.5*PWaveTimeStep()))*10
para=int(iter_simu/100)

O.dt=0.5*PWaveTimeStep()

O.engines=O.engines+[
	PyRunner(iterPeriod=1,initRun=True,dead=False,command='Apply_load()',label='Applyload'),
	PyRunner(iterPeriod=1,dead=False,command='Apply_confiningpressure()',label='Applycpressure'),
	NewtonIntegrator(damping=0.7,gravity=(0,0,g),label='Newton'),
	PyRunner(initRun=True,iterPeriod=1,command='dataCollector()'),
	PyRunner(initRun=True,iterPeriod=1000,command='saveDonnees()'),
	VTKRecorder(iterPeriod=para,initRun=True,fileName='paraview/'+O.tags['description']+'_',recorders=['spheres','velocity','intr']),
	]	


plot.plot(noShow=False, subPlots=True)

#O.run(5000,True)
#moveWall(v)
#O.run(iter_simu)
