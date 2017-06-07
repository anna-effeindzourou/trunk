# -*- coding: utf-8
from yade import ymport, utils,pack,export, geom, plot,qt
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
def hMax(n):
	idHMax=0
	hMax=-1000000.0
	for i in O.bodies:
		h=i.state.pos[n]
		if (h>hMax) and (type(i.shape)==Sphere):
			hMax=h
			idHMax=i.id
	hMax=hMax+O.bodies[idHMax].shape.radius
	return (hMax)	
    
    
def hMin(n):
	idHMin=0
	hMin=100000.0
	for i in O.bodies:
		h=i.state.pos[n]
		if (h<hMin)and (type(i.shape)==Sphere):
			hMin=h
			idHMin=i.id
	hMin=hMin-O.bodies[idHMin].shape.radius
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
 
def sup():
	for i in O.bodies:
		if (type(i.shape)==Sphere) and (i.state.pos[2]>0.098):
			O.bodies.erase(i.id)    
  
def scalar(u,v):
	ps=u[0]*v[0]+u[1]*v[1]+u[2]*v[2]
	return ps

def cross(u,v):
	ps=Vector3(u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2] ,u[0]*v[1]-u[1]*v[0])
	return ps  
  

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
	]),
	InteractionLoop([
		Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_Wall_Sphere_ScGeom(),
		Ig2_Wall_PFacet_ScGeom()
		],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True),
		Ip2_FrictMat_FrictMat_MindlinPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(membrane=True),
		Law2_ScGridCoGeom_MindlinPhys_Mindlin(),
		Law2_ScGeom_MindlinPhys_Mindlin(),
		
		]
	),
]


######################
###  PROPERTIES    ###
######################
radius=0.0025*rm #r_s=0.001
sigma=80e3

#### Parameters of a rectangular grid ###
L=0.205 #length [m]
l=L/4.	#half width	(radius) [m]

nbL=60#number of nodes for the length	[#] doit etre paire 32
nbl=64	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!! 40

r=radius
color=[1,0,0]
nodesIds=[]
pfIds=[]
pfIds2=[]

top_boundary=[]
bottom_boundary=[]

####################
###  MATERIAL    ###
####################
poisson=0.28
E=2e11 #!#1e11,7.9e10

density=7.8e12##7.8e10
density_mem=1e12

frac1=1
frac2=1
Emem=E/frac1 #t1, 10
Ecyl=E/frac2#t1, 500

poisson_mem=0.49
frictionAngle=0.096
poro=38
threshold=0.009
O.tags['description']='triax_E_'+str(E)+'_Emem_factor_'+str(frac1)+'_Ecyl_factor'+str(frac2)+'_sigma_'+str(sigma)+'_threshold_'+str(threshold)+'_timestep_'+str(20)

print O.tags['description']
O.materials.append(CohFrictMat(young=Emem,poisson=poisson_mem,density=density_mem,frictionAngle=0,normalCohesion=1e100,shearCohesion=1e100,momentRotationLaw=False,alphaKr=0,label='NodeMat'))
O.materials.append(CohFrictMat(young=Emem,poisson=poisson_mem,density=density_mem,frictionAngle=0,normalCohesion=1e100,shearCohesion=1e100,momentRotationLaw=False,alphaKr=0,label='2NodeMat'))
O.materials.append(FrictMat(young=Emem,poisson=poisson_mem,density=density_mem,frictionAngle=0,label='memMat'))

#O.materials.append(FrictMat(young=100*E,poisson=poisson,density=density,frictionAngle=0,label='Wallmat'))

O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngle,label='Smat'))

##############################
###  SAMPLE GENERATION     ###
##############################
L=0.203 #high [m]



r_s=0.0025*1.5
#l=0.0505	# radius [m]
kw={'color':[0.8,0.8,0.8],'wire':False,'dynamic':True,'material':3}
predicate=inCylinder(centerBottom=Vector3(0,0,0), centerTop=Vector3(0,0,L), radius=l)	
sp=SpherePack()
sp=pack.regularHexa(predicate, radius=r_s, gap=0.00,**kw)
pile2=O.bodies.append(sp)
#################################
####  MEMBRANE GENERATION     ###
#################################
#Create all nodes first :
for i in range(0,nbL+1):
	for j in range(0,nbl):
		z=i*L/float(nbL)
		#print z
		y=(l+r)*sin(2*pi*j/float(nbl))
		x=(l+r)*cos(2*pi*j/float(nbl))
		
		if(i==0) or (i==nbL):
			nodesIds.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='2NodeMat',color=Vector3(0,0,0))) )
		else:
			nodesIds.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )
L=z

typ=0
wire=True
###Create PFacets
for i in range(0,nbL):
	for j in range(0,nbl-1):
			if  ((i%2==0) and  (j%2==0))or ((i%2!=0) and  (j%2!=0)):
				pfacetCreator3(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
				pfacetCreator3(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j+1],nodesIds[(i)*nbl+j+1],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
			else:
				#pass
				pfacetCreator3(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i)*nbl+j+1],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
				pfacetCreator3(nodesIds[i*nbl+j+1],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
for i in range(0,nbL,1):
	color=color
	j=nbl
	if  ((i%2==0) and  (j%2==0))or ((i%2!=0) and  (j%2!=0)):
		pfacetCreator3( nodesIds[i*nbl+j],nodesIds[(i-1)*nbl+j],nodesIds[(i+1)*nbl+j-1],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
		pfacetCreator3( nodesIds[(i)*nbl+j-1],nodesIds[(i+1)*nbl+j-1],nodesIds[(i-1)*nbl+j],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)

	else:
		pfacetCreator3( nodesIds[i*nbl+j],nodesIds[(i)*nbl+j-1],nodesIds[(i+1)*nbl+j-1],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
		pfacetCreator3( nodesIds[i*nbl+j],nodesIds[(i-1)*nbl+j],nodesIds[(i)*nbl+j-1],pfIds=pfIds,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
###########################
####### WALL GENERATION  ##
###########################

topringmem=[]
bottomringmem=[]

for i in range(0,nbl,1):
	topringmem.append(nodesIds[len(nodesIds)-1-i])
	bottomringmem.append(nodesIds[i])

topring=topringmem
#for j in range(0,nbl):
	#z=L+2*r
	#y=(l+r)*sin(2*pi*j/float(nbl))
	#x=(l+r)*cos(2*pi*j/float(nbl))
	#topring.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='2NodeMat',color=Vector3(1,1,0))) )
	
nodesIds.append( O.bodies.append(gridNode([0,0,z],r,wire=False,fixed=False,material='2NodeMat',color=Vector3(1,0,0))) )

pfIds2=[]
for i in range(0,nbl-1,1): 
	pfacetCreator3( nodesIds[len(nodesIds)-1],topring[i],topring[i+1],pfIds=pfIds2,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
pfacetCreator3( nodesIds[len(nodesIds)-1],topring[i+1],topring[0],pfIds=pfIds2,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
topring.append(nodesIds[len(nodesIds)-1])

pfIdsb2=[]
bottomring=bottomringmem
#for j in range(0,nbl):
	#z=-2*r
	#y=(l+r)*sin(2*pi*j/float(nbl))
	#x=(l+r)*cos(2*pi*j/float(nbl))
	#bottomring.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='2NodeMat',color=Vector3(1,1,0))) )
#z=-2*r	
nodesIds.append( O.bodies.append(gridNode([0,0,0],r,wire=False,fixed=False,material='2NodeMat',color=Vector3(1,0,0))) )

for i in range(0,nbl-1,1): 
	pfacetCreator3( nodesIds[len(nodesIds)-1],bottomring[i],bottomring[i+1],pfIds=pfIdsb2,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)
pfacetCreator3( nodesIds[len(nodesIds)-1],bottomring[0],bottomring[i+1],pfIds=pfIdsb2,wire=wire,material='memMat',Ecyl=Ecyl,color=color,mask=-1,fixed=False)

bottomring.append(nodesIds[len(nodesIds)-1])

tbodies=[]
for i in topring:
	tbodies.append(O.bodies[i])

topPlate,clumpIds=O.bodies.appendClumped(tbodies)
#O.bodies[topPlate].state.blockedDOFs='xyzXYZ'
#topPlate=topring[-1]

bbodies=[]
for i in bottomring:
	bbodies.append(O.bodies[i])

bottomPlate,clumpIdsb=O.bodies.appendClumped(bbodies)
#O.bodies[bottomPlate].state.blockedDOFs='xyzXYZ'

#bottomPlate=bottomring[-1]
pos0=O.bodies[topPlate].state.pos[2]


eps=L*0.01
middle=[]
for i in nodesIds:
	if(O.bodies[i].state.pos[2]>L/2.-eps)and(O.bodies[i].state.pos[2]<L/2.+eps):
		 middle.append(i)
		 O.bodies[i].shape.color=Vector3(1,1,1)
#print middle
m0=middle[0]
mm=middle[nbl/2]
O.bodies[m0].shape.color=Vector3(0,0,0)
O.bodies[mm].shape.color=Vector3(0,0,0)


###################
#### APPLY LOAD  ##
###################
normalVEL=0
S0=pi*l**2
normalSTRESS=sigma
shearing=False
sigmaN=0
sigmaNb=0

#### APPLY CONFINING PRESSURE  
					
step=1e-02
lim=int(1/step)
print 'lim = ', lim
threshold_sigma=0.1


normalVEL=-5e-04
target_porosity = 0.39
	


move=False
conf=False
confdone=False
load=False	
topload=False
def Apply_confiningT():
	global move,sigmaN,Fnt,topload,conf
	Fnt=O.forces.f(topPlate)[2]
	sigmaN=Fnt/S0
	Fnb=abs(O.forces.f(bottomPlate)[2])
	sigmaNb=Fnb/S0
	
	if(unbalancedForce()<threshold) and (conf==False):
		Newton.damping=0.2
		Applycpressure.dead=False
		confblock()
		conf=True
		move=True
		print 'start applying confining pressure, O.realtime (min) = ', O.realtime/60.

	if (move==True)and(confdone==True) and (unbalancedForce()<threshold):
		unblock()
		porosity_cor.dead=False
		O.pause()
		saveData()
		ApplycpressureT.dead=True
		print 'topPlate velocity = ', O.bodies[topPlate].state.vel[2]
		print 'end confining pressure,start porosity correction, O.realtime (min) = ', O.realtime/60.
count1=1
c2=0
p0=0
sigma0=0
#step=1e-05
def Apply_confiningpressure():
	global confdone,move,sigmaN,count1,c2,p0,sigma0	
	if(load==False):			
		  area= pi*l**2
		  Fapplied=area*sigma
		  normal=Vector3(0,0,-1)
		  if(step*count1<1):
			F=Fapplied*step*count1	
		  else:
			confdone=True
			F=Fapplied
		  O.forces.addF(topPlate,F*normal,permanent=False)
		  
	for i in pfIds: 
	#i=35
		n=[O.bodies[i].shape.node1.id,O.bodies[i].shape.node2.id,O.bodies[i].shape.node3.id]
		n1= O.bodies[i].shape.node1.state.pos
		n2= O.bodies[i].shape.node2.state.pos
		n3= O.bodies[i].shape.node3.state.pos
		
		e0 =n3 - n1
		e1 =n2 - n1
		e2 =n2 - n3
		P=(n1+n2+n3)/3
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

		V  = 2*pi*pow(0.0025,3)
		#print 'V = ',V
		dt = O.dt
		v0 = (2*dt*Fapplied)/(density*V)
		#print 'v0 = ',v0
		
		normal = -cross(e0,e1)
		normal=normal/normal.norm()

		#F=Fapplied
		p1n=Fapplied*p1
		p2n=Fapplied*p2
		p3n=Fapplied*p3
		
		if(step*count1<1):
			
			#print count1
			p1incr=p1n*step*count1
			p2incr=p2n*step*count1
			p3incr=p3n*step*count1
			if(i==pfIds[0]):
				p0=p1incr
				sigma0=p1incr/area
				
				#print count1
				
		else:
			confdone=True
			p1incr=p1n
			p2incr=p2n
			p3incr=p3n
			if(i==pfIds[0]):
				p0=p1incr
				sigma0=p1incr/area

		
		p1normal=p1incr*Vector3(-n1[0],-n1[1],0)
		p2normal=p2incr*Vector3(-n2[0],-n2[1],0)
		p3normal=p3incr*Vector3(-n3[0],-n3[1],0)
		

		pnormal=[p1normal,p2normal,p3normal]
		
		
		for j in range(0,3):
			if (n[j] in topring) or (n[j] in bottomring):
				pass
			else:
				O.forces.addF(n[j],pnormal[j],permanent=False)
	c2+=1			
	if(c2==1)and(count1*step<1):	
		count1+=1
		c2=0
			
			
#### MOVE TOP AND BOTTOM WALL 
O.bodies[topPlate].state.pos[2]=hMax(2)
postop=O.bodies[topPlate].state.pos[2]
posbot=O.bodies[bottomPlate].state.pos[2]
displ0=O.bodies[topPlate].state.pos[2]-O.bodies[bottomPlate].state.pos[2]


def porosity_correction():
	global frictionAngle,move,displ0,load,sigmaN
	zmax=hMax(2)
	zmin=hMin(2)
	V = S0*(zmax-zmin)
	poro=porosity(V)	
	
	#print 'start applying correcting porosity, O.iter = ', O.iter
	if(poro>target_porosity):
		frictionAngle = 0.95*frictionAngle
		for i in pile2:
			try:
				O.bodies[i].mat.frictionAngle=frictionAngle
			except:
				pass
	if (unbalancedForce()<threshold) and (poro<target_porosity) :
		print 'end correcting porosity, O.realtime (min) = ', O.realtime/60.
		move=False
		load=True
		moveWall(normalVEL)
		displ0=O.bodies[topPlate].state.pos[2]-O.bodies[bottomPlate].state.pos[2]
		print 'start of loading, O.realtime (min) = ', O.realtime/60.,', displ0 = ',displ0,', sigmaN = ',sigmaN
		porosity_cor.dead=True
		#O.run(1000,True)
		



def moveWall(v):
	O.bodies[topPlate].state.vel=(0,0,v)
	#for i in range(0,len(topring),1):
		#O.bodies[topring[i]].state.vel[2]=v
	#for i in range(0,len(topringmem),1):
		#O.bodies[topringmem[i]].state.vel[2]=v
#O.bodies[topring[i]].state.pos
###########################
##### MEMBRANE CONTROL   ##
########################### 	
epsr0=O.bodies[m0].state.pos[0]-O.bodies[mm].state.pos[0]
epsa0=O.bodies[topPlate].state.pos[2]-O.bodies[bottomPlate].state.pos[2]	
	
def blocklimit():
	for i in nodesIds:
		O.bodies[i].state.blockedDOFs='xyzXYZ'
#blocklimit()
######################
#####   SAVE DATA   ##
###################### 	
def dataCollector():
	zmax=hMax(2)
	zmin=hMin(2)
	V = S0*(zmax-zmin)
	#poro=porosity(V)	
	#F=O.forces.f(O.bodies[pfIds[45]].shape.node1.id)
	#print F
	S=pi*l**2
	Fnt=O.forces.f(topPlate)[2]
	sigmaN=Fnt/S0
	Fnb=abs(O.forces.f(bottomPlate)[2])
	sigmaNb=Fnb/S0
	pos=O.bodies[topPlate].state.pos[2]
	q=(sigmaNb-sigma)
	cui=(sigmaNb-sigma)/(sigmaNb+sigma)
	p=(sigmaNb+2*sigma)/2
	displ=O.bodies[topPlate].state.pos[2]-O.bodies[bottomPlate].state.pos[2]
	epsr=epsr0-(O.bodies[m0].state.pos[0]-O.bodies[mm].state.pos[0])
	epsa=epsa0-(O.bodies[topPlate].state.pos[2]-O.bodies[bottomPlate].state.pos[2])
	epsv=epsa+2*epsr
	if((displ0-displ)>0.008) and (load==True):
		O.pause()
		#print 'Real time = ', O.realtime
		print 'end of loading, O.realtime (min) = ', O.realtime/60.
		O.bodies[bottomPlate].state.vel=(0,0,0) 
		O.bodies[topPlate].state.vel=(0,0,0)
		saveData()
		#O.exitNoBacktrace()
	plot.addData(t=O.time,t2=O.time,displ=displ,d=(displ0-displ)/displ0,pos=pos,Fnt=Fnt,Fnb=Fnb,sigmaN=sigmaN,v=O.bodies[topPlate].state.vel[2],cui=cui,sigmaNb=sigmaNb,unbF=unbalancedForce(),p=p,q=q,p0=p0,sigma0=sigma0,epsv=epsv)
    
plot.plots={'t':('unbF','p0'),'t2':('epsv')}
plot.plot(subPlots=True)
    
def saveData():
	plot.saveDataTxt('data/'+O.tags['description']+'.dat.bz2',vars=('t','pos','displ','Fnt','Fnb','sigmaN','sigmaNb','cui','unbF','p','q','v','p0','sigma0','epsv'))

	
#def saveDonnees():
	#saveData()			
def block():
	O.bodies[topPlate].state.blockedDOFs='xyzXYZ'	
	O.bodies[bottomPlate].state.blockedDOFs='xyzXYZ'	
		
def confblock():
	O.bodies[topPlate].state.blockedDOFs='xy'	
	O.bodies[bottomPlate].state.blockedDOFs='xyz'		

def unblock():

	O.bodies[topPlate].state.blockedDOFs='xyz'	
	#O.bodies[bottomPlate].state.blockedDOFs='xyz'	

		
block()

def checkvelocity(val):
	tot=0
	for i in O.bodies:
		if(abs(i.state.vel[0])>val)or(abs(i.state.vel[1])>val) or(abs(i.state.vel[2])>val) or(abs(i.state.angVel[0])>val)or(abs(i.state.angVel[1])>val) or(abs(i.state.angVel[2])>val):
			tot+=1
	return tot
  
def checkvelocityid(val):
	for i in O.bodies:
		if(abs(i.state.vel[0])>val)or(abs(i.state.vel[1])>val) or(abs(i.state.vel[2])>val) or(abs(i.state.angVel[0])>val)or(abs(i.state.angVel[1])>val) or(abs(i.state.angVel[2])>val):
			print i.id, i.state.vel.norm(),i.state.angVel.norm()

def initvel():
	print 'unbalancedForce = ', unbalancedForce(),', checkvelocity(0.001) = ',checkvelocity(0.001),', sigmaN = ',sigmaN

###########################
##### ENGINE DEFINITION  ##
########################### 
g=0
Gl1_PFacet.wire=True

timestep_coeff=0.5
O.dt=timestep_coeff*PWaveTimeStep()
para=int(0.5/O.dt)
para=100
coll=100
print 'para = ',para
O.engines=O.engines+[
	PyRunner(iterPeriod=1,dead=True,command='Apply_confiningpressure()',label='Applycpressure'),
	PyRunner(iterPeriod=100,dead=False,command='Apply_confiningT()',label='ApplycpressureT'),
	PyRunner(iterPeriod=1000,dead=True,command='porosity_correction()',label='porosity_cor'),
	NewtonIntegrator(damping=0.9,gravity=(0,0,g),label='Newton'),GlobalStiffnessTimeStepper(timestepSafetyCoefficient=timestep_coeff),
	PyRunner(initRun=True,iterPeriod=coll,command='dataCollector()'),
	PyRunner(initRun=False,iterPeriod=para,command='saveData()'),
	VTKRecorder(iterPeriod=para,dead=True,initRun=True,fileName='paraview/'+O.tags['description']+'_',recorders=['spheres','velocity','intr','materialId','stress','bstresses'],label='parav'),
	]


parav.dead=False
qt.GLViewer.ortho=True
#  VISUALIZATION  

qt.Controller()
qtv = qt.View()
qtr = qt.Renderer()
qtr.light2=True        # back
qtr.lightPos=Vector3(1200,1500,500)
qtr.bgColor=[1,1,1]


O.run()