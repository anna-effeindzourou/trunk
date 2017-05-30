# enoding: utf-8
from yade import utils,pack,geom,qt,export,ymport

#try:
    #os.mkdir('data_pullout')
#except:
    #pass
#try:
    #os.mkdir('paraview')
#except:
    #pass



def hMax(n):
	idHMax=0
	hMax=-1000000.0
	for i in O.bodies:
		h=i.state.pos[n]
		if (h>hMax):
			hMax=h
			idHMax=i.id
	return (hMax)	


def hMaxS(n):
	idHMax=0
	hMax=-1000000.0
	for i in O.bodies:
		if (type(i.shape)==Sphere):
			h=i.state.pos[n]
			if (h>hMax):
				hMax=h
				idHMax=i.id
	return (hMax)	

    
def hMin(n):
    idHMin=0
    hMin=100000.0
    for i in O.bodies:
	h=i.state.pos[n]
	if (h<hMin):
	   hMin=h
	   idHMin=i.id
    return (hMin) 
 
 
O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_Sphere_Aabb(),
		Bo1_Box_Aabb(),Bo1_Wall_Aabb(),
		Bo1_PFacet_Aabb(),
	]),
	InteractionLoop(
		[
		Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		Ig2_Wall_Sphere_ScGeom(),
		Ig2_Box_Sphere_ScGeom()
		],
		[
		Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=False,label='ipf'),
		Ip2_FrictMat_FrictMat_FrictPhys()
		],
		[
		Law2_ScGeom_FrictPhys_CundallStrack(),
		Law2_ScGeom6D_CohFrictPhys_CohesionMoment() ,
		Law2_ScGridCoGeom_FrictPhys_CundallStrack()
		]
	),
]


O.engines[1].avoidSelfInteractionMask=4
angle=20
O.tags['description']='mem_'+str(angle)
print O.tags['description']
O.materials.append(CohFrictMat(young=5e6,poisson=0.3,density=2650,frictionAngle=radians(angle),normalCohesion=3e100,shearCohesion=3e100,momentRotationLaw=False,label='gridNodeMat'))
O.materials.append(FrictMat(young=5e6,poisson=0.3,density=2650,frictionAngle=radians(angle),label='gridCoMat'))

### Parameters of a rectangular grid ###
L=0.05
H=L
nH=6
nh=nH-1
nL=2*nH
nS=nL
if(nL%2==0): nS=nL-1
z=0
y=0
x=0

xy=0.05


r=0.001
color=[255./255.,102./255.,0./255.]
zm=-2*r


#STEP 1: GENERATION OF THE MEMBRANE
nodesIds=[]
nodesIds1=[]
cylIds=[]
PFacetIds=[]
		
zm=0
for i in range(0,nL+1):
	for j in range(0,nH+1):	
		#print i*float(0.1/nL),j*float(0.1/nH)
		if(i%2==0):
			nodesIds.append( O.bodies.append(gridNode([i*float((xy+2*r)/nL),j*float((xy+2*r)/nH),zm],r,wire=False,fixed=False,material='gridNodeMat',color=color)) )
		else:
			if(i!=nL)and (j!=nH):
				x=(i)*float((xy+2*r)/nL)
				y=(j+0.5)*float((xy+2*r)/nH)
				nodesIds1.append( O.bodies.append(gridNode([x,y,zm],r,wire=False,fixed=False,material='gridNodeMat',color=color)) )
			
for i in range(0,(nL+1)/2+1,1):
	for j in range(0,nH,1):	
		cylIds.append(O.bodies.append( gridConnection(nodesIds[i*(nH+1)+j],nodesIds[i*(nH+1)+1+j],r,color=color,mask=5,material='gridCoMat') ))
		
for i in range(0,(nL+1)/2,1):
	for j in range(0,nH+1,1):	
		cylIds.append(O.bodies.append( gridConnection(nodesIds[i*(nH+1)+j],nodesIds[(i+1)*(nH+1)+j],r,color=color,mask=5,material='gridCoMat') ))
for i in range(0,(nL+1)/2,1):
	for j in range(0,nH,1):	
		cylIds.append(O.bodies.append( gridConnection(nodesIds[j+i*(nH+1)],nodesIds1[j+i*(nH)],r,color=color,mask=5,material='gridCoMat') ))
		cylIds.append(O.bodies.append( gridConnection(nodesIds[j+1+i*(nH+1)],nodesIds1[j+i*(nH)],r,color=color,mask=5,material='gridCoMat') ))
		cylIds.append(O.bodies.append( gridConnection(nodesIds[j+1+(i+1)*(nH+1)],nodesIds1[j+i*(nH)],r,color=color,mask=5,material='gridCoMat') ))
		cylIds.append(O.bodies.append( gridConnection(nodesIds[j+(i+1)*(nH+1)],nodesIds1[j+i*(nH)],r,color=color,mask=5,material='gridCoMat') ))
		
for i in range(0,(nL+1)/2,1):
	for j in range(0,nH,1):	
		O.bodies.append(pfacet(nodesIds[j+i*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+(i+1)*(nH+1)],wire=False,color=color,mask=5,material='gridCoMat') )
		O.bodies.append(pfacet(nodesIds[j+1+i*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+i*(nH+1)],wire=False,color=color,mask=5,material='gridCoMat') )
		O.bodies.append(pfacet(nodesIds[j+1+(i+1)*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+1+i*(nH+1)],wire=False,color=color,mask=5,material='gridCoMat') )
		O.bodies.append(pfacet(nodesIds[j+(i+1)*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+1+(i+1)*(nH+1)],wire=False,color=color,mask=5,material='gridCoMat') )	




O.engines=O.engines+[NewtonIntegrator(damping=0.8,gravity=[0,0,-9.81],label='newton'),PyRunner(command='history()',iterPeriod=1000)]
					 #, qt.SnapshotEngine(iterPeriod=1000,fileBase=O.tags['description']+'-i_'+str(O.iter)+'-',label='snapshooter')]
#####################



hminx0 = hMin(0)
hmaxx0 = hMax(0)

hminy0 = hMin(1)
hmaxy0 = hMax(1)

hminz0 = hMin(2)
hmaxz0 = hMax(2)
r_s=2*r
xy=hmaxx0

zmin=-(r+2*r_s)
zmax=r+2*r_s
zmean=(zmax+zmin)*0.5

O.materials.append(FrictMat(young=5e6,poisson=0.3,density=2650,frictionAngle=radians(20),label='sphereMat'))

walls=utils.aabbWalls((Vector3(-0.001,-0.001,zmin),Vector3(hmaxx0,hmaxx0,zmax)),thickness=2.*r,oversizeFactor=1.,material='sphereMat')
wallIds=O.bodies.append(walls)
O.bodies.erase(walls[5].id)

startxyz=Vector3(-0.001,-0.001,0)

color=Vector3(0.9,0.9,0.9)
kw={'color':[0.8,0.8,0.8],'wire':False,'dynamic':True,'material':'sphereMat'}


predicate=inAlignedBox(Vector3(-0.001,-0.001,zmax),Vector3(xy+r_s,xy+r_s,zmax+2*r))   
sp=SpherePack()
sp=pack.regularHexa(predicate, radius=r_s, gap=0.001*0.10,**kw)
O.bodies.append(sp)


#predicate=inAlignedBox(Vector3(-0.001,-0.001,zmin),Vector3(xy+r_s,xy+r_s,zmean-2*r))   
#sp=SpherePack()
#sp=pack.regularHexa(predicate, radius=r_s, gap=0.001*0.10,**kw)
#O.bodies.append(sp)

Gl1_Sphere.stripes=True


topPlate=utils.wall(position=hMaxS(2)+r_s,sense=0, axis=2,color=Vector3(1,0,0),material='sphereMat')
O.bodies.append(topPlate)


top_boundary=[]
top_boundary_m1=[]
bottom_boundary=[]
	
O.dt=0.2*PWaveTimeStep()


for i in range(10,221,21):			
		top_boundary.append(i)
for i in range(0,221,21):			
		bottom_boundary.append(i)
	

		
def moveGrid(v,d):
	newton.damping=d
	for i in bottom_boundary:
		O.bodies[i].state.blockedDOFs='xyzXYZ'
		O.bodies[i].state.vel=Vector3(0,-v,0.)

O.engines=O.engines+[PyRunner(command='pulloutover()',iterPeriod=100),PyRunner(command='Apply_load()',iterPeriod=100)]
#################
## APPLY LOAD  ##
#################
top=False

normalVEL=1e-03
normalSTRESS=80e03
sigmaN=0
Fn=0
S0=xy

def Apply_load():
	global sigmaN, Fn, top, load
	Fn=abs(O.forces.f(topPlate.id)[1])
	 
	sigmaN=Fn/S0 #constant normal load
		
	if abs(topPlate.state.vel[1])<abs(normalVEL):
		topPlate.state.vel[1]+=normalVEL/1000
		
	if (sigmaN>(0.9*normalSTRESS))and (top==False):
		topPlate.state.vel[1]=normalVEL*((normalSTRESS-sigmaN)/normalSTRESS)
		
	if abs((normalSTRESS-sigmaN)/normalSTRESS)<0.001:
		top=True
		topPlate.state.vel[1]=0
		
		#if(unbalancedForce()<0.01):
			#O.pause()
			#O.tags['description']='packing_Hexa_CuiOSul_'+str(normalSTRESS)
			#saveData()
			#pos = [] #Cette liste contiendra ma map en 2D
			#for i in range(n_sp):
				#pos.append([0] * 4) 
			#pos=spheres_pos(pos)
			#save(pos)
			#load=load+1
			#top=False
			#O.run()

 #plot some results


def history():
	global x,y,z,ytop
	F=Vector3(0,0,0)
	for i in bottom_boundary:
		F+=O.forces.f(i)
	#print y
	Fnorm=F.norm()

	y=O.bodies[bottom_boundary[0]].state.pos[1]
	xtop=O.bodies[top_boundary[0]].state.pos[0]
	ytop=O.bodies[top_boundary[0]].state.pos[1]
	ztop=O.bodies[top_boundary[0]].state.pos[2]
	Fx=F[0]
	Fy=F[1]
	Fz=F[2]
	#print Fnorm
	plot.addData(y=y,xtop=xtop,ytop=ytop,ztop=ztop,Fnorm=Fnorm,Fx=Fx,Fy=Fy,Fz=Fz,t=O.time,i=O.iter,u=unbalancedForce())

def move():
	for i in nodesIds:
		O.bodies[i].state.blockedDOFs='xyXYZ'
	
from math import *
from yade import plot


test1=0
test2=0
def pulloutover():
	global test1, test2
	if((-0.05<y<-0.049)and (test1==0)):
		test1=1
		O.pause()
		O.save(O.tags['description']+str(y)+'.xml.bz2',quiet=True)
		O.run()
	if((-0.051<y<-0.05)and (test2==0)):
		test2=2
		O.pause()
		O.save(O.tags['description']+str(y)+'.xml.bz2',quiet=True)
		O.run()
	if(y<-0.001):
		O.pause()
		saveData()
		O.save(O.tags['description']+'_final.xml.bz2',quiet=True)
		O.exitNoBacktrace()
	#print 'ok'

def saveData():
  plot.saveDataTxt('data_pullout/'+O.tags['description']+'.dat',vars=('t','i','xtop','ytop','ztop','y','Fnorm','Fx','Fy','Fz','u'))

O.saveTmp()

#O.run(40000)
#O.save(O.tags['description']+'_init.xml.bz2',quiet=True)
#moveGrid(0.05,0.8)
#O.run()

plot.plots={'y':('Fnorm','Fy'),'t':('u')}
plot.plot(subPlots=True)

#  VISUALIZATION  

qt.Controller()
qtv = qt.View()
qtr = qt.Renderer()
qtr.light2=True        # back
qtr.lightPos=Vector3(1200,1500,500)
qtr.bgColor=[1,1,1]

## view from side (for full model)
#qtv.screenSize=Vector2i(1920,1145)

#qtv.upVector= Vector3(-0.10856143151406383,0.06871364753099884,0.9917120803088985)


#qtv.lookAt=Vector3(-0.505406128227311,0.42266320021406945,-0.13413842322923802)

#qtv.viewDir=Vector3(-0.8277689260283403,0.5461659825561284,-0.1284574816831906)

#qtv.eyePosition= Vector3(0.3223627978010293,-0.1235027823420589,-0.005680941546047429)

qtv.timeDisp=''


