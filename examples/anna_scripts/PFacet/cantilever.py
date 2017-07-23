# enoding: utf-8
from yade import utils,pack,geom,qt,export,ymport

try:
	os.mkdir('data')
except:
	pass
try:
	os.mkdir('paraview')
except:
	pass


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
def scalar(u,v):
	ps=u[0]*v[0]+u[1]*v[1]+u[2]*v[2]
	return ps

def cross(u,v):
	ps=Vector3(u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2] ,u[0]*v[1]-u[1]*v[0])
	return ps   
 
O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_Sphere_Aabb(),
		Bo1_Box_Aabb(),
		Bo1_PFacet_Aabb(),
	]),
	InteractionLoop(
		[
		Ig2_GridNode_GridNode_GridNodeGeom6D(),
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_Sphere_PFacet_ScGridCoGeom(),
		Ig2_Box_Sphere_ScGeom()
		],
		[
		Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=False,label='ipf'),
		Ip2_FrictMat_FrictMat_FrictPhys()
		],
		[
		Law2_ScGeom_FrictPhys_CundallStrack(),
		Law2_ScGeom6D_CohFrictPhys_CohesionMoment(membrane=True) ,
		Law2_ScGridCoGeom_FrictPhys_CundallStrack()
		]
	),
]


#O.engines[1].avoidSelfInteractionMask=4
angle=20

O.materials.append(CohFrictMat(young=5e6,poisson=0.3,density=2650,frictionAngle=radians(angle),normalCohesion=3e100,shearCohesion=3e100,momentRotationLaw=False,label='gridNodeMat'))
O.materials.append(FrictMat(young=5e6,poisson=0.3,density=2650,frictionAngle=radians(angle),label='gridCoMat'))

### Parameters of a rectangular grid ###

d=2

r=0.0005

L=0.05
H=0.1
nH=10
nh=nH-1
nL=8
nS=nL
if(nL%2==0): nS=nL-1
z=0
y=0
x=0

color=[255./255.,102./255.,0./255.]
zm=-0.046

boundary_cond=1
displacement=0

O.tags['description']='tension_boundary_cond_'+str(boundary_cond)+'_displacement_'+str(displacement)
print O.tags['description']

#STEP 1: GENERATION OF THE MEMBRANE
nodesIds=[]
nodesIds1=[]
cylIds=[]
PFacetIds=[]
	
for i in range(0,nL+1):
	for j in range(0,nH+1):	
		#print i*float(0.1/nL),j*float(0.1/nH)
		if(i%2==0):
			nodesIds.append( O.bodies.append(gridNode([i*float((L+2*r)/nL),j*float((H+2*r)/nH),zm],r,wire=False,fixed=False,material='gridNodeMat',color=color)) )
		else:
			if(i!=nL)and (j!=nH):
				
				nodesIds1.append( O.bodies.append(gridNode([i*float((L+2*r)/nL),(j+.5)*float((H+2*r)/nH),zm],r,wire=False,fixed=False,material='gridNodeMat',color=color)) )
		
for i in range(0,(nL+1)/2,1):
	for j in range(0,nH,1):	
		pfacetCreator3(nodesIds[j+i*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+(i+1)*(nH+1)],pfIds=PFacetIds,wire=False,mask=5,material='gridCoMat',color=color,fixed=False)
		pfacetCreator3(nodesIds[j+1+i*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+i*(nH+1)],pfIds=PFacetIds,wire=False,mask=5,material='gridCoMat',color=color,fixed=False)
		pfacetCreator3(nodesIds[j+1+(i+1)*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+1+i*(nH+1)],pfIds=PFacetIds,wire=False,mask=5,material='gridCoMat',color=color,fixed=False)
		pfacetCreator3(nodesIds[j+(i+1)*(nH+1)],nodesIds1[j+i*(nH)],nodesIds[j+1+(i+1)*(nH+1)],pfIds=PFacetIds,wire=False,mask=5,material='gridCoMat',color=color,fixed=False)		



top_boundary=[]
bottom_boundary=[]
	
#top_boundary=[10, 31, 52, 73, 94, 115]

n=nodesIds[len(nodesIds)-1]
for i in range(nH,n+21,21):			
		top_boundary.append(i)
		O.bodies[i].shape.color=Vector3(0,0,1)
		O.bodies[i].state.blockedDOFs='xyz'

for i in range(0,n,21):			
		bottom_boundary.append(i)
		O.bodies[i].shape.color=Vector3(1,0,0)
		O.bodies[i].state.blockedDOFs='xyz'


color=Vector3(0.1,0.1,0.1)


p1incr=0
p2incr=0
p3incr=0
sigma=-3e5
n1=0
n2=0
n3=0
step=1e-6
def Apply_confiningpressure():
	global p1incr,p2incr,p3incr,p1n,p2n,p3n
	for i in PFacetIds: 
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
		#print 'Fapplied = ',Fapplied,' ,normal = ',normal
		
		
		F=Fapplied
		#print 'F=',F
		p1n=F*p1
		p2n=F*p2
		p3n=F*p3
		

		
		
		if(abs(p1incr)<abs(p1n)):
			p1incr+=p1n*step
			p2incr+=p2n*step
			p3incr+=p3n*step
		#print 'p1incr=',p1incr,', p1n=',p1n
		
		p1normal=p1incr*normal
		p2normal=p2incr*normal
		p3normal=p3incr*normal
		
		O.forces.addF(O.bodies[i].shape.node1.id,p1normal,1)
		O.forces.addF(O.bodies[i].shape.node2.id,p2normal,1)
		O.forces.addF(O.bodies[i].shape.node3.id,p3normal,1)



O.dt=0.2*PWaveTimeStep()

 #plot some results



def history():
    global y0,y1,y2,y3,y4,Fnorm
    if(displacement==0):
	F=Vector3(0,0,0)
	for i in bottom_boundary:
	    F+=O.forces.f(i)
  
    Fnorm=F.norm()

    y0=O.bodies[bottom_boundary[0]].state.pos[1]
    y1=O.bodies[bottom_boundary[1]].state.pos[1]
    y2=O.bodies[bottom_boundary[2]].state.pos[1]
    y3=O.bodies[bottom_boundary[3]].state.pos[1]
    y4=O.bodies[bottom_boundary[4]].state.pos[1]
    
    Fx=F[0]
    Fy=F[1]
    Fz=F[2]
    #print Fnorm
    plot.addData(y0=y0,y1=y1,y2=y2,y3=y3,y4=y4,Fnorm=Fnorm,Fx=Fx,Fy=Fy,Fz=Fz)


   

	
from math import *
from yade import plot


        
def saveData():
  plot.saveDataTxt('data/'+O.tags['description']+'.dat',vars=('y0','y1','y2','y3','y4','Fnorm','Fx','Fy','Fz'))

O.saveTmp()




plot.plots={'y0':('Fnorm','Fy'),'y1':('Fnorm','Fy'),'y2':('Fnorm','Fy'),'y3':('Fnorm','Fy'),'y4':('Fnorm','Fy')}
plot.plot(subPlots=True)


O.engines=O.engines+[PyRunner(iterPeriod=50000,dead=False,command='Apply_confiningpressure()',label='Applycpressure'),NewtonIntegrator(gravity=(0,0,0),damping=0.1,label='newton'),PyRunner(initRun=True,iterPeriod=1,command='history()'),VTKRecorder(iterPeriod=500,dead=False,initRun=True,fileName='paraview/'+O.tags['description']+'_',recorders=['spheres','velocity','intr'],label='VTK')
]

#  VISUALIZATION  

qt.Controller()
qtv = qt.View()
qtr = qt.Renderer()
qtr.light2=True        # back
qtr.lightPos=Vector3(1200,1500,500)
qtr.bgColor=[1,1,1]

Gl1_PFacet.wire=True
