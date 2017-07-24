# -*- coding: utf-8 -*-
# encoding: utf-8
from yade import pack,export,qt
from yade.gridpfacet import *
from yade import utils,plot
## definition of some colors for colored text output in terminal
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
BLACK = '\033[0m'

#### short description of script
print BLUE+'This script generates the particles for the net based on a chained link'+BLACK
print BLUE+'and stores them in a file named '+RED+'net-pack.spheres'+BLUE+'.'+BLACK

#### definition for net generation
def diamondNet( radius, cornerCoord=[0,0,0], xLength=1., yLength=0.5, mosx=0.083, mosy=0.143, startAtCorner=False, **kw ):
	"""Definition of the particles for a chained link wire net in the x-y-plane for the WireMatPM.

	:param radius: radius of the particle
	:param cornerCoord: coordinates of the lower left corner of the net
	:param xLenght: net length in x-direction
	:param yLenght: net length in y-direction
	:param mosx: mesh opening size in x-direction
	:param mosy: mesh opening size in y-direction
	:param startAtCorner: if true the generation starts with a double-twist at the lower left corner
	:param isSymmetric: defines if the net is symmetric with respect to the y-axis

	:return: set of spheres which defines the net (net) and exact dimensions of the net (lx,ly).
	
	note::
	This packing works for the WireMatPM only. The particles at the corner are always generated first. For examples on how to use this packing see examples/WireMatPM. In order to create the proper interactions for the net the interaction radius has to be adapted in the simulation.

	"""
	# check input dimension
	if(xLength<mosx): raise ValueError("xLength must be greater than x!");
	if(yLength<mosy/2.): raise ValueError("yLength must be greater than y/2!");
	xstart = cornerCoord[0]
	ystart = cornerCoord[1]
	z = cornerCoord[2]

	# number of particles in x-direction and real length lx
	nx = int( xLength/mosx ) + 2
	lx = (nx-1)*mosx
	# number of particles in y-direction and real length ly
	ny = int( yLength/(0.5*mosy) ) + 1
	print nx,ny
	ly = (ny-1)*0.5*mosy
	net = []
	netcyl = []
	netcylId = []
	# generate corner particles
	if startAtCorner:
		jump=1
		for i in range(ny):
			y = ystart + i*mosy/2.
			for j in range(nx):
				x = xstart + j*mosx
				net+=[gridNode((x,y,z),radius=radius,**kw)]
			# set values for next section
			xstart = xstart - 0.5*mosx*pow(-1,i+jump)
			nx = int(nx + 1*pow(-1,i+jump))
	else:
		nx-=1
		xstart+=0.5*mosx 
		for i in range(ny):
			y = ystart + i*mosy/2.
			for j in range(nx):
				x = xstart + j*mosx
				net+=[gridNode((x,y,z),radius=radius,**kw)]
			# set values for next section
			xstart = xstart - 0.5*mosx*pow(-1,i)
			nx = int(nx + 1*pow(-1,i))
	netId=O.bodies.append(net)
	
	#print len(netId),nx
	#for i in range(0,len(netId)-nx):
	cx=nx-1
	#if (nx%2==0): cx=nx-1
	#else: cx=nx
	
	print 'nx: ',nx,', ny: ',ny
	#nx=nx+1
	#cx=cx+1
	for j in range(0,ny-1):#ny-1
		#j=4
		k=j/2
		k=k-j
		#print k,j
		for i in range(k,cx+k):
			netcyl+=[gridConnection(i+j*(ny-1),i+nx+j*(ny-1),radius=radius,color=Vector3(1,0,0),material=meshMat)]
		if(j<=1):
			for i in range(0,cx): #j=0 0,cx; j=1 1, cx+1;j=2 1, cx+1;
				netcyl+=[gridConnection(i+j*(ny-1),i+(nx-1)+j*(ny-1),radius=radius,color=Vector3(1,0,0),material=meshMat)]
		else:
			#print k,j
			if(j%2==0): k=k
			else: k=k+1
			print k,j
			
			for i in range(k,cx+k):
				netcyl+=[gridConnection(i+j*(ny-1),i+(nx-1)+j*(ny-1),radius=radius,color=Vector3(1,0,0))]
			
	netcylId=O.bodies.append(netcyl)	
	return [netId,netcylId,lx,ly]







###################
#####   Box   #####
###################
## material properties
density = 6e4
young = 1e8
poisson = 1
fricangle = 10
cohesion = 1e19
multYoungInternal = 1	## multiplicator for internal Young, i.e., moment resistance of connections

meshMatInternal = O.materials.append( CohFrictMat( young=multYoungInternal*young,poisson=poisson,density=density,frictionAngle=radians(fricangle),normalCohesion=cohesion,shearCohesion=cohesion,momentRotationLaw=True,alphaKr=0, label='boxMatInternal' ) )

meshMat = O.materials.append( FrictMat( young=young,poisson=poisson,frictionAngle=radians(fricangle),density=density, label='boxMat' ) )
strainStressValues=[(0.0019230769,2.5e8),(0.0192,3.2195e8),(0.05,3.8292e8),(0.15,5.1219e8),(0.25,5.5854e8),(0.3,5.6585e8),(0.35,5.6585e8)]
d=3./1000
radius = d*4.

particleVolume = 4./3.*pow(radius,3)*pi
particleMass = 3.9/1000.
density = particleMass/particleVolume
young = strainStressValues[0][1] / strainStressValues[0][0]

meshMatInternal = O.materials.append( WireMat( young=young,poisson=poisson,frictionAngle=radians(30),density=density,isDoubleTwist=False,diameter=d,strainStressValues=strainStressValues ) )


##############################
#####   Common Engines   #####
##############################
O.engines=[
	ForceResetter(),
	InsertionSortCollider(
		[

			],
		sortThenCollide=False
	),
	InteractionLoop(
		[Ig2_GridNode_GridNode_GridNodeGeom6D(),
		],
		[Ip2_WireMat_WireMat_WirePhys(linkThresholdIteration=1,label='interactionPhys'),
		],
		[Law2_ScGeom_WirePhys_WirePM(linkThresholdIteration=1,label='interactionLaw'),
		]
	),
	NewtonIntegrator(gravity=(0,0,0),damping=0.5,label='newton'),
]


#### define parameters for the net
# mesh geometry
mosx = 0.083
mosy = 0.143
# wire diameter
d = 3.0/1000.
# net dimension
cornerCoord=[0,0,0]
Lx = 1.079
Ly = 1.001


#### properties of particles
radius = d*4.
#kw = {'color':[1,1,0],'wire':True,'highlight':False,'fixed':False,'material':netMat}
kw = {'color':[1,1,0],'wire':True,'material': meshMatInternal,'highlight':False,'fixed':False}


##### create packing
[netpack,netcyl,lx,ly] = diamondNet( radius=radius, cornerCoord=cornerCoord, xLength=Lx, yLength=Ly, mosx=mosx, mosy=mosy, startAtCorner=False, **kw )
print lx,ly

sideNodes = []
topNodes = []
botNodes = []
for i in netpack:
	if(O.bodies[i].state.pos[0]==0)or(O.bodies[i].state.pos[0]==lx):
		O.bodies[i].shape.color=Vector3(1,0,0)
		sideNodes.append(i)
		O.bodies[i].state.blockedDOFs='xyz'
	if(O.bodies[i].state.pos[1]==0):
		O.bodies[i].shape.color=Vector3(0,1,0)
		botNodes.append(i)
		O.bodies[i].state.blockedDOFs='xyz'
	if(O.bodies[i].state.pos[1]==ly):
		O.bodies[i].shape.color=Vector3(1,1,1)
		topNodes.append(i)
		O.bodies[i].state.blockedDOFs='xyz'
for i in botNodes:
	O.bodies[i].state.vel=Vector3(0,-1,0)
		
#Fy=0
O.engines+=[PyRunner(iterPeriod=100,command='history()')]
		

plot.plots={'pos':'Fy'}
plot.plot(subPlots=False)
posp=O.bodies[botNodes[0]].state.pos[1]-O.bodies[topNodes[0]].state.pos[1]
def history():
	global posp
	Fy=0
	for i in range(0,len(botNodes)):
		Fy+=O.forces.f(i)[1]
		
	pos=O.bodies[botNodes[0]].state.pos[1]-O.bodies[topNodes[0]].state.pos[1]
	dpos=pos-posp
	posp=pos
	
  	plot.addData(pos=-pos,dpos=dpos,
		Fy=Fy,
		t=O.time)
plot.plot(subPlots=False)


O.dt=0.1*PWaveTimeStep()

## to see it
v=qt.Controller()
v=qt.View()

#export.text('net-pack.spheres')
