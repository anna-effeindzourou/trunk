# -*- coding: utf-8 -*-
from yade import ymport, plot, qt
import numpy

#import yade.log
#yade.log.setLevel('WireMat',yade.log.TRACE)
#yade.log.setLevel('Law2_ScGeom_WirePhys_WirePM',yade.log.TRACE)
#yade.log.setLevel('Ip2_WireMat_WireMat_WirePhys',yade.log.TRACE)



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
	#print nx,ny
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
				net+=[sphere((x,y,z),radius=radius,**kw)]
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
				net+=[sphere((x,y,z),radius=radius,**kw)]
			# set values for next section
			xstart = xstart - 0.5*mosx*pow(-1,i)
			nx = int(nx + 1*pow(-1,i))
	#netId=O.bodies.append(net)

			
	#print(nx,ny)
	#print(lx,ly)

	return [net,lx,ly]


#### define parameters for the net
# wire diameter
d = 3.0/1000.
# particle radius
radius = d*4.
# define piecewise lineare stress-strain curve [Pa]
#strainStressValues=[(0.01,1.50E+09),(0.015,1.85E+09),(0.02,1.95E+09),(0.027,1.95E+09),(0.035,1.50E+09)]
strainStressValues=[(0.0125,1.67E+09),(0.02,1.77E+09),(0.027,1.77E+09),(0.035,1.32E+09)]
# mesh geometry
mosx = 0.083
mosy = 0.143
# net dimension
cornerCoord=[0,0,0]
Lx = 1.079
Ly = 1.001

# elastic material properties
particleVolume = 4./3.*pow(radius,3)*pi
#particleMass = 1.65/175.
particleMass = 1.65*50

density = particleMass/particleVolume
young = strainStressValues[0][1] / strainStressValues[0][0]
poisson = 1
fricangle=30

#### material definition
lambdaF= 0.05
lambdau= 0.35
axis=1
print 'lambdaF: ',lambdaF,' ,lambdau: ',lambdau

meshMatInternal = O.materials.append( WireMat( young=young,poisson=poisson,frictionAngle=radians(fricangle),density=density,isDoubleTwist=False,diameter=d,strainStressValues=strainStressValues,lambdaF=lambdaF, lambdau=lambdau,seed=12,type=2) )

O.tags['description']='t21_tensile_test_axis_'+str(axis)+'_lambdaF_'+str(lambdaF)+'_lambdau_'+str(lambdau)+'_'

#### properties of particles
kw = {'color':[1,1,0],'wire':True,'highlight':False,'fixed':False,'material':meshMatInternal}


##### create packing
[netpack,lx,ly] = diamondNet( radius=radius, cornerCoord=cornerCoord, xLength=Lx, yLength=Ly, mosx=mosx, mosy=mosy, startAtCorner=False, **kw )
netpackId=O.bodies.append(netpack)



#### define engines to create link
interactionRadius=3.45
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=interactionRadius,label='aabb')]), 
	InteractionLoop(
	[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=interactionRadius,label='Ig2ssGeom')],
	[Ip2_WireMat_WireMat_WirePhys(linkThresholdIteration=1,label='wire_wire')],
	[Law2_ScGeom_WirePhys_WirePM(linkThresholdIteration=1,label='Law_1')]
	),
		GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.5,label='ts'), 

	NewtonIntegrator(damping=0.),
	
]


#### time step definition for first time step to create links
O.step()



#set the boundary conditions
leftNodes = []
rightNodes = []
posIds = []
negIds = []

for i in netpackId:
	if(O.bodies[i].state.pos[0]==0):
		O.bodies[i].shape.color=Vector3(0,0,1)
		leftNodes.append(i)
		O.bodies[i].state.blockedDOFs='xyz'
	if(O.bodies[i].state.pos[0]==lx):
		rightNodes.append(i)
		O.bodies[i].shape.color=Vector3(0,0,1)
		O.bodies[i].state.blockedDOFs='xyz'
	if(O.bodies[i].state.pos[1]==0):
		O.bodies[i].shape.color=Vector3(0,1,0)
		negIds.append(i)
		O.bodies[i].state.blockedDOFs='yz'
	if(O.bodies[i].state.pos[1]==ly):
		O.bodies[i].shape.color=Vector3(1,1,1)
		posIds.append(i)
		O.bodies[i].state.blockedDOFs='yz'

#set the tension velocity
#vel=0.005

vel=17e-05

L=0
sideIds=leftNodes+rightNodes
topbotIds=posIds+negIds
axis=1
if(axis==1):
	for i in posIds:
		O.bodies[i].state.vel=vel*Vector3(0,1,0)
	L=Lx
	for i in sideIds:
		O.bodies[i].state.blockedDOFs='xz'
		
	O.bodies[posIds[7]].state.blockedDOFs='xyz'
	O.bodies[negIds[7]].state.blockedDOFs='xyz'
else:
	for i in rightNodes:
		O.bodies[i].state.vel=vel*Vector3(1,0,0)
	L=Ly
	for i in topbotIds:
		O.bodies[i].state.blockedDOFs='yz'


################################
#####   Plot some results   ####
################################

plot.plots={'un':('Fn','Fnn')}
plot.plot(noShow=False, subPlots=False)

def addPlotData():
	Fn = 0.
	Fnn = 0.
	if(axis==1):
		for i in posIds:
			Fn +=	O.forces.f(i)[axis]
			
		Fn=abs(Fn)	
		sigma=Fn/(1000*L)
		
		un = (O.bodies[posIds[0]].state.pos[axis] - O.bodies[posIds[0]].state.refPos[axis])
		eps=un/L
	
	else:
		for i in rightNodes:
			Fn += O.forces.f(i)[axis]
		Fn=abs(Fn)	
		sigma=Fn/(1000*L)
		
		un = (O.bodies[rightNodes[0]].state.pos[axis] - O.bodies[rightNodes[0]].state.refPos[axis])
		eps=un/L
	if un*1000 > 300:
		O.pause()
	#print eps	
	plot.addData( eps=eps*1000, un=un*1000, Fn=Fn,sigma=sigma )

O.engines+=[PyRunner(initRun=True,iterPeriod=1000,command='addPlotData()'),
	 VTKRecorder(iterPeriod=1000,initRun=True,fileName='paraview/'+O.tags['description']+'-',recorders=['spheres','velocity','intr'])]
def saveData():
	plot.saveDataTxt('data/'+O.tags['description']+'.dat',vars=('un','eps','Fn','sigma'))
	
#### define engines for simulation
v = qt.Controller()
v = qt.View()
rr = qt.Renderer()
rr.intrAllWire = True


#### initializes now the interaction detection factor
aabb.aabbEnlargeFactor=-1.
Ig2ssGeom.interactionDetectionFactor=-1.

