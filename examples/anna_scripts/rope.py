from yade import export,plot
from yade.gridpfacet import *


def hMinSpheres(n):
	idHMin=0
	hMin=100000.0
	for i in O.bodies:
		h=i.state.pos[n]
		if (type(i.shape)==Sphere):
			if (h<hMin):
				hMin=h
				idHMin=i.id
	return (hMin)    
	
def hMaxSpheres(n):
	idHMax=0
	hMax=-1000000.0
	for i in O.bodies:
		h=i.state.pos[n]
		if (type(i.shape)==Sphere):
			if (h>hMax):
				hMax=h
				idHMax=i.id
	return (hMax)	


try:
	os.mkdir('data')
except:
	pass
try:
	os.mkdir('paraview')
except:
	pass

betaR=5
damping=0.2
O.engines=[
	ForceResetter(),
	InsertionSortCollider(
		[
		],
		sortThenCollide=False
	),
	InteractionLoop(
		[
			Ig2_GridNode_GridNode_GridNodeGeom6D()
		],
		[
			Ip2_WireMat_WireMat_WirePhys(linkThresholdIteration=1,label='interactionPhys')
		],
		[
			Law2_ScGeom_WirePhys_WirePM(linkThresholdIteration=1,label='interactionLaw')
		]
	),
	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.5,label='ts'), 
	NewtonIntegrator(gravity=(0,0,-9.81),damping=damping,label='newton'),
	PyRunner(initRun=True,iterPeriod=1,command='plotData()')

]
O.engines[1].avoidSelfInteractionMask=4

## material properties
poisson = 0.2
fricangle = 27
d=3./1000
radius = d*2
particleVolume = 4./3.*pow(radius,3)*pi
particleperm2=202
meshmassperm2=1.680
particlemass= (meshmassperm2/particleperm2)
print particlemass
density = particlemass/particleVolume
print density
coeff_dech_cyl=11
strainStressValues=[(0.2,1.0E+1000),(120,1.0E+1000)] 

young = (strainStressValues[0][1] / strainStressValues[0][0])

#cylinderMat.young=250e6
ropeMat = O.materials.append( WireMat( young=young,poisson=poisson,frictionAngle=radians(fricangle),density=density,isDoubleTwist=False,diameter=8e-03,strainStressValues=strainStressValues) )
cylinderMat = O.materials.append(NormalInelasticMat(young=young ,poisson=poisson,frictionAngle=radians(0),density=density,coeff_dech=coeff_dech_cyl))
z=0.25
nx=41
rope=[]
fixed=True
R=0.5+30e-03
r=4e-03
dtheta=2.*pi/nx
for i in range(0,nx):
	theta=dtheta*i
	x=R*cos(theta)
	y=R*sin(theta)
	u=gridNode(center=(x,y,z),radius=r,fixed=fixed,material=ropeMat)
	rope.append(u)
color=Vector3(1,1,1)
ropeIds=O.bodies.append(rope)		

#for i in ropeIds:
	#O.bodies[i].state.blockedDOFs='zXYZ'

ropering=[]
for i in range (0,len(ropeIds)-1):
	ropering.append(O.bodies.append(gridConnection(ropeIds[i],ropeIds[i+1],radius=4e-03,color=color,material=cylinderMat,mask=1)))
ropering.append(O.bodies.append(gridConnection(ropeIds[0],ropeIds[len(ropeIds)-1],radius=4e-03,color=color,material=cylinderMat,mask=1)))


def connected():
	n=0
	for i in ropering:
		id1=O.bodies[i].shape.node1.id
		id2=O.bodies[i].shape.node2.id
		n+=O.interactions[id1,id2].isReal
	return n

##################################
#### PLOT + SIMULATION STORAGE ###
##################################
plot.plots={'t':('UnbF',None,('hMax','r--'),None,('vz','b--'))}
plot.plot(subPlots=False)

def plotData():
	print connected()

Gl1_PFacet.wire=True
#from yade import qt,plot
			
#qr=qt.Renderer()
##qr.wire=True
##### to see it
#v=qt.Controller()

os=int(0.5/O.dt)
print 'os: ',os
#O.run(os)
#impact()

#O.run(3*os,True)