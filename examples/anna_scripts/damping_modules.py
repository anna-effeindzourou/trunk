# -*- coding: utf-8
from yade import ymport, utils,pack,export,qt
import gts,os
from yade import geom
#import matplotlib
from yade import plot
#from pylab import *
#import os.path, locale



#################################
#####     FUNCTIONS          ####
#################################
def hMax(n):
    idHMax=0
    hMax=-1000000.0
    for i in O.bodies:
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
  
def limitfinder():
	for b in O.bodies:
		if(b.state.pos[2]>=L-2*radius):
			if isinstance(b.shape,GridNode):
				top_boundary.append(b.id)
				b.shape.color=(1,0,0)
				b.state.blockedDOFs='xyzXYZ'
		if(b.state.pos[2]<0.1*radius ):
			if isinstance(b.shape,GridNode):
				bottom_boundary.append(b.id)
				b.state.blockedDOFs='xyzXYZ'
				b.shape.color=(1,0,0)	
				
##############################
#####     SCRIPT          ####
##############################

      
############################
### DEFINING PARAMETERS  ###
############################

 #GEOMETRIC :dimension of the rectangular box
a=1.5	 #lenght 
h=1.5	 # depth  
d=0.4	 # hight 



#MATERIAL PROPERTIES
frictionAngle=radians(17)
Kn=170e6 #normal stiffness Plassiard thesis
Ks=35e6 #tangential stiffness Plassiard thesis


#PARTICLE SIZE

Rmax=0.0507/float(2)
rmin=0.025/float(2)
Rf=(Rmax-rmin)/float(Rmax+rmin) # dispersion (Rs±Rf*Rs)
Rs=rmin/float(1-Rf)# mean particle radius
nSpheres=2000# number of particles

#ENGINES PARAMETERS
g=-9.81
damping_coeff = 0.0

#gravelPhys
krot=1.8
eta=1.8
alpha=1
c=20e3
#MATERIAL'S PARAMETERS
#gravel
poisson_g=Ks/Kn
young_g=Kn
density=1650
##PSD PARAMETER (STEP)
n=10
##

#TIME STEP
#dt=1e-5
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
		Ig2_Wall_Sphere_ScGeom(),Ig2_Wall_PFacet_ScGeom(),
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





#######################
####  PROPERTIES    ###
#######################

Ey=12.5e8
poisson=1
density=2.60e4
Kr=6250
frictionAngle=radians(00)
radius=.01
#sigma=-12e6

##### Parameters of a rectangular grid ###
L=0.6 #length [m]
l=0.5	#half width	(radius) [m]
nbL=12#number of nodes for the length	[#] doit etre paire
nbl=18	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!

##nbL=1 #number of nodes for the length	[#] doit etre paire
##nbl=4	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!

r=radius
color=[155./255.,155./255.,100./255.]

oriBody = Quaternion(Vector3(0,0,1),(pi/2))
nodesIds=[]
nodesIds1=[]
cylIds=[]
pfIds=[]


#top_boundary=[]
#bottom_boundary=[]

#####################
####  MATERIAL    ###
#####################
O.materials.append(CohFrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=0,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=True,alphaKr=0,label='NodeMat'))
O.materials.append(FrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=0,label='Pmat'))

#O.materials.append(FrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=26.57,label='Smat'))
###############################
####  SAMPLE GENERATION     ###
###############################

##pile=ymport.text('spheres.txt')
##pile2=O.bodies.append(pile)
##sup()

#################################
####  MEMBRANE GENERATION     ###
#################################
mesh=2
def cylinder(X,Y,Z,nodesIds=[],pfIds=[]):
	#Create all nodes first :
	for i in range(0,nbL+1):
		for j in range(0,nbl):
			z=i*L/float(nbL)-Z+r
			y=l*sin(2*pi*j/float(nbl))+Y
			x=l*cos(2*pi*j/float(nbl))+X
			nodesIds.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )

	##Create connection between the nodes
	for i in range(0,nbL+1):
		for j in range(0,nbl-1):
			O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[i*nbl+j+1],r,color=color) )
	for i in range(0,nbL,1):	
			for j in range(0,nbl):
				O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],r,color=color) )
	for i in range(-1,nbL):
		j=nbl
		O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color) )
		
	for i in range(0,nbL):
		for j in range(0,nbl-1):
				if (j%2==0):
					O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j+1],r,color=color) )
				else:
					O.bodies.append( gridConnection(nodesIds[(i+1)*nbl+j],nodesIds[i*nbl+j+1],r,color=color) )
	for i in range(0,nbL):
		j=nbl
		#O.bodies[nodesIds[(i-1)*nbl+j]].shape.color=Vector3(155./255.,155./255.,1.)
		O.bodies[nodesIds[(i)*nbl+j-1]].shape.color=Vector3(1,0,0)
		O.bodies.append( gridConnection(nodesIds[(i-1)*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color) )		


	###Create PFacets
	##wire=True
	for i in range(0,nbL):
		for j in range(0,nbl-1):
				if (j%2==0):
					pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],color=color,mask=5)))
					pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j+1],nodesIds[(i)*nbl+j+1],color=color,mask=5)))
				else:
					pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i)*nbl+j+1],color=color,mask=5)))
					pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j+1],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],color=color,mask=5)))
		
	for i in range(0,nbL,1):
		j=nbl
		pfIds.append(O.bodies.append(pfacet( nodesIds[i*nbl+j],nodesIds[(i-1)*nbl+j],nodesIds[(i+1)*nbl+j-1],color=color )))
		pfIds.append(O.bodies.append(pfacet( nodesIds[(i)*nbl+j-1],nodesIds[(i+1)*nbl+j-1],nodesIds[(i-1)*nbl+j],color=color )))
		
mod=0.9
zdamp=0.39
cylinder(0.5,-mod,zdamp,nodesIds=[],pfIds=[])
cylinder(-0.5,-mod,zdamp,nodesIds=[],pfIds=[])

cylinder(0,0,zdamp,nodesIds=[],pfIds=[])

cylinder(0.5,mod,zdamp,nodesIds=[],pfIds=[])
cylinder(-0.5,mod,zdamp,nodesIds=[],pfIds=[])

cylinder(1,0,zdamp,nodesIds=[],pfIds=[])
cylinder(-1,0,zdamp,nodesIds=[],pfIds=[])
#cylinder(-1,1,0.4,nodesIds=[],pfIds=[])



############################
###### ENGINE DEFINITION  ##
############################  	
O.dt=0.5*PWaveTimeStep()
O.engines=O.engines+[
	NewtonIntegrator(damping=0.2,gravity=(0,0,g),label='Newton'),
	]	