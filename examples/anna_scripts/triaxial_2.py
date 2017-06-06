# -*- coding: utf-8
from yade import ymport, utils,pack,export,qt
import gts,os
from yade import geom
#import matplotlib
from yade import plot
#from pylab import *
#import os.path, locale
			
##############################
#####     SCRIPT          ####
##############################

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

Ey=12.5e8
poisson=1
density=2.60e4
Kr=6250
frictionAngle=radians(00)
radius=.00063
sigma=-12e6

#### Parameters of a rectangular grid ###
L=0.203 #length [m]
l=0.101/2.	#half width	(radius) [m]
nbL=24#number of nodes for the length	[#] doit etre paire
nbl=36	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!

#nbL=1 #number of nodes for the length	[#] doit etre paire
#nbl=4	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!

r=radius
color=[155./255.,155./255.,100./255.]

oriBody = Quaternion(Vector3(0,0,1),(pi/2))
nodesIds=[]
nodesIds1=[]
cylIds=[]
pfIds=[]


top_boundary=[]
bottom_boundary=[]

####################
###  MATERIAL    ###
####################
O.materials.append(CohFrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=0,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=True,alphaKr=0,label='NodeMat'))
O.materials.append(FrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=0,label='Pmat'))

O.materials.append(FrictMat(young=Ey,poisson=poisson,density=density,frictionAngle=26.57,label='Smat'))
##############################
###  SAMPLE GENERATION     ###
##############################

pile=ymport.text('spheres.txt')
pile2=O.bodies.append(pile)
#sup()

################################
###  MEMBRANE GENERATION     ###
################################
mesh=2
if(mesh==1):
	#Create all nodes first :
	for i in range(0,nbL+1):
		for j in range(0,nbl):
			if (i%2==0):
				z=i*L/float(nbL)
				y=l*sin(2*pi*j/float(nbl))
				x=l*cos(2*pi*j/float(nbl))
				nodesIds.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )
			else:
				z=i*L/float(nbL)
				y=0.5*l*(sin(2*pi*j/float(nbl))+sin(2*pi*(j+1)/float(nbl)))
				x=0.5*l*(cos(2*pi*j/float(nbl))+cos(2*pi*(j+1)/float(nbl)))
				
				
				nodesIds1.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )

	##Create connection between the nodes
	for i in range(0,(nbL+1)/2+1,1):
		for j in range(0,nbl-1):
			O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[i*nbl+j+1],r,color=color) )
	for i in range(0,(nbL+1)/2,1):	
			for j in range(0,nbl):
				O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],r,color=color) )
	for i in range(-1,(nbL+1)/2):
		j=nbl
		O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color) )
		
	for i in range(0,(nbL+1)/2,1):
		for j in range(0,nbl,1):	
			cylIds.append(O.bodies.append( gridConnection(nodesIds[j+i*(nbl)],nodesIds1[j+i*(nbl)],r,color=color) ))
			cylIds.append(O.bodies.append( gridConnection(nodesIds[j+(i+1)*(nbl)],nodesIds1[j+i*(nbl)],r,color=color) ))
	for i in range(0,(nbL+1)/2,1):
		for j in range(0,nbl-1,1):
			print j,nbl
			cylIds.append(O.bodies.append( gridConnection(nodesIds[j+1+(i+1)*(nbl)],nodesIds1[j+i*(nbl)],r,color=color) ))
			cylIds.append(O.bodies.append( gridConnection(nodesIds[j+1+i*(nbl)],nodesIds1[j+i*(nbl)],r,color=color) ))
	for i in range(0,(nbL+1)/2,1):
			j=nbl-1	
			cylIds.append(O.bodies.append( gridConnection(nodesIds[i*nbl+j+1],nodesIds1[j+i*(nbl)],r,color=color) ))
			cylIds.append(O.bodies.append( gridConnection(nodesIds[(i-1)*nbl+j+1],nodesIds1[j+i*(nbl)],r,color=color) ))
	##Create PFacets
	#wire=True
	for i in range(0,(nbL+1)/2,1):
		for j in range(0,nbl):
			pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds1[i*nbl+j],color=color,mask=5)))
	for i in range(0,(nbL+1)/2,1):
		for j in range(0,nbl-1):		
			pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds1[i*nbl+j],nodesIds[i*nbl+j+1],color=color,mask=5)))
			pfIds.append(O.bodies.append(pfacet(nodesIds[(i+1)*nbl+j+1],nodesIds[i*nbl+j+1],nodesIds1[i*nbl+j],color=color,mask=5)))
			pfIds.append(O.bodies.append(pfacet(nodesIds[(i+1)*nbl+j+1],nodesIds1[i*nbl+j],nodesIds[(i+1)*nbl+j],color=color,mask=5)))
		
	for i in range(0,(nbL+1)/2,1):
		j=nbl-1	
		pfIds.append(O.bodies.append(pfacet( nodesIds[(i-1)*nbl+j+1],nodesIds1[i*nbl+j],nodesIds[i*nbl+j+1],color=color )))
		pfIds.append(O.bodies.append(pfacet(nodesIds[(i-1)*nbl+j+1],nodesIds[i*nbl+j],nodesIds1[i*nbl+j],color=color,mask=5)))
		pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j+1],nodesIds1[i*nbl+j],nodesIds[(i+1)*nbl+j],color=color,mask=5)))
		
if(mesh==2):
	#Create all nodes first :
	for i in range(0,nbL+1):
		for j in range(0,nbl):
			z=i*L/float(nbL)
			y=l*sin(2*pi*j/float(nbl))
			x=l*cos(2*pi*j/float(nbl))
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



########################
#### WALL GENERATION  ##
########################
topPlate=utils.wall(position=L+L/20.+radius,sense=0, axis=2,color=Vector3(1,0,0))
O.bodies.append(topPlate)
bottomPlate=utils.wall(position=-radius,sense=0, axis=2,color=Vector3(1,0,0))
O.bodies.append(bottomPlate)

###################
#### APPLY LOAD  ##
###################

#### APPLY CONFINING PRESSURE  

def Apply_confiningpressure():
	for i in pfIds:
		e0 =O.bodies[i].shape.node3.state.pos - O.bodies[i].shape.node1.state.pos
		e1 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node1.state.pos
		e2 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node3.state.pos
		P=(O.bodies[i].shape.node1.state.pos+O.bodies[i].shape.node2.state.pos+O.bodies[i].shape.node3.state.pos)/3
		#print e0,e1,e2
		#nodesIds.append( O.bodies.append(gridNode([P[0],P[1],P[2]],r,wire=False,fixed=True,material='NodeMat',color=color)) )
		#print 'P=',P
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
			
		O.forces.addF(O.bodies[i].shape.node1.id,p1normal)
		O.forces.addF(O.bodies[i].shape.node2.id,p2normal)
		O.forces.addF(O.bodies[i].shape.node3.id,p3normal)
		
		
		
		
#### MOVE TOP AND BOTTOM WALL 

def moveWall():
	topPlate.state.vel=(0,0,-0.05)
	#bottomPlate.state.vel=(0,0,0.05)
#g=-9.81
g=0
#limitfinder()
###########################
##### ENGINE DEFINITION  ##
###########################  	
O.dt=0.5*PWaveTimeStep()
O.engines=O.engines+[
	#PyRunner(iterPeriod=1,dead=False,command='Apply_confiningpressure()'),
	NewtonIntegrator(damping=0.7,gravity=(0,0,g),label='Newton'),
	]	