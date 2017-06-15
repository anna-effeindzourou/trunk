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
	hMax=hMax+O.bodies[idHMax].shape.radius
	return (hMax)	
    
    
def hMin(n):
	idHMin=0
	hMin=100000.0
	for i in O.bodies:
		h=i.state.pos[n]
		if (h<hMin):
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
  
def limitfinder():
	for b in O.bodies:
		if(b.state.pos[2]>=L-2*radius):
			if isinstance(b.shape,GridNode):
				top_boundary.append(b.id)
				b.shape.color=(1,0,0)
				b.state.blockedDOFs='z'
		if(b.state.pos[2]<0.1*radius ):
			if isinstance(b.shape,GridNode):
				bottom_boundary.append(b.id)
				b.state.blockedDOFs='z'
				b.shape.color=(1,0,0)	
				
##############################
#####     SCRIPT          ####
##############################
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
radius=.0008
sigma=-3e6

#### Parameters of a rectangular grid ###
L=0.205 #length [m]
l=0.101/2.	#half width	(radius) [m]
nbL=36#number of nodes for the length	[#] doit etre paire
nbl=44	#number of nodes for the perimeter	[#] ABSOLUMENT MULTIPLE de 4 !!!

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
poisson=0.28
E=2*7.9e10*(1+poisson)
density=7.8e10
Et=0
frictionAngle=0.096
frictionAngleW=0.228
O.materials.append(CohFrictMat(young=E*0.1,poisson=poisson,density=density,frictionAngle=frictionAngle,normalCohesion=1e19,shearCohesion=1e19,momentRotationLaw=False,alphaKr=0,label='NodeMat'))
O.materials.append(FrictMat(young=E*0.1,poisson=poisson,density=density,frictionAngle=frictionAngle,label='Pmat'))

O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngle,label='Smat'))


##############################
###  SAMPLE GENERATION     ###
##############################
kw={'color':[1,1,1],'wire':False,'dynamic':True,'material':2}
pile=ymport.text('spheres.txt',**kw)
pile2=O.bodies.append(pile)
#sup()
print hMin(2), hMax(2)
zmin=hMin(2)
zmax=hMax(2)

#L=hMax(2)
#################################
####  MEMBRANE GENERATION     ###
#################################
#mesh=2
#if(mesh==1):
	##Create all nodes first :
	#for i in range(0,nbL+1):
		#for j in range(0,nbl):
			#if (i%2==0):
				#z=i*L/float(nbL)
				#y=l*sin(2*pi*j/float(nbl))
				#x=l*cos(2*pi*j/float(nbl))
				#nodesIds.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )
			#else:
				#z=i*L/float(nbL)
				#y=0.5*l*(sin(2*pi*j/float(nbl))+sin(2*pi*(j+1)/float(nbl)))
				#x=0.5*l*(cos(2*pi*j/float(nbl))+cos(2*pi*(j+1)/float(nbl)))
				
				
				#nodesIds1.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )

	###Create connection between the nodes
	#for i in range(0,(nbL+1)/2+1,1):
		#for j in range(0,nbl-1):
			#O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[i*nbl+j+1],r,color=color,material='Pmat',Et=Et) )
	#for i in range(0,(nbL+1)/2,1):	
			#for j in range(0,nbl):
				#O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],r,color=color,material='Pmat',Et=Et) )
	#for i in range(-1,(nbL+1)/2):
		#j=nbl
		#O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color,material='Pmat',Et=Et) )
		
	#for i in range(0,(nbL+1)/2,1):
		#for j in range(0,nbl,1):	
			#cylIds.append(O.bodies.append( gridConnection(nodesIds[j+i*(nbl)],nodesIds1[j+i*(nbl)],r,color=color,material='Pmat',Et=Et) ))
			#cylIds.append(O.bodies.append( gridConnection(nodesIds[j+(i+1)*(nbl)],nodesIds1[j+i*(nbl)],r,color=color,material='Pmat',Et=Et) ))
	#for i in range(0,(nbL+1)/2,1):
		#for j in range(0,nbl-1,1):
			#print j,nbl
			#cylIds.append(O.bodies.append( gridConnection(nodesIds[j+1+(i+1)*(nbl)],nodesIds1[j+i*(nbl)],r,color=color,material='Pmat',Et=Et) ))
			#cylIds.append(O.bodies.append( gridConnection(nodesIds[j+1+i*(nbl)],nodesIds1[j+i*(nbl)],r,color=color,material='Pmat',Et=Et) ))
	#for i in range(0,(nbL+1)/2,1):
			#j=nbl-1	
			#cylIds.append(O.bodies.append( gridConnection(nodesIds[i*nbl+j+1],nodesIds1[j+i*(nbl)],r,color=color,material='Pmat',Et=Et) ))
			#cylIds.append(O.bodies.append( gridConnection(nodesIds[(i-1)*nbl+j+1],nodesIds1[j+i*(nbl)],r,color=color,material='Pmat',Et=Et) ))
	###Create PFacets
	##wire=True
	#for i in range(0,(nbL+1)/2,1):
		#for j in range(0,nbl):
			#pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds1[i*nbl+j],color=color,mask=5,material='Pmat')))
	#for i in range(0,(nbL+1)/2,1):
		#for j in range(0,nbl-1):		
			#pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds1[i*nbl+j],nodesIds[i*nbl+j+1],color=color,mask=5,material='Pmat')))
			#pfIds.append(O.bodies.append(pfacet(nodesIds[(i+1)*nbl+j+1],nodesIds[i*nbl+j+1],nodesIds1[i*nbl+j],color=color,mask=5,material='Pmat')))
			#pfIds.append(O.bodies.append(pfacet(nodesIds[(i+1)*nbl+j+1],nodesIds1[i*nbl+j],nodesIds[(i+1)*nbl+j],color=color,mask=5,material='Pmat')))
		
	#for i in range(0,(nbL+1)/2,1):
		#j=nbl-1	
		#pfIds.append(O.bodies.append(pfacet( nodesIds[(i-1)*nbl+j+1],nodesIds1[i*nbl+j],nodesIds[i*nbl+j+1],color=color,mask=5,material='Pmat' )))
		#pfIds.append(O.bodies.append(pfacet(nodesIds[(i-1)*nbl+j+1],nodesIds[i*nbl+j],nodesIds1[i*nbl+j],color=color,mask=5,material='Pmat')))
		#pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j+1],nodesIds1[i*nbl+j],nodesIds[(i+1)*nbl+j],color=color,mask=5,material='Pmat')))
		
#if(mesh==2):
	##Create all nodes first :
	#for i in range(0,nbL+1):
		#for j in range(0,nbl):
			#z=i*L/float(nbL)
			#y=l*sin(2*pi*j/float(nbl))
			#x=l*cos(2*pi*j/float(nbl))
			#nodesIds.append( O.bodies.append(gridNode([x,y,z],r,wire=False,fixed=False,material='NodeMat',color=color)) )

	###Create connection between the nodes
	#for i in range(0,nbL+1):
		#for j in range(0,nbl-1):
			#O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[i*nbl+j+1],r,color=color,mask=5,material='Pmat',Et=Et) )
	#for i in range(0,nbL,1):	
			#for j in range(0,nbl):
				#O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],r,color=color,mask=5,material='Pmat',Et=Et) )
	#for i in range(-1,nbL):
		#j=nbl
		#O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color,mask=5,material='Pmat',Et=Et) )
		
	#for i in range(0,nbL):
		#for j in range(0,nbl-1):
				#if (j%2==0):
					#O.bodies.append( gridConnection(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j+1],r,color=color,mask=5,material='Pmat',Et=Et) )
				#else:
					#O.bodies.append( gridConnection(nodesIds[(i+1)*nbl+j],nodesIds[i*nbl+j+1],r,color=color,mask=5,material='Pmat',Et=Et) )
	#for i in range(0,nbL):
		#j=nbl
		##O.bodies[nodesIds[(i-1)*nbl+j]].shape.color=Vector3(155./255.,155./255.,1.)
		##O.bodies[nodesIds[(i)*nbl+j-1]].shape.color=Vector3(1,0,0)
		#O.bodies.append( gridConnection(nodesIds[(i-1)*nbl+j],nodesIds[(i+1)*nbl+j-1],r,color=color,mask=5,material='Pmat',Et=Et) )		


	####Create PFacets
	###wire=True
	#for i in range(0,nbL):
		#for j in range(0,nbl-1):
				#if (j%2==0):
					#pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],color=color,mask=5,material='Pmat')))
					#pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j+1],nodesIds[(i)*nbl+j+1],color=color,mask=5,material='Pmat')))
				#else:
					#pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j],nodesIds[(i+1)*nbl+j],nodesIds[(i)*nbl+j+1],color=color,mask=5,material='Pmat')))
					#pfIds.append(O.bodies.append(pfacet(nodesIds[i*nbl+j+1],nodesIds[(i+1)*nbl+j],nodesIds[(i+1)*nbl+j+1],color=color,mask=5,material='Pmat')))
		
	#for i in range(0,nbL,1):
		#j=nbl
		#pfIds.append(O.bodies.append(pfacet( nodesIds[i*nbl+j],nodesIds[(i-1)*nbl+j],nodesIds[(i+1)*nbl+j-1],color=color,material='Pmat' )))
		#pfIds.append(O.bodies.append(pfacet( nodesIds[(i)*nbl+j-1],nodesIds[(i+1)*nbl+j-1],nodesIds[(i-1)*nbl+j],color=color,material='Pmat' )))


#limitfinder()

#########################
##### WALL GENERATION  ##
#########################
#O.materials.append(FrictMat(young=E,poisson=poisson,density=density,frictionAngle=frictionAngleW,label='Wmat'))
#topPlate=utils.wall(position=hMax(2)+radius,sense=0, axis=2,color=Vector3(1,0,0),material='Wmat')
#O.bodies.append(topPlate)
#bottomPlate=utils.wall(position=-hMin(2)-radius,sense=0, axis=2,color=Vector3(1,0,0),material='Wmat')
#O.bodies.append(bottomPlate)

####################
##### APPLY LOAD  ##
####################

##### APPLY CONFINING PRESSURE  

#def Apply_confiningpressure():
	##print 'Apply_confiningpressure'
	#for i in pfIds:
		#e0 =O.bodies[i].shape.node3.state.pos - O.bodies[i].shape.node1.state.pos
		#e1 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node1.state.pos
		#e2 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node3.state.pos
		#P=(O.bodies[i].shape.node1.state.pos+O.bodies[i].shape.node2.state.pos+O.bodies[i].shape.node3.state.pos)/3
		##print e0,e1,e2
		##nodesIds.append( O.bodies.append(gridNode([P[0],P[1],P[2]],r,wire=False,fixed=True,material='NodeMat',color=color)) )
		##print 'P=',P
		#v0 = e0
		#v1 = e1
		#v2 = P - O.bodies[i].shape.node1.state.pos
		
		###// Compute dot products
		#dot00 = scalar(v0,v0)
		#dot01 = scalar(v0,v1)
		#dot02 = scalar(v0,v2)
		#dot11 = scalar(v1,v1)
		#dot12 = scalar(v1,v2)
		
		###// Compute the barycentric coordinates of the projection P
		#invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
		#p1 = (dot11 * dot02 - dot01 * dot12) * invDenom
		#p2 = (dot00 * dot12 - dot01 * dot02) * invDenom
		#p3 = 1-p1-p2

		#a = sqrt(scalar(e0,e0))
		#b = sqrt(scalar(e1,e1))
		#c = sqrt(scalar(e2,e2))

		#s=0.5*(a+b+c)
		#area= sqrt(s*(s-a)*(s-b)*(s-c))
		
		#Fapplied=area*sigma
		
		#normal = cross(e0,e1)	
		
		#normal=normal/normal.norm()
		
		#F=Fapplied
		#p1normal=F*p1*normal
		#p2normal=F*p2*normal
		#p3normal=F*p3*normal
			
		#O.forces.addF(O.bodies[i].shape.node1.id,p1normal,permanent=False)
		#O.forces.addF(O.bodies[i].shape.node2.id,p2normal,permanent=False)
		#O.forces.addF(O.bodies[i].shape.node3.id,p3normal,permanent=False)
##Apply_confiningpressure()

#sigma3=0
#def check_confiningpressure():
	#global sigma3
	#sigma3=0
	#for i in pfIds:
		#e0 =O.bodies[i].shape.node3.state.pos - O.bodies[i].shape.node1.state.pos
		#e1 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node1.state.pos
		#e2 =O.bodies[i].shape.node2.state.pos - O.bodies[i].shape.node3.state.pos

		#a = sqrt(scalar(e0,e0))
		#b = sqrt(scalar(e1,e1))
		#c = sqrt(scalar(e2,e2))

		#s=0.5*(a+b+c)
		#area= sqrt(s*(s-a)*(s-b)*(s-c))
		#F=(O.forces.f(O.bodies[i].shape.node1.id) + O.forces.f(O.bodies[i].shape.node2.id)+O.forces.f(O.bodies[i].shape.node3.id)).norm()
		#sigma3=sigma3+F/area
	##print sigma3
	#return sigma3
#pos=topPlate.state.pos[2]
#def dataCollector():
	#global pos
	#if(pos<0.1845):
			#O.wait()
			#saveData()
			
	#S=pi*l**2
	#Fnt=O.forces.f(topPlate.id)[2]
	#Fnb=O.forces.f(bottomPlate.id)[2]
	#sigma1=Fnt/S
	#sigma3=check_confiningpressure()
	#pos=topPlate.state.pos[2]
	#plot.addData(t=O.time,pos=pos,Fnt=Fnt,Fnb=Fnb,sigma1=sigma1,sigma3=sigma3,unbF=unbalancedForce())
    
    
#def saveData():
	#plot.saveDataTxt('triaxial_res.dat',vars=('t','pos','Fnt','Fnb','sigma1','sigma3','unbF'))

#plot.plots={'t':('sigma1',Et,'sigma3')}




		
##### MOVE TOP AND BOTTOM WALL 
#v=1.7e-03
##v=1
#def moveWall(v):
	#topPlate.state.vel=(0,0,-v)
	##bottomPlate.state.vel=(0,0,v)
##g=-9.81
#g=0
##moveWall(v)
##limitfinder()
############################
###### ENGINE DEFINITION  ##
############################  	
#O.dt=0.5*PWaveTimeStep()
#O.engines=O.engines+[
	#PyRunner(iterPeriod=1,dead=False,command='Apply_confiningpressure()'),
	#NewtonIntegrator(damping=0.7,gravity=(0,0,g),label='Newton'),
	#PyRunner(initRun=True,iterPeriod=1,command='dataCollector()'),
	#VTKRecorder(iterPeriod=500,initRun=True,fileName='paraview/'+'triaxial_res-',recorders=['spheres','velocity']),
	#]	

#if not isBatch:
	##  VISUALIZATION  
	#from yade import qt
	#qt.Controller()
	##qtv = qt.View()
	##qtr = qt.Renderer()
	#plot.plot(noShow=False, subPlots=True)
	#O.run(5000)
	##moveWall(v)
#else:
  #O.run(5000,True)
  #moveWall(v)
  #O.wait()
  #saveData()