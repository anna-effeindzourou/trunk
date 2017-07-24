#!/usr/bin/python
# -*- coding: utf-8 -*-
from yade import qt
from yade.gridpfacet import *

def ringPack(intR=1., r=0.1,mosx=0,mosz=0, h=0.1, **kw):
	"""Definition of the particles for a ring/tube along the z-axis.

	:param intR: internal radius of the ring/tube
	:param r: approximate radius of the particles which will be used to generate the ring/tube
	:param h: height in +z-direction (min z will be 0, i.e. first ring will be generated a +r)
	:param isOrtho: if true the spheres are generated in a regular orthogonal grid, if false in a regular hexagonal grid
	:param \*\*kw: passed to :yref:`yade.utils.sphere`

	:return: set of spheres which defines the ring/tube and exact dimensions of the ring tube (r,h).
	
	note::
	This function will use the input value for r as first approximation to calculate the number of particles in one ring. The calculated number is rounded to the next integer and a new effective r is calculated. The calculated value for r is always smaller then the input value.

	"""

	# check input dimension
	if(intR<=r): raise ValueError("intR must be greater than r!");

	ring=[]
	nx=int(2*pi*intR/mosx)
	print nx
	dtheta=2.*pi/nx
	R=intR+r

	nz=int(h/mosz) + 2
	for j in range(0,nz):
		#z=r+2.*r*j*sin(pi/3)
		z=mosz*j
		print z
		if(j%2==0):
			for i in range(0,nx):
				theta=dtheta*i
				x=R*cos(theta)
				y=R*sin(theta)
				ring.append(gridNode(center=(x,y,z),radius=r,**kw))
			
		else:
			for i in range(0,nx):
				theta=dtheta*(i+.5)
				x=R*cos(theta)
				y=R*sin(theta)
				ring.append(gridNode(center=(x,y,z),radius=r,**kw))
	return [ring,r,z,nx,nz]    


###################
#####   Box   #####
###################
## material properties
poisson = 1
fricangle = 30

#Lx = 1.079

strainStressValues=[(10e-03,20e03),(20e-03,60e03),(30e-03,100e03),(40e-03,130e03),(58.93e-03,165.7e03)]
d=3./1000
radius = d*4.

particleVolume = 4./3.*pow(radius,3)*pi
particleMass = 1.65/1740.
density = particleMass/particleVolume
young = strainStressValues[0][1] / strainStressValues[0][0]

meshMatInternal = O.materials.append( WireMat( young=young,poisson=poisson,frictionAngle=radians(fricangle),density=density,isDoubleTwist=False,diameter=d,strainStressValues=strainStressValues ) )

meshMat = O.materials.append( FrictMat( young=young,poisson=poisson,frictionAngle=radians(fricangle),density=density, label='boxMat' ) )



#### properties of particles
kw = {'color':[1,0,0],'wire':False,'highlight':False,'fixed':False,}
#kw = {'color':[1,1,0],'wire':False,'highlight':False,'fixed':False,'material':0,'mask':0}
d=3./1000
radius = d*4.
#### generate membrane
[ring,r,h,nx,nz] = ringPack( intR=0.5, r=radius, h=0.6,mosx=0.083,mosz=0.0715, **kw )
O.bodies.append(ring)




netcyl=[]
for j in range(0,nz-1):#ny-1
	cx=nx
	if(j%2==0):
			for i in range(1,cx):	
				netcyl+=[gridConnection(i+j*(nx),i+nx+j*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]
				netcyl+=[gridConnection(i+j*(nx),i+nx-1+j*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]
	else:
		for i in range(-1,cx):
			netcyl+=[gridConnection(i+j*(nx),i+nx+1+j*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]
		for i in range(0,cx):	
			netcyl+=[gridConnection(i+j*(nx),i+nx+j*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]	

i=0
for j in range(1,nz-1):
	if(j%2==1):
		netcyl+=[gridConnection(i+j*(nx),i-1+nx+(j+2)*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]

i=0
j=0	
netcyl+=[gridConnection(i+j*(nx),i-1+nx+(j+1)*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]
netcyl+=[gridConnection(i+j*(nx),i+nx+(j)*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]	
netcyl+=[gridConnection(i+nx+(j+1)*(nx),i+nx+(j+2)*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]

netcylghost=[]
for j in range(1,nz-1):#ny-1
	for i in range(0,cx-1):		netcylghost+=[gridConnection(i+j*(nx),i+1+j*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]
			
for j in range(1,nz-1):
	i=cx-1
	netcylghost+=[gridConnection(i+j*(nx),i+1+(j-1)*(nx),radius=radius,color=Vector3(1,1,1),material=meshMat)]


limit=[]
j=0
for i in range(0,cx-1):	limit+=[gridConnection(i+j*(nx),i+1+j*(nx),radius=radius,color=Vector3(1,0,1),material=meshMat)]
i=cx-1
limit+=[gridConnection(i+j*(nx),i+1+(j-1)*(nx),radius=radius,color=Vector3(1,0,1),material=meshMat)]



j=nz-1
for i in range(0,cx-1):		limit+=[gridConnection(i+j*(nx),i+1+j*(nx),radius=radius,color=Vector3(1,0,1),material=meshMat)]
i=cx-1
limit+=[gridConnection(i+j*(nx),i+1+(j-1)*(nx),radius=radius,color=Vector3(1,0,1),material=meshMat)]








netcylId=O.bodies.append(netcyl)	
netcylIdghost=O.bodies.append(netcylghost)	
limitId=O.bodies.append(limit)					
				

#### to see it
v=qt.Controller()
v=qt.View()
