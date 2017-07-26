#!/usr/bin/python
# -*- coding: utf-8 -*-
from yade import qt

def ringPack(intR=1., r=0.1, h=0.1, isOrtho=True, **kw):
	"""Definition of the particles for a ring/tube along the z-axis.

	:param inR: internal radius of the ring/tube
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
	n = int( pi/asin(r/(intR+r)) ) + 1
	dtheta=2.*pi/n
	r = intR*sin(dtheta/2.)/(1.-sin(dtheta/2.))
	R=intR+r

	nz=int(h/(2.*r)) + 1

	for j in range(0,nz):
		z=r+2.*r*j
		for i in range(0,n):
			theta=dtheta*i
			x=R*cos(theta)
			y=R*sin(theta)
			ring.append(utils.sphere(center=(x,y,z),radius=r,**kw))

	return [ring,r,z]


#### properties of particles
kw = {'color':[1,0,0],'wire':False,'highlight':False,'fixed':False,}
#kw = {'color':[1,1,0],'wire':False,'highlight':False,'fixed':False,'material':0,'mask':0}

#### generate membrane
[ring,r,h] = ringPack( intR=1, r=0.1, h=4., isOrtho=True, **kw )
O.bodies.append(ring)

#### to see it
v=qt.Controller()
v=qt.View()
