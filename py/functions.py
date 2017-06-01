# encoding: utf-8
# 2015 © Anna Effeindzourou <anna.effeindzourou@gmail.com>


# -*- coding: utf-8
from yade import ymport, utils,pack,export,qt



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

def hMax(n):
	idHMax=0
	hMax=-1000000.0
	for i in O.bodies:
		h = i.state.pos[n]
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



def hMinSpheres(n):
	idHMin=0
	hMin=100000.0
	for i in O.bodies:
		h=i.state.pos[n]-i.shape.radius
		if (type(i.shape)==Sphere):
			if (h<hMin):
				hMin=h
				idHMin=i.id
	res=[idHMin,hMin]
	return (res) 
	
def hMaxSpheres(n):
	idHMax=0
	hMax=-1000000.0
	for i in O.bodies:
		h=i.state.pos[n]+i.shape.radius
		if (type(i.shape)==Sphere):
			if (h>hMax):
				hMax=h
				idHMax=i.id
	res=[idHMax,hMax]
	return (res)


#Function in order to calculate rmin (minimum radius) and rmax (maximum radius)
def RMinMax():
	rmax=0
	rmin=10
	r=0
	for i in O.bodies:
		if(type(i.shape)==Sphere):
			r=i.shape.radius
			if(r>rmax): rmax=r
			if(r<rmin): rmin=r
	l=[rmin,rmax]
	return (l)
 
def sup():
	for i in O.bodies:
		if (type(i.shape)==Sphere) and (i.state.pos[2]>0.098): O.bodies.erase(i.id)    
  
def scalar(u,v):
	ps=u[0]*v[0]+u[1]*v[1]+u[2]*v[2]
	return ps

def cross(u,v):
	ps=Vector3(u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])
	return ps  
  
def limitfinder(DOFs):
	for b in O.bodies:
		if(b.state.pos[2]>=L-2*radius):
			if isinstance(b.shape,GridNode):
				top_boundary.append(b.id)
				b.shape.color=(1,0,0)
				b.state.blockedDOFs=DOFs
		if(b.state.pos[2]<0.1*radius):
			if isinstance(b.shape,GridNode):
				bottom_boundary.append(b.id)
				b.state.blockedDOFs = DOFs
				b.shape.color=(1,0,0)	
				
				
				
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
  


			