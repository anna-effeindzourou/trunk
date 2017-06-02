#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from yade import plot,pack
import time, sys, os, copy

from yade import ymport, utils,pack,export,plot
import gts,os
#import getpass
#user=getpass.getuser()
#sys.path.append('/home/'+user+'/YADE-git/install/bin')
#from Functions import *
r_c=0.1
zmin=0
zmax=3*r_c
r_c=0.1
n=400
zmin=0
zmax=3*r_c
hexa=True
mask=2
r_scyl=(2*pi*r_c)/(2*n)
##########################
##	PARAMETERS	##
##########################

####membrane material
young_wang=0.6e6
poisson_wang=.25
cohesion_wang=20e3
frictionAngle_mem=0
Tensile_wang=20e100
density=1000

##sample parameters
frictionAngle_pp = 0.096
G=79e6
poisson=0.28
Kn=2*G*(1+poisson)

frictionAngle_pw = 0
#Kn_pw=1.7*10e9
Kn_pp=Kn
Kn_pw=Kn_pp
density=78000
g=0

####
intRadius=1.5
##########################
##	MATERIAL	##
##########################
membrane=O.materials.append(CFpmMat(young=young_wang,poisson=poisson_wang,density=density,frictionAngle=frictionAngle_mem))

wallmat=O.materials.append(FrictMat(density=density,young=Kn_pw,poisson=poisson,frictionAngle=frictionAngle_pw))
#membrane=O.materials.append(FrictMat(density=density,young=Kn_pw,poisson=poisson,frictionAngle=frictionAngle_pw))
samplemat=O.materials.append(FrictMat(density=density,young=Kn_pp,poisson=poisson,frictionAngle=frictionAngle_pp))

#membrane=O.materials.append(CFpmMat(young=young_wang,poisson=poisson_wang,density=density,frictionAngle=frictionAngle_wang))




sp=yade._packSpheres.SpherePack()
sp.load('fichier.txt')
sp.toSimulation()
#cyl=O.bodies.append(facetCylinder(center=Vector3(0,0,zmax/2.),radius=r_c+0.003, height=zmax-zmin))
zmax=hMaxSpheres(2)
r_c=r_c+4*r_scyl
m=mem(r_c, n,zmin-0*r_scyl, zmax+0*r_scyl,hexa,membrane,mask)
#m=utils.sphere(center=(0,0,0),color=(1,0,0),radius=r_scyl,material=membrane,mask=10)
m_id=O.bodies.append(m)
n_sp=len(m_id)
Min_Max=MinMax()
##walls
top=O.bodies.append(utils.wall(position=(0,0,zmax+Min_Max[1]), axis=2,material=wallmat,mask=8))
bottom=O.bodies.append(utils.wall(position=(0,0,zmin-Min_Max[1]), axis=2, material=wallmat,mask=8))


#####################
##     FUNCTIONS   ##
#####################
def shearVel(n):
  O.bodies[top].state.vel[2]=-n
  O.bodies[bottom].state.vel[2]=n
########################
## ENGINE DEFINITION  ##
########################   
O.dt=PWaveTimeStep()
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=intRadius,label='aabb'),Bo1_Wall_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intRadius,label='ss2sc'),
		Ig2_Sphere_Sphere_ScGeom6D(label='ss2d3dg'),
		Ig2_Wall_Sphere_ScGeom()],
		[
		  #Ip2_FrictMat_FrictMat_GravelPhys(alpha=1.0, krot=1.8, eta=1.8,psi_unloading=2, c=0,T=10000,b=1,power=True, label='part_part'),
		  Ip2_FrictMat_FrictMat_FrictPhys(),
		  Ip2_CFpmMat_CFpmMat_CFpmPhys(Alpha=0.098, Beta=1, cohesion=cohesion_wang, eta=1, strengthSoftening=1,tensileStrength=Tensile_wang, useAlphaBeta=False)
		],
		[
		  Law2_ScGeom_FrictPhys_CundallStrack(),
		  #Law2_ScGeom_GravelPhys_Plassiard(includeMoment=True, label='law'),
		  Law2_ScGeom_CFpmPhys_CohesiveFrictionalPM()
		  
		  ]),
	#PyRunner(initRun=True,iterPeriod=1,command='ajout_force()'),

	NewtonIntegrator(damping=.7,gravity=(0,0,g),label='Newton'),
	#ForceEngine(force=Vector3(Fx,Fy,12), ids=[0]),
	#PyRunner(initRun=True,iterPeriod=1,command='dataCollector()'),
	#VTKRecorder(iterPeriod=1000,initRun=True,fileName='paraview/'+O.tags['description']+'-',recorders=['spheres','velocity','intr'])
	]