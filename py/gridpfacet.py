# encoding: utf-8
#
# 2015 © Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>
# 2015 © Anna Effeindzourou <anna.effeindzourou@gmail.com>
# 2015 © François Kneib <francois.kneib@gmail.com>
# 2015 © Klaus Thoeni <klaus.thoeni@gmail.com>

"""
Helper functions for creating cylinders, grids and membranes. For more details on this type of elements see [Effeindzourou2016]_, [Effeindzourou2015a]_, [Bourrier2013]_,.

For examples using :yref:`GridConnections<GridConnection>`, see

* :ysrc:`examples/grids/CohesiveGridConnectionSphere.py`
* :ysrc:`examples/grids/GridConnection_Spring.py`
* :ysrc:`examples/grids/Simple_Grid_Falling.py`
* :ysrc:`examples/grids/Simple_GridConnection_Falling.py`

For examples using :yref:`PFacets<PFacet>`, see

* :ysrc:`examples/pfacet/gts-pfacet.py`
* :ysrc:`examples/pfacet/mesh-pfacet.py`
* :ysrc:`examples/pfacet/pfacetcreators.py`

"""

import math,random,doctest,geom,numpy
from yade.wrapper import *
try: # use psyco if available
	import psyco
	psyco.full()
except ImportError: pass

from yade import utils
from yade._utils import createInteraction

from minieigen import *

def gridNode(center,radius,dynamic=None,fixed=False,wire=False,color=None,highlight=False,material=-1):
	"""
	Create a :yref:`GridNode` which is needed to set up :yref:`GridConnections<GridConnection>`.

	See documentation of :yref:`yade.utils.sphere` for meaning of parameters.

	:return: Body object with the :yref:`gridNode` :yref:`shape<Body.shape>`.
	"""
	b=Body()
	b.shape=GridNode(radius=radius,color=color if color else utils.randomColor(),wire=wire,highlight=highlight)
	#V=(4./3)*math.pi*radius**3	# will be overwritten by the connection
	V=0.
	geomInert=(2./5.)*V*radius**2	# will be overwritten by the connection
	utils._commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=center,dynamic=dynamic,fixed=fixed)
	b.aspherical=False
	b.bounded=False
	b.mask=0	# avoid contact detection with the nodes. Manual interaction will be set for them in "gridConnection" below.
	return b


def gridConnection(id1,id2,radius,wire=False,color=None,highlight=False,material=-1,mask=1,cellDist=None):
	"""
	Create a :yref:`GridConnection` by connecting two :yref:`GridNodes<GridNode>`.

	:param id1,id2: the two :yref:`GridNodes<GridNode>` forming the cylinder. 
	:param float radius: radius of the cylinder. Note that the radius needs to be the same as the one for the :yref:`GridNodes<GridNode>`.
	:param Vector3 cellDist: for periodic boundary conditions, see :yref:`Interaction.cellDist`. Note: periodic boundary conditions for gridConnections are not yet implemented! 

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:return: Body object with the :yref:`GridConnection` :yref:`shape<Body.shape>`.

	.. note:: The material of the :yref:`GridNodes<GridNode>` will be used to set the constitutive behaviour of the internal connection, i.e., the constitutive behaviour of the cylinder. The material of the :yref:`GridConnection` is used for interactions with other (external) bodies.
	"""
	b=Body()
	b.shape=GridConnection(radius=radius,color=color if color else utils.randomColor(),wire=wire,highlight=highlight)
	sph1=O.bodies[id1] ; sph2=O.bodies[id2]
	i=createInteraction(id1,id2)
	nodeMat=sph1.material
	b.shape.node1=sph1 ; b.shape.node2=sph2
	sph1.shape.addConnection(b) ; sph2.shape.addConnection(b)
	if(O.periodic):
		if(cellDist!=None):
			i.cellDist=cellDist
		segt=sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos
	else: segt=sph2.state.pos - sph1.state.pos
	L=segt.norm()
	V=0.5*L*math.pi*radius**2
	geomInert=(2./5.)*V*radius**2
	utils._commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=sph1.state.pos,dynamic=False,fixed=True)
	sph1.state.mass = sph1.state.mass + V*nodeMat.density
	sph2.state.mass = sph2.state.mass + V*nodeMat.density
	for k in [0,1,2]:
		sph1.state.inertia[k] = sph1.state.inertia[k] + geomInert*nodeMat.density
		sph2.state.inertia[k] = sph2.state.inertia[k] + geomInert*nodeMat.density
	b.aspherical=False
	if O.periodic:
		i.phys.unp= -(sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius
		b.shape.periodic=True
		b.shape.cellDist=i.cellDist
	else:
		i.phys.unp= -(sph2.state.pos - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius	
	i.geom.connectionBody=b
	I=math.pi*(2.*radius)**4/64.
	E=nodeMat.young
	i.phys.kn=E*math.pi*(radius**2)/L
	i.phys.kr=E*I/L
	i.phys.ks=12.*E*I/(L**3)
	G=E/(2.*(1+nodeMat.poisson))
	i.phys.ktw=2.*I*G/L
	b.mask=mask
	return b


#TODO: find a better way of handling the Id lists for checking duplicated gridNodes or gridConnections with the same coordinates etc. It would be better to handle this globally, maybe implement something like O.bodies.getGridNodes
def cylinder(begin=Vector3(0,0,0),end=Vector3(1.,0.,0.),radius=0.2,nodesIds=[],cylIds=[],dynamic=None,fixed=False,wire=False,color=None,highlight=False,intMaterial=-1,extMaterial=-1,mask=1):
	"""
	Create a cylinder with given parameters. The shape corresponds to the Minkowski sum of line-segment and sphere, hence, the cylinder has rounded vertices. The cylinder (:yref:`GridConnection<GridConnection>`) and its corresponding nodes (yref:`GridNodes<GridNode>`) are automatically added to the simulation. The lists with nodes and cylinder ids will be updated automatically.

	:param Vector3 begin: first point of the Minkowski sum in the global coordinate system.
	:param Vector3 end: last point of the Minkowski sum in the global coordinate system.
	:param Real radius: radius of sphere in the Minkowski sum.
	:param list nodesIds: list with ids of already existing :yref:`GridNodes<GridNode>`. New ids will be added.
	:param list cylIds: list with ids of already existing :yref:`GridConnections<GridConnection>`. New id will be added.
	:param intMaterial: :yref:`Body.material` used to create the interaction physics between the two GridNodes
	:param extMaterial: :yref:`Body.material` used to create the interaction physics between the Cylinder (GridConnection) and other bodies (e.g., spheres interaction with the cylinder)
	
	See :yref:`yade.utils.sphere`'s documentation for meaning of other parameters.

	"""
	id1 = O.bodies.append( gridNode(begin,radius,dynamic=dynamic,fixed=fixed,wire=wire,color=color,highlight=highlight,material=intMaterial) )
	nodesIds.append(id1)
	id2 = O.bodies.append( gridNode(end,radius,dynamic=dynamic,fixed=fixed,wire=wire,color=color,highlight=highlight,material=intMaterial) )
	nodesIds.append(id2)
	cylIds.append(O.bodies.append( gridConnection(id1,id2,radius=radius,wire=wire,color=color,highlight=highlight,material=extMaterial,mask=mask,cellDist=None) ))


def cylinderConnection(vertices,radius=0.2,nodesIds=[],cylIds=[],dynamic=None,fixed=False,wire=False,color=None,highlight=False,intMaterial=-1,extMaterial=-1,mask=1):
	"""
	Create a chain of cylinders with given parameters. The cylinders (:yref:`GridConnection<GridConnection>`) and its corresponding nodes (yref:`GridNodes<GridNode>`) are automatically added to the simulation. The lists with nodes and cylinder ids will be updated automatically.

	:param [[Vector3]] vertices: coordinates of vertices to connect in the global coordinate system.
	
	See :yref:`yade.gridpfacet.cylinder` documentation for meaning of other parameters.

	"""
	# create all gridNodes first
	nodesIdsCC=[]
	for i in vertices:
		nodesIdsCC.append( O.bodies.append(gridNode(i,radius=radius,
						dynamic=dynamic,fixed=fixed,wire=wire,color=color,highlight=highlight,material=intMaterial)) )
	nodesIds.extend(nodesIdsCC)
	# now create connection between the gridNodes
	for i,j in zip( nodesIdsCC[:-1], nodesIdsCC[1:]):
		cylIds.append( O.bodies.append( gridConnection(i,j,radius=radius,
						wire=wire,color=color,highlight=highlight,material=intMaterial,mask=mask,cellDist=None)) )


def pfacet(id1,id2,id3,wire=True,color=None,highlight=False,material=-1,mask=1,cellDist=None):
	"""
	Create a :yref:`PFacet<PFacet>` element from 3 :yref:`GridNodes<GridNode>` which are already connected via 3 :yref:`GridConnections<GridConnection>`:
	
	:param id1,id2,id3: already with :yref:`GridConnections<GridConnection>` connected :yref:`GridNodes<GridNode>`
	:param bool wire: if ``True``, top and bottom facet are shown as skeleton; otherwise facets are filled.
	:param Vector3-or-None color: color of the PFacet; random color will be assigned if ``None``.
	:param Vector3 cellDist: for periodic boundary conditions, see :yref:`Interaction.cellDist`. Note: periodic boundary conditions are not yet implemented for PFacets! 

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:return: Body object with the :yref:`PFacet<PFacet>` :yref:`shape<Body.shape>`.

	.. note:: :yref:`GridNodes<GridNode>` and :yref:`GridConnections<GridConnection>` need to have the same radius. This is also the radius used to create the :yref:`PFacet<PFacet>`

	"""
	b=Body()
	GridN1=O.bodies[id1]; GridN2=O.bodies[id2]; GridN3=O.bodies[id3] 
	b.shape=PFacet(color=color if color else randomColor(),wire=wire,highlight=highlight,node1=GridN1,node2=GridN2,node3=GridN3)
	GridN1.bounded=False; GridN2.bounded=False; GridN3.bounded=False
	GridC1=O.bodies[O.interactions[id1,id2].geom.connectionBody.id]
	GridC2=O.bodies[O.interactions[id2,id3].geom.connectionBody.id]
	GridC3=O.bodies[O.interactions[id1,id3].geom.connectionBody.id]
	GridC1.bounded=False
	GridC2.bounded=False
	GridC3.bounded=False
	
	b.shape.conn1=GridC1
	b.shape.conn2=GridC2 
	b.shape.conn3=GridC3
	
	b.shape.radius=GridN1.shape.radius
	GridC1.shape.addPFacet(b) 
	GridC2.shape.addPFacet(b)
	GridC3.shape.addPFacet(b)
	GridN1.shape.addPFacet(b); GridN2.shape.addPFacet(b); GridN3.shape.addPFacet(b)
	
	V=0
	
	utils._commonBodySetup(b,V,Vector3(0,0,0),material,pos=GridN1.state.pos,dynamic=False,fixed=True)
	b.aspherical=False # mass and inertia are lumped into the GridNodes
	b.mask=mask
	return b


#TODO: find a better way of handling the Id lists for checking duplicated gridNodes or gridConnections with the same coordinates etc. It would be better to handle this globally, maybe implement something like O.bodies.getGridNodes
def pfacetCreator1(vertices,radius,nodesIds=[],cylIds=[],pfIds=[],wire=False,fixed=True,materialNodes=-1,material=-1,color=None):
	"""
	Create a :yref:`PFacet<PFacet>` element from 3 vertices and automatically append to simulation. The function uses the vertices to create :yref:`GridNodes<GridNode>` and automatically checks for existing nodes.
	
	:param [Vector3,Vector3,Vector3] vertices: coordinates of vertices in the global coordinate system.
	:param float radius: radius used to create the :yref:`PFacets<PFacet>`.
	:param list nodesIds: list with ids of already existing :yref:`GridNodes<GridNode>`. New ids will be added.
	:param list cylIds: list with ids of already existing :yref:`GridConnections<GridConnection>`. New ids will be added.
	:param list pfIds: list with ids of already existing :yref:`PFacets<PFacet>`. New ids will be added. 
	:param materialNodes: specify :yref:`Body.material` of :yref:`GridNodes<GridNode>`. This material is used to make the internal connections.
	:param material: specify :yref:`Body.material` of :yref:`PFacets<PFacet>`. This material is used for interactions with external bodies.
	
	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.
	"""
	n=len(nodesIds)
	k=[0,0,0]
	f=[0,0,0]
	u=0
	nod=0
	for i in vertices:
		u=0
		for j in nodesIds:
			if(i==O.bodies[j].state.pos):
				f[nod]=j
				k[nod]=1
				u+=1
		nod+=1
		test=True
		#if(u==0):
		for GN in nodesIds:
			if(i==O.bodies[GN].state.pos):
				  u=1
		if(u==0):
			nodesIds.append( O.bodies.append(gridNode(i,radius,wire=wire,fixed=fixed,material=materialNodes,color=color)) )

	if(k==[0,0,0]):
		pfacetCreator3(nodesIds[n],nodesIds[n+1],nodesIds[n+2],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )
	if(k==[1,0,0]):
		pfacetCreator3(f[0],nodesIds[n],nodesIds[n+1],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )
	if(k==[0,1,0]):
		pfacetCreator3(nodesIds[n],f[1],nodesIds[n+1],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )
	if(k==[0,0,1]):
		pfacetCreator3(nodesIds[n],nodesIds[n+1],f[2],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )
	if(k==[1,1,0]):
		pfacetCreator3(f[0],f[1],nodesIds[n],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )
	if(k==[0,1,1]):
		pfacetCreator3(nodesIds[n],f[1],f[2],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )
	if(k==[1,0,1]):
		pfacetCreator3(f[0],nodesIds[n],f[2],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )
	if(k==[1,1,1]):
		pfacetCreator3(f[0],f[1],f[2],cylIds=cylIds,pfIds=pfIds,wire=wire,material=material,color=color,fixed=fixed )


def pfacetCreator2(id1,id2,vertex,radius,nodesIds=[],wire=True,materialNodes=-1,material=-1,color=None,fixed=True):
	"""
	Create a :yref:`PFacet<PFacet>` element from 2 already existing and connected :yref:`GridNodes<GridNode>` and one vertex. The element is automatically appended to the simulation.
	
	:param int id1,id2: ids of already with :yref:`GridConnection` connected :yref:`GridNodes<GridNode>`.
	:param Vector3 vertex: coordinates of the vertex in the global coordinate system.
	
	See documentation of :yref:`yade.gridpfacet.pfacetCreator1` for meaning of other parameters.
	"""  
	n=len(nodesIds)
	nodesIds.append( O.bodies.append(gridNode(vertex,radius,wire=wire,fixed=fixed,material=materialNodes,color=color)) )
	O.bodies.append(gridConnection(id1,nodesIds[n],radius=radius,material=materialNodes,color=color,wire=wire))
	O.bodies.append(gridConnection(id2,nodesIds[n],radius=radius,material=materialNodes,color=color,wire=wire))
	O.bodies.append(pfacet(id1,id2,nodesIds[n],wire=wire,material=material,color=color))


def pfacetCreator3(id1,id2,id3,cylIds=[],pfIds=[],wire=True,material=-1,color=None,fixed=True,mask=-1):
	"""
	Create a :yref:`PFacet` element from 3 already existing :yref:`GridNodes<GridNode>` which are not yet connected. The element is automatically appended to the simulation.
	
	:param int id1,id2,id3: id of the 3 :yref:`GridNodes<GridNode>` forming the :yref:`PFacet`.
	
	See documentation of :yref:`yade.gridpfacet.pfacetCreator1` for meaning of other parameters.
	"""  
	radius=O.bodies[id1].shape.radius
	try:
		cylIds.append(O.bodies.append(gridConnection(id1,id2,radius=radius,material=material,color=color,wire=wire,mask=mask)))
	except:
		pass
	try:
		cylIds.append(O.bodies.append(gridConnection(id2,id3,radius=radius,material=material,color=color,wire=wire,mask=mask)))
	except:
		pass	      
	try:
		cylIds.append(O.bodies.append(gridConnection(id3,id1,radius=radius,material=material,color=color,wire=wire,mask=mask)))
	except:
		pass	      
	pfIds.append(O.bodies.append(pfacet(id1,id2,id3,wire=wire,material=material,color=color,mask=mask)))


def pfacetCreator4(id1,id2,id3,pfIds=[],wire=True,material=-1,color=None,fixed=True,mask=-1):
	"""
	Create a :yref:`PFacet<PFacet>` element from 3 already existing :yref:`GridConnections<GridConnection>`. The element is automatically appended to the simulation.
	
	:param int id1,id2,id3: id of the 3 :yref:`GridConnections<GridConnection>` forming the :yref:`PFacet`.
	
	See documentation of :yref:`yade.gridpfacet.pfacetCreator1` for meaning of other parameters.
	"""  
	radius=O.bodies[id1].shape.radius
	GridN=[]
	GridN1=O.bodies[id1].shape.node1.id
	GridN2=O.bodies[id1].shape.node2.id
	GridN.append(GridN1)
	GridN.append(GridN2)
	
	GridN1=O.bodies[id2].shape.node1.id
	if(GridN1 not in GridN):
		GridN.append(GridN1)
		
	GridN2=O.bodies[id2].shape.node2.id
	if(GridN2 not in GridN):
		GridN.append(GridN2)
		
	GridN1=O.bodies[id3].shape.node1.id
	if(GridN1 not in GridN):
		GridN.append(GridN1)
		
	GridN2=O.bodies[id3].shape.node2.id
	if(GridN2 not in GridN):
		GridN.append(GridN2)
	pfIds.append(O.bodies.append(pfacet(GridN[0],GridN[1],GridN[2],wire=wire,material=material,color=color,mask=mask)))


def gtsPFacet(meshfile,shift=Vector3.Zero,scale=1.0,
							radius=1,wire=True,fixed=True,materialNodes=-1,material=-1,color=None):
	"""
	Imports mesh geometry from .gts file and automatically creates connected :yref:`PFacet3<PFacet>` elements. For an example see :ysrc:`examples/pfacet/gts-pfacet.py`.

	:param string filename: .gts file to read.
	:param [float,float,float] shift: [X,Y,Z] parameter shifts the mesh.
	:param float scale: factor scales the mesh.
	:param float radius: radius used to create the :yref:`PFacets<PFacet>`.
	:param materialNodes: specify :yref:`Body.material` of :yref:`GridNodes<GridNode>`. This material is used to make the internal connections.
	:param material: specify :yref:`Body.material` of :yref:`PFacets<PFacet>`. This material is used for interactions with external bodies.

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:returns: lists of :yref:`GridNode<GridNode>` ids `nodesIds`, :yref:`GridConnection<GridConnection>` ids `cylIds`, and :yref:`PFacet<PFacet>` ids `pfIds`
	"""
	import gts,yade.pack
	surf=gts.read(open(meshfile))
	surf.scale(scale,scale,scale)
	surf.translate(shift[0],shift[1],shift[2]) 
	nodesIds=[]; cylIds=[]; pfIds=[]
	for face in surf.faces():
		a=face.vertices()[0].coords()
		b=face.vertices()[1].coords()
		c=face.vertices()[2].coords()
		pfacetCreator1([a,b,c],radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,wire=wire,fixed=fixed,materialNodes=materialNodes,material=material,color=color)
		#print a,b,c
	return nodesIds,cylIds,pfIds


def gmshPFacet(meshfile="file.mesh",shift=Vector3.Zero,scale=1.0,orientation=Quaternion.Identity,
							 radius=1.0,wire=True,fixed=True,materialNodes=-1,material=-1,color=None):
	"""
	Imports mesh geometry from .mesh file and automatically creates connected PFacet elements. For an example see :ysrc:`examples/pfacet/mesh-pfacet.py`.

	:param string filename: .gts file to read.
	:param [float,float,float] shift: [X,Y,Z] parameter shifts the mesh.
	:param float scale: factor scales the mesh.
	:param quaternion orientation: orientation of the imported geometry.
	:param float radius: radius used to create the :yref:`PFacets<PFacet>`.
	:param materialNodes: specify :yref:`Body.material` of :yref:`GridNodes<GridNode>`. This material is used to make the internal connections.
	:param material: specify :yref:`Body.material` of :yref:`PFacets<PFacet>`. This material is used for interactions with external bodies.

	See documentation of :yref:`yade.utils.sphere` for meaning of other parameters.

	:returns: lists of :yref:`GridNode<GridNode>` ids `nodesIds`, :yref:`GridConnection<GridConnection>` ids `cylIds`, and :yref:`PFacet<PFacet>` ids `pfIds`
	
	mesh files can easily be created with `GMSH <http://www.geuz.org/gmsh/>`_.
	
	Additional examples of mesh-files can be downloaded from 
	http://www-roc.inria.fr/gamma/download/download.php
	"""	
	infile = open(meshfile,"r")
	lines = infile.readlines()
	infile.close()

	nodelistVector3=[]
	findVerticesString=0
	
	while (lines[findVerticesString].split()[0]<>'Vertices'): # find the string with the number of Vertices
		findVerticesString+=1
	findVerticesString+=1
	numNodes = int(lines[findVerticesString].split()[0])
	
	for i in range(numNodes):
		nodelistVector3.append(Vector3(0.0,0.0,0.0))
	id = 0
	
	for line in lines[findVerticesString+1:numNodes+findVerticesString+1]:
		data = line.split()
		nodelistVector3[id] = orientation*Vector3(float(data[0])*scale,float(data[1])*scale,float(data[2])*scale)+shift
		id += 1
	
	findTriangleString=findVerticesString+numNodes
	while (lines[findTriangleString].split()[0]<>'Triangles'): # find the string with the number of Triangles
		findTriangleString+=1
	findTriangleString+=1
	numTriangles = int(lines[findTriangleString].split()[0])

	triList = []
	for i in range(numTriangles):
		triList.append([0,0,0,0])
	
	tid = 0
	for line in lines[findTriangleString+1:findTriangleString+numTriangles+1]:
		data = line.split()
		id1 = int(data[0])-1
		id2 = int(data[1])-1
		id3 = int(data[2])-1
		triList[tid][0] = tid
		triList[tid][1] = id1
		triList[tid][2] = id2
		triList[tid][3] = id3
		tid += 1
		
	nodesIds=[]; cylIds=[]; pfIds=[]
	for i in triList:
		a=nodelistVector3[i[1]]
		b=nodelistVector3[i[2]]
		c=nodelistVector3[i[3]]
		#print 'i',i
		#print 'a',a
		#print 'b',b
		#print 'c',c
		try:
			pfacetCreator1([a,b,c],radius=radius,nodesIds=nodesIds,cylIds=cylIds,pfIds=pfIds,
														wire=wire,fixed=fixed,materialNodes=materialNodes,material=material,color=color)
		except:
			pass
	return nodesIds,cylIds,pfIds

