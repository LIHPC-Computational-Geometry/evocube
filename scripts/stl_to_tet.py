# surface STL -> Tet .mesh

# https://gitlab.onelab.info/gmsh/gmsh/-/issues/1598

import gmsh
import sys

gmsh.initialize()

stl = gmsh.merge(sys.argv[1]) #3D STL file of a cylinder
s = gmsh.model.getEntities(2)
surf = gmsh.model.geo.addSurfaceLoop([e[1] for e in s])

gmsh.model.geo.add_volume([surf])
gmsh.model.geo.synchronize()


#gmsh.option.setNumber('Mesh.Algorithm', 1)
#gmsh.option.setNumber('Mesh.MeshSizeMax', 1.0)
#gmsh.option.setNumber('Mesh.MeshSizeMin', 1.0)
gmsh.model.mesh.generate(3)
gmsh.write(sys.argv[2])