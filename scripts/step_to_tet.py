# from https://gitlab.com/franck.ledoux/mambo/-/blob/master/Scripts/gen_tet_mesh.py

import gmsh
import sys

model = gmsh.model
factory = model.occ
    
def check_param():
    #===========================================================================================
    # Parsing of the input parameters
    #===========================================================================================
    if(len(sys.argv)<3 or len(sys.argv)>4):
        print("Wrong usage, it should be:")
        print("\t python gen_tet_mesh.py  in.step out.vtk [mesh_size]")
        print("where:")
        print("\t - in.step is the input geometry file (step format),")
        print("\t - out.vtk contains the generated tet mesh (vtk format),")
        print("\t - mesh_size is the expected element size factor in ]0,1].")
        exit(1)
    #==========================================================================================
    # Check step file suffix
    #===========================================================================================
    
    step_file = sys.argv[1]
    #if(not step_file.endswith(".step") and not step_file.endswith(".stp")):
        #print("ERROR: the input geometry file must be at the step format (.stp or .step)")
        #exit(1)
            
    #===========================================================================================
    # Check vtk file suffix
    #===========================================================================================
    
    mesh_file = sys.argv[2]
    #if(not mesh_file.endswith(".vtk")):
        #print("ERROR: the output mesh file must be a vtk legacy file (.vtk)")
        #exit(1)
    #===========================================================================================
    # Check mesh size value
    #===========================================================================================
    mesh_size = 0.2
    if(len(sys.argv)==4):
        mesh_size = float(sys.argv[3])
        if(mesh_size<=0 or mesh_size>1):
            print("ERROR: the mesh size must be in ]0,1]")
            exit(1)
            
    params = [step_file, mesh_file, mesh_size]
    return params

def process(step_file, mesh_file, mesh_size):
    #======================================================================================
    # Conversion process
    #======================================================================================
    print("> Geometry import of "+step_file)
    model.add("step_to_tet")
    factory.importShapes(step_file)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", mesh_size)
    print("> Tet mesh generation with element size factor "+str(mesh_size))
    
    factory.synchronize()
    model.mesh.generate(3)
    print("> Writing mesh into file "+mesh_file)
    gmsh.write(mesh_file)


if __name__=="__main__":
    gmsh.initialize(sys.argv)
    params = check_param()
    process(params[0],params[1],params[2])
    gmsh.finalize()
