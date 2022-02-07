# OBJ -> STL
# Use this script from Blender>Scripting
# don't forget to change input path


import bpy
import os

path = '/data/'  # set this path

for root, dirs, files in os.walk(path):
    for f in files:
        print(f)
        if f.endswith('.obj') :
            mesh_file = os.path.join(path, f)
            obj_file = os.path.splitext(mesh_file)[0] + ".stl"

            bpy.ops.object.select_all(action='SELECT')
            bpy.ops.object.delete()

            bpy.ops.import_scene.obj(filepath=mesh_file) # change this line

            bpy.ops.object.select_all(action='SELECT')

            bpy.ops.export_mesh.stl(filepath=obj_file)