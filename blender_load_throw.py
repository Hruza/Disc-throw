import bpy
import numpy as np
import sys
import subprocess
import os
import mathutils
import math

sys.path.append("C:/Users/hruza/Documents/Programování/Disc throw")
os.chdir("C:/Users/hruza/Documents/Programování/Disc throw")

import disc_throw

def rotate_euler(obj, fi,theta, psi):
    obj.rotation_euler = mathutils.Euler((0,0,0))
    obj.rotation_euler.rotate_axis("Z",fi)
    obj.rotation_euler.rotate_axis("X",theta)
    obj.rotation_euler.rotate_axis("Z",psi)

def rotate_vector(obj, direction):
    direction=mathutils.Vector(direction)
    loc_obj = obj.matrix_world.to_translation()
    # point the cameras '-Z' and use its 'Y' as up
    rot_quat = (-direction).to_track_quat('-Z', 'Y')
    # assume we're using euler rotation
    obj.rotation_euler = rot_quat.to_euler()
    
def add_vector(obj, direction, solution):
    arrow.location = (solution[0, i], solution[1, i], solution[2, i])
    rotate_vector(arrow,direction)
    arrow.scale[2]=np.linalg.norm(direction)
    arrow.keyframe_insert(data_path="location", index=-1)
    arrow.keyframe_insert(data_path="rotation_euler", index=-1)
    arrow.keyframe_insert(data_path="scale", index=-1)
    

solution = np.load("solutions/25_7_-30.npy")
disc = bpy.data.objects["Jade"]
arrow = bpy.data.objects["Arrow"]
arrow2 = bpy.data.objects["Arrow2"]

points, data_F, data_M = disc_throw.load_data()

forces = np.array([disc_throw.get_forces(solution[3,i],solution[4,i],solution[5,i],solution[6,i],solution[7,i],solution[8,i],solution[11,i],points,data_F,data_M)[0] for i in range(len(solution[0,:])) ])
moments=np.array([disc_throw.get_forces(solution[3,i],solution[4,i],solution[5,i],solution[6,i],solution[7,i],solution[8,i],solution[11,i],points,data_F,data_M)[1] for i in range(len(solution[0,:])) ])


frame_number=0
# apply solution to jadefor i in range(50):
for i in range(len(solution[0])):
    bpy.context.scene.frame_set(frame_number)
    disc.location = (solution[0,i],solution[1,i],solution[2,i])
    rotate_euler(disc,solution[6,i],solution[7,i],solution[8,i])
    disc.keyframe_insert(data_path="location",index=-1)
    disc.keyframe_insert(data_path="rotation_euler",index=-1)
    
    add_vector(arrow,forces[i],solution)
    
    #add moment viasualisation 
    
    
    frame_number+=1