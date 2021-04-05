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
    obj.location = (solution[0, i], solution[1, i], solution[2, i])
    rotate_vector(obj,direction)
    obj.scale[2]=np.linalg.norm(direction)
    obj.keyframe_insert(data_path="location", index=-1)
    obj.keyframe_insert(data_path="rotation_euler", index=-1)
    obj.keyframe_insert(data_path="scale", index=-1)
    

solution = np.load("solutions/30_1_-30.npy")
disc = bpy.data.objects["Jade"]
arrow = [bpy.data.objects["Arrow"+str(i)] for i in range(6) ]

points, data_F, data_M = disc_throw.load_data()

data = np.array([disc_throw.get_forces(solution[3,i],solution[4,i],solution[5,i],solution[6,i],solution[7,i],solution[8,i],solution[11,i],points,data_F,data_M,debug=True) for i in range(len(solution[0,:])) ])

corot_to_global=np.array([np.transpose( np.array([data[j,i] for i in [6,7,8]]) ) for j in range(len(solution[0])) ])
local_to_global=np.array([np.transpose( np.array([data[j,i] for i in [3,4,5]]) ) for j in range(len(solution[0])) ])

frame_number=0
# apply solution to jadefor i in range(50):
for i in range(len(solution[0])):
    bpy.context.scene.frame_set(frame_number)
    disc.location = (solution[0,i],solution[1,i],solution[2,i])
    rotate_euler(disc,solution[6,i],solution[7,i],solution[8,i])
    disc.keyframe_insert(data_path="location",index=-1)
    disc.keyframe_insert(data_path="rotation_euler",index=-1)
    
    add_vector(arrow[0],data[i,0],solution)
    add_vector(arrow[1],100*np.matmul(corot_to_global[i],data[i,1]),solution)
    
    add_vector(arrow[2],data[i,3],solution)
    add_vector(arrow[3],data[i,4],solution)
    add_vector(arrow[4],data[i,5],solution)
    
    add_vector(arrow[5],np.array([0,0,-9.81]),solution)
    frame_number+=1