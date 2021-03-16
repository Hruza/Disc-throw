import bpy
import numpy
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
    obj.rotation_euler.rotate_axis("Z",-fi)
    obj.rotation_euler.rotate_axis("Y",-theta)
    obj.rotation_euler.rotate_axis("Z",-psi)

v0 = 15  # initial velocity
angl = 15  # angle of throw
init_rotation = 10;

disc = bpy.data.objects["Jade"]
solution = disc_throw.compute(v0,angl,init_rotation)


frame_number=0
# apply solution to jadefor i in range(50):
for i in range(len(solution.t)):
    bpy.context.scene.frame_set(frame_number)
    disc.location = (solution.y[0,i],solution.y[1,i],solution.y[2,i])
    rotate_euler(disc,solution.y[6,i],solution.y[7,i],solution.y[8,i])
    disc.keyframe_insert(data_path="location",index=-1)
    disc.keyframe_insert(data_path="rotation_euler",index=-1)
    frame_number+=1