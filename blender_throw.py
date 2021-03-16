import bpy
import numpy
import sys
import subprocess
import os

sys.path.append("C:/Users/hruza/Documents/Programování/Disc throw")
os.chdir("C:/Users/hruza/Documents/Programování/Disc throw")
import disc_throw

v0 = 10  # initial velocity
angl = 15  # angle of throw
init_rotation = 10;

jade = bpy.data.objects["Jade"]
solution = disc_throw.compute(v0,angl,init_rotation)

# apply solution to jade