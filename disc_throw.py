import numpy as np
from stl import mesh
from scipy.integrate import solve_ivp
from numpy import cos, sin
import matplotlib.pyplot as plt
from os import listdir
from scipy.interpolate import griddata

# generates rotated meshes
def generate_angles(mesh):
    for theta in [-1, -0.5, -1 / 6, -1 / 12, 1 / 12]:
        mesh.rotate([0, 1, 0], -np.pi * theta / 2)
        mesh.save(f'jade\jade{theta * 90}.stl')
        mesh.rotate([0, 1, 0], np.pi * theta / 2)


# force function
def Fdrag(vx, vy, vz):
    # dummy force
    coef = 0.1
    v = np.array([vx, vy, vz])
    return -coef * v ** 2


# fuction to stop simulation when disk is at y = 0
def fallEarth(t, nezVec):
    x, y, z, vx, vy, vz, om1, om2, om3, dom1, dom2, dom3 = nezVec
    return z


# returns angle between two vectors
def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


# translate ODE parameters to CFD parameters
def get_forces(vx, vy, vz, phi, theta, om3):
    omega = om3
    dir = np.array([vx, vy, vz])
    conversion_matrix = np.array([[cos(phi), sin(phi), 0],
                                  [-sin(phi) * cos(theta), cos(phi) * cos(theta), sin(theta)],
                                  [sin(phi) * sin(theta), -sin(theta) * cos(phi), cos(theta)]])
    normal = np.matmul(conversion_matrix, np.array([0, 0, 1]))
    angle = angle_between(dir, normal) - (np.pi / 2)
    return get_model(omega, angle, dir)


#takes values from CFD model
def get_model(omega, angle, velocity):
    F = Fdrag(velocity[0], velocity[1], velocity[2])
    M = np.array([0, 0, 0])
    return F, M


def model(t, nezVec):
    # unpack unknowns
    x, y, z, vx, vy, vz, phi, theta, psi, om1, om2, om3 = nezVec

    # calculate forces and acceleration
    F, M = get_forces(vx, vy, vz, phi, theta, om3)
    a = F / m

    # RHS of eq.
    # change of position
    dxdt = vx
    dydt = vy
    dzdt = vz
    dvxdt = a[0]
    dvydt = a[1]
    dvzdt = a[2] - g

    # change of angular velocity
    dom1 = (1 / Ixy) * ((Ixy - Iz) * om2 * om3 + M[0])
    dom2 = (1 / Ixy) * ((Iz - Ixy) * om1 * om3 + M[1])
    dom3 = M[2] / Iz
    # change of rotation
    dphi = 1 / cos(theta) * (cos(phi) * sin(theta) * om1 + sin(phi) * sin(theta) * om2 + cos(theta) * om3)
    dtheta = 1 / cos(theta) * (- sin(phi) * cos(theta) * om1 + cos(phi) * sin(theta) * om2)
    dpsi = 1 / cos(theta) * (sin(phi) * om2 + cos(phi) * om1)
    return dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt, dphi, dtheta, dpsi, dom1, dom2, dom3


# parameters
v0 = 10  # initial velocity
angl = 15  # angle of throw
init_rotation = 10;
g = 9.81

# getting information from stl file
disc = mesh.Mesh.from_file('jade.stl')
volume, m, cog, inertia = disc.get_mass_properties_with_density(1250)
Ixy = inertia[0, 0]
Iz = inertia[2, 2]

# load data
CFD_data = dict()
for file in listdir("forcesDir"):
    CFD_data[file.split('.')[0]] = np.load("forcesDir/"+file)


# initial conditions
x, y, z = 0, 0, 0  # m  -- positions
vx, vy, vz = v0 * np.cos(angl * np.pi / 180), 0,v0 * np.sin(angl * np.pi / 180)  # ms -- velocities
om1, om2, om3 = 0, 0, init_rotation
phi, theta, psi = 0, angl * np.pi / 180, 0
init_cond = x, y, z, vx, vy, vz, phi, theta, psi, om1, om2, om3

t = np.linspace(0, 1000, 2)  # s -- time
fallEarth.direction = -1
fallEarth.terminal = True
solution = solve_ivp(model, t, init_cond, events=fallEarth, method='BDF')
plt.plot(solution.y[0, :], solution.y[2, :])
plt.plot(solution.y[0, :], solution.y[7,:])
plt.show()
