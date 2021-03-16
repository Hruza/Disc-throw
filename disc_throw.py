import numpy as np
from stl import mesh
from scipy.integrate import solve_ivp
from numpy import cos, sin
import matplotlib.pyplot as plt
from os import listdir
from scipy.interpolate import griddata
from scipy.spatial.transform import Rotation as rot

g = 9.81


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
    return min(z,1-t)


# returns angle between two vectors
def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


# translate ODE parameters to CFD parameters
def get_forces(vx, vy, vz, phi, theta, om3, points, data_F, data_M):
    omega = om3
    dir = np.array([vx, vy, vz])

    rotation = rot.from_rotvec(np.array([0,0,1])*phi)
    rotation = rot.from_rotvec(rotation.apply([1,0,0])*theta)
    normal = rotation.apply(np.array([0, 0, 1]))

    angle = angle_between(dir, normal) - (np.pi / 2)
    F, M = get_model(omega, angle, dir, points, data_F, data_M)
    M[0] *= -1
    localX = (dir / np.linalg.norm(dir))
    localY = np.cross( normal,localX)
    localZ = np.cross( localX,localY)
    F = -(F[0] * localX) + (F[1] * localY) + (F[2] * localZ)
    return F, M


# takes values from CFD model
def get_model(omega, angle, velocity, points, data_F, data_M):
    U = np.linalg.norm(velocity)
    F = np.array([griddata(points, f, [U, omega, angle]) for f in data_F])
    M = np.array([griddata(points, m, [U, omega, angle]) for m in data_M])
    return F, M


def model(t, nezVec, m, Ixy, Iz, points, data_F, data_M):
    # unpack unknowns
    x, y, z, vx, vy, vz, phi, theta, psi, om1, om2, om3 = nezVec

    # calculate forces and acceleration
    F, M = get_forces(vx, vy, vz, phi, theta, om3, points, data_F, data_M)
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


def compute(v0, angle, init_rotation):
    # getting information from stl file
    disc = mesh.Mesh.from_file('jade.stl')
    volume, m, cog, inertia = disc.get_mass_properties_with_density(1250)
    Ixy = inertia[0, 0]
    Iz = inertia[2, 2]

    # load data
    CFD_data = dict()
    for file in listdir("forcesDir"):
        CFD_data[file.split('.')[0]] = np.load("forcesDir/" + file)

    points = np.array([[U, omg, al] for U in CFD_data["U"] for omg in CFD_data["omg"] for al in CFD_data["al"]])
    data_F = [CFD_data["Fx"].flatten(), CFD_data["Fy"].flatten(), CFD_data["Fz"].flatten()]
    data_M = [CFD_data["Mx"].flatten(), CFD_data["My"].flatten(), CFD_data["Mz"].flatten()]

    # initial conditions
    x, y, z = 0, 0, 0  # m  -- positions
    vx, vy, vz = v0 * np.cos(angl * np.pi / 180), 0, v0 * np.sin(angl * np.pi / 180)  # ms -- velocities
    om1, om2, om3 = 0, 0, init_rotation
    phi, theta, psi = -np.pi/2, angle * np.pi / 180, 0
    init_cond = x, y, z, vx, vy, vz, phi, theta, psi, om1, om2, om3

    t = np.linspace(0, 1000, 2)  # s -- time
    fallEarth.direction = -1
    fallEarth.terminal = True

    odr_model = lambda t, nezVec: model(t, nezVec, m, Ixy, Iz, points, data_F, data_M)

    solution = solve_ivp(odr_model, t, init_cond, events=fallEarth, method='BDF')
    return solution


# parameters
v0 = 10  # initial velocity
angl = 15  # angle of throw
init_rotation = 10;

solution = compute(v0, angl, init_rotation)

plt.plot(solution.y[0, :], solution.y[2, :])
#plt.plot(solution.y[0, :], solution.y[1, :])
plt.show()
