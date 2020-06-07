import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import scipy
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
import math

def plot_cube_plane(cube_definition,plane_definition):
    cube_definition_array = [
        np.array(list(item))
        for item in cube_definition
    ]

    points = []
    points += cube_definition_array
    vectors = [
        cube_definition_array[1] - cube_definition_array[0],
        cube_definition_array[2] - cube_definition_array[0],
        cube_definition_array[3] - cube_definition_array[0]
    ]

    points += [cube_definition_array[0] + vectors[0] + vectors[1]]
    points += [cube_definition_array[0] + vectors[0] + vectors[2]]
    points += [cube_definition_array[0] + vectors[1] + vectors[2]]
    points += [cube_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

    points = np.array(points)

    edges = [
        [points[0], points[3], points[5], points[1]],
        [points[1], points[5], points[7], points[4]],
        [points[4], points[2], points[6], points[7]],
        [points[2], points[6], points[3], points[0]],
        [points[0], points[2], points[4], points[1]],
        [points[3], points[6], points[7], points[5]]
    ]

    edges2 = [plane_definition]
    # edges2 = plane_definition

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
    faces.set_facecolor((0,0,1,0.1))

    faces2 = Poly3DCollection(edges2, linewidths=1, edgecolors='k')
    faces.set_facecolor((0,0,1,0.1))

    ax.add_collection3d(faces)
    ax.add_collection3d(faces2)

    # Plot the points themselves to force the scaling of the axes
    ax.scatter(points[:,0], points[:,1], points[:,2], s=0)

    ax.set_aspect('equal')

def plane(n,distance):
	x1 = n[0]*distance
	x2 = n[1]*distance
	x3 = n[2]*distance

	return([x1,x2,x3])

def find_x(yz,n,distance):
	x = - ( yz[0] * n[1] + yz[1] * n[2] - distance ) / (n[0]+1e-12)
	return([x,yz[0],yz[1]])

def find_y(xz,n,distance):
	y = - ( xz[0] * n[0] + xz[1] * n[2] - distance ) / (n[1]+1e-12)
	return([xz[0],y,xz[1]])

def find_z(xy,n,distance):
	z = - ( xy[0] * n[0] + xy[1] * n[1] - distance ) / (n[2]+1e-12)
	return([xy[0],xy[1],z])

def sort_clockwise(p,n):
    l = len(p)
    pp = np.transpose(p)
    pc = np.array([np.mean(pp[0]), np.mean(pp[1]), np.mean(pp[2])])
    vec =  pc -p

    ang = []
    for i in range(len(p)-1):
        ang.append(angle(vec[0],vec[i+1],n))

    for i in range(0,len(p)-2):
        for j in range(i+1,len(p)-1):
            if ang[i] < ang[j]:
                a = ang[i]
                ang[i] = ang[j]
                ang[j] = a
                a = p[i+1]
                p[i+1] = p[j+1]
                p[j+1] = a

    return p

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def crossproduct(v1, v2):
  v31 = v1[1] * v2[2] - v1[2] * v2[1]
  v32 = v1[2] * v2[0] - v1[0] * v2[2]
  v33 = v1[0] * v2[1] - v1[1] * v2[0]
  return [v31,v32,v33]

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2,n):
    angle = math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
    if (dotproduct(n,crossproduct(v1,v2))) < 0.0:
        angle = np.pi * 2.0 - angle 

    return angle

def volume_and_centroid(p,cube,n,distance):
	
    pp = p
    p_pos = []
    p_neg = []
    d_pos = []
    d_neg = []
    for i in range(8):
        d = (n[0]*cube[i][0] + n[1]*cube[i][1] + n[2]*cube[i][2]-distance) / np.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
        if d > 0:
            p_pos.append(cube[i])
            d_pos.append(d)
        else:
            if d<0:
                p_neg.append(cube[i])
                d_neg.append(d)

    if np.max(np.abs(d_pos)) >= np.max(np.abs(d_neg)):
        for i in range(len(p_neg)):
            pp.append(p_neg[i])
    else:
        for i in range(len(p_pos)):
            pp.append(p_pos[i])

    pp = np.array(pp)


    dt = Delaunay(pp)

    tets = dt.points[dt.simplices]

    tetra_volume = tetrahedron_volume(tets[:, 0], tets[:, 1], 
                                tets[:, 2], tets[:, 3])

    tetra_centroid = tetrahedron_centroid(tets[:, 0], tets[:, 1], 
                                tets[:, 2], tets[:, 3])

    volume = np.sum(tetra_volume)

    centroid = np.transpose(tetra_centroid).dot(np.transpose(tetra_volume))

    return np.append(centroid,volume)

def tetrahedron_volume(a, b, c, d):

    volume = np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6.0
    centroid = np.abs(a+b+c+d) / 4.0

    return volume

def tetrahedron_centroid(a, b, c, d):

    centroid = np.abs(a+b+c+d) / 4.0

    return centroid

def calculate_volume_and_centroid(cube,nn,d):

    n1 = nn[0]
    n2 = nn[1]
    n3 = nn[2]
    n  = [n1,n2,n3]
    distance = d

    xx = [ [-0.5,-0.5], [-0.5,0.5], [0.5,-0.5], [0.5,0.5] ]
    xx = [ [-0.5,-0.5], [0.5,-0.5], [-0.5,0.5], [0.5,0.5] ]


    points = []

    for i in range(4):
        p = find_x(xx[i],n,distance)
        if p not in points:
            if p[0] >= -0.5 and p[0] <= 0.5:
                points.append(p)
        p = find_y(xx[i],n,distance)
        if p not in points:
            if p[1] >= -0.5 and p[1] <= 0.5:
                points.append(p)
        p = find_z(xx[i],n,distance)
        if p not in points:
            if p[2] >= -0.5 and p[2] <= 0.5:
                points.append(p)

    if len(points) <= 3:
        return 

    # cube_definition = [
    #     (-0.5,-0.5,-0.5), (-0.5,0.5,-0.5), (0.5,-0.5,-0.5), (-0.5,-0.5,0.5)
    # ]
    # points = sort_clockwise(points, n)
    # plot_cube_plane(cube_definition,points)
    # plt.show()



    volume_centroid = volume_and_centroid(points,cube,n,distance)

    return volume_centroid

