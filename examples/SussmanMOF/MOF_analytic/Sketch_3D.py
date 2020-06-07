import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

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

    # edges2 = [plane_definition]
    edges2 = plane_definition

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


cube_definition = [
    (0,0,0), (0,1,0), (1,0,0), (0,0,1)
]

#  tetrahedrons
Tri_intersection = []
Tri_intersection.append([[0,0.5,0], [0.5,1,0], [0,1,0.5]])
Tri_intersection.append([[0.5,0,0], [0,0.5,0], [0,0,0.5]])
Tri_intersection.append([[0.5,0,0], [1,0.5,0], [1,0,0.5]])
Tri_intersection.append([[1,0.5,0], [0.5,1,0], [1,1,0.5]])
Tri_intersection.append([[0,0.5,1], [0.5,0,1], [0,0,0.5]])
Tri_intersection.append([[1,0.5,1], [0.5,0,1], [1,0,0.5]])
Tri_intersection.append([[1,1,0.5], [1,0.5,1], [0.5,1,1]])
Tri_intersection.append([[0,0.5,1], [0.5,1,1], [0,1,0.5]])

#  quad
Quad1_intersection = []
Quad1_intersection.append([ [0,0.3,0], [0.3,0,0], [0.7,0,1], [0,0.7,1] ])
Quad1_intersection.append([ [0.3,0,0], [1,0.7,0], [1,0.3,1], [0.7,0,1] ])
Quad1_intersection.append([ [1,0.7,0], [0.7,1,0], [0.3,1,1], [1,0.3,1] ])
Quad1_intersection.append([ [0,0.3,0], [0.7,1,0], [0.3,1,1], [0,0.7,1] ])

Quad1_intersection.append([ [0.3,0,0], [0,0,0.3], [0,1,0.7], [0.7,1,0] ])
Quad1_intersection.append([ [0,0,0.3], [0.7,0,1], [0.3,1,1], [0,1,0.7] ])
Quad1_intersection.append([ [0.7,0,1], [1,0,0.7], [1,1,0.3], [0.3,1,1] ])
Quad1_intersection.append([ [1,0,0.7], [0.3,0,0], [0.7,1,0], [1,1,0.3] ])

Quad1_intersection.append([ [0,0,0.3], [0,0.3,0], [1,0.7,0], [1,0,0.7] ])
Quad1_intersection.append([ [0,0.3,0], [0,1,0.7], [1,1,0.3], [1,0.7,0] ])
Quad1_intersection.append([ [0,1,0.7], [0,0.7,1], [1,0.3,1], [1,1,0.3] ])
Quad1_intersection.append([ [0,0.7,1], [0,0,0.3], [1,0,0.7], [1,0.3,1] ])

#  quad2
Quad2_intersection = []
Quad2_intersection.append([ [0.2,0,0], [0.2,0,1], [0.5,1,1], [0.5,1,0] ])
Quad2_intersection.append([ [0,0.5,0], [0,0.5,1], [1,0.2,1], [1,0.2,0] ])
Quad2_intersection.append([ [0.5,0,0], [0.5,0,1], [0.8,1,1], [0.8,1,0] ])
Quad2_intersection.append([ [0,0.8,0], [0,0.8,1], [1,0.5,1], [1,0.5,0] ])
Quad2_intersection.append([ [0,0.2,0], [1,0.2,0], [1,0.5,1], [0,0.5,1] ])
Quad2_intersection.append([ [0,0,0.5], [1,0,0.5], [1,1,0.2], [0,1,0.2] ])
Quad2_intersection.append([ [0,0.5,0], [1,0.5,0], [1,0.8,1], [0,0.8,1] ])
Quad2_intersection.append([ [0,0,0.8], [1,0,0.8], [1,1,0.5], [0,1,0.5] ])
Quad2_intersection.append([ [0,0,0.2], [0,1,0.2], [1,1,0.5], [1,0,0.5] ])
Quad2_intersection.append([ [0.5,0,0], [0.5,1,0], [0.2,1,1], [0.2,0,1] ])
Quad2_intersection.append([ [0,0,0.5], [0,1,0.5], [1,1,0.8], [1,0,0.8] ])
Quad2_intersection.append([ [0.8,0,0], [0.8,1,0], [0.5,1,1], [0.5,0,1] ])



# Pentagon
Panta_intersection = []
Panta_intersection.append([ [0,0.7,0], [0.7,0,0], [1,0,0.2], [1,1,13.0/15.0], [0,1,0.2]])
Panta_intersection.append([ [0.3,0,0], [1,0.7,0], [1,1,0.2], [0,1,13.0/15.0], [0,0,0.2]])
Panta_intersection.append([ [1,0.3,0], [0.3,1,0], [0,1,0.2], [0,0,13.0/15.0], [1,0,0.2]])
Panta_intersection.append([ [0.7,1,0], [0,0.3,0], [0,0,0.2], [1,0,13.0/15.0], [1,1,0.2]])
Panta_intersection.append([ [0,0.7,1], [0.7,0,1], [1,0,0.8], [1,1,1-13.0/15.0], [0,1,0.8]])
Panta_intersection.append([ [0.3,0,1], [1,0.7,1], [1,1,0.8], [0,1,1-13.0/15.0], [0,0,0.8]])
Panta_intersection.append([ [1,0.3,1], [0.3,1,1], [0,1,0.8], [0,0,1-13.0/15.0], [1,0,0.8]])
Panta_intersection.append([ [0.7,1,1], [0,0.3,1], [0,0,0.8], [1,0,1-13.0/15.0], [1,1,0.8]])

Panta_intersection.append([ [0,0,0.7], [0,0.7,0], [0.2,1,0], [13.0/15.0,1,1], [0.2,0,1]])
Panta_intersection.append([ [0,0.3,0], [0,1,0.7], [0.2,1,1], [13.0/15.0,0,1], [0.2,0,0]])
Panta_intersection.append([ [0,1,0.3], [0,0.3,1], [0.2,0,1], [13.0/15.0,0,0], [0.2,1,0]])
Panta_intersection.append([ [0,0.7,1], [0,0,0.3], [0.2,0,0], [13.0/15.0,1,0], [0.2,1,1]])
Panta_intersection.append([ [1,0,0.7], [1,0.7,0], [0.8,1,0], [1-13.0/15.0,1,1], [0.8,0,1]])
Panta_intersection.append([ [1,0.3,0], [1,1,0.7], [0.8,1,1], [1-13.0/15.0,0,1], [0.8,0,0]])
Panta_intersection.append([ [1,1,0.3], [1,0.3,1], [0.8,0,1], [1-13.0/15.0,0,0], [0.8,1,0]])
Panta_intersection.append([ [1,0.7,1], [1,0,0.3], [0.8,0,0], [1-13.0/15.0,1,0], [0.8,1,1]])

Panta_intersection.append([ [0.7,0,0], [0,0,0.7], [0,0.2,1], [1,13.0/15.0,1], [1,0.2,0]])
Panta_intersection.append([ [0,0,0.3], [0.7,0,1], [1,0.2,1], [1,13.0/15.0,0], [0,0.2,0]])
Panta_intersection.append([ [0.3,0,1], [1,0,0.3], [1,0.2,0], [0,13.0/15.0,0], [0,0.2,1]])
Panta_intersection.append([ [1,0,0.7], [0.3,0,0], [0,0.2,0], [0,13.0/15.0,1], [1,0.2,1]])
Panta_intersection.append([ [0.7,1,0], [0,1,0.7], [0,0.8,1], [1,1-13.0/15.0,1], [1,0.8,0]])
Panta_intersection.append([ [0,1,0.3], [0.7,1,1], [1,0.8,1], [1,1-13.0/15.0,0], [0,0.8,0]])
Panta_intersection.append([ [0.3,1,1], [1,1,0.3], [1,0.8,0], [0,1-13.0/15.0,0], [0,0.8,1]])
Panta_intersection.append([ [1,1,0.7], [0.3,1,0], [0,0.8,0], [0,1-13.0/15.0,1], [1,0.8,1]])


# Hexagon
Hex_intersection = []
Hex_intersection.append([[0.7,0,0], [1,0,0.3], [1,0.7,1],[0.7,1,1],[0,1,0.3],[0,0.7,0]])
Hex_intersection.append([[0.3,0,0], [1,0.7,0], [1,1,0.3],[0.3,1,1],[0,0.7,1],[0,0,0.3]])
Hex_intersection.append([[0.3,1,0], [0,1,0.3], [0,0.3,1],[0.3,0,1],[1,0,0.3],[1,0.3,0]])
Hex_intersection.append([[0.7,1,0], [0,0.3,0], [0,0,0.3],[0.7,0,1],[1,0.3,1],[1,1,0.3]])

Hex_intersection.append([[0.3,0,0], [1,0,0.7], [1,0.3,1],[0.3,1,1],[0,1,0.7],[0,0.3,0]])
Hex_intersection.append([[0.7,0,0], [1,0.3,0], [1,1,0.7],[0.7,1,1],[0,0.3,1],[0,0,0.7]])
Hex_intersection.append([[0.7,1,0], [0,1,0.7], [0,0.7,1],[0.7,0,1],[1,0,0.7],[1,0.7,0]])
Hex_intersection.append([[0.3,1,0], [0,0.7,0], [0,0,0.7],[0.3,0,1],[1,0.7,1],[1,1,0.7]])

# Trainsision
Trans_intersection = []
Trans_intersection.append([ [0,0,0], [0,1,1], [1,1,0] ])
Trans_intersection.append([ [0,0,0], [0,0,1], [1,1,1], [1,1,0] ])
Trans_intersection.append([ [0,0,0.5], [1,0,0], [1,1,0.5], [0,1,1]])

kk = 24

plot_cube_plane(cube_definition,Panta_intersection[kk-1:kk])
plt.show()