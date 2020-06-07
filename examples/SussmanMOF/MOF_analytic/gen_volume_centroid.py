from calculate_volume_and_centroid import calculate_volume_and_centroid
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import scipy
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
import math
import mayavi
from mayavi.mlab import *

n1 = np.sqrt(1.0/3.0)
n2 = np.sqrt(1.0/3.0)
n3 = np.sqrt(1-n1*n1-n2*n2)
d  = 0.0
nn = [n1,n2,n3]

cube = [ [0.5,0.5,0.5],[0.5,0.5,-0.5], [0.5,-0.5,0.5], [0.5,- 0.5,-0.5],
         [-0.5,0.5,0.5],[-0.5,0.5,-0.5], [-0.5,-0.5,0.5], [-0.5,- 0.5,-0.5],
]

vectors = np.random.normal(size=(10,3))
distances = np.linspace(-np.sqrt(2.0/3.0),np.sqrt(2.0/3.0),10)

centroids = []
volume = []

for i in range(len(vectors)):
	for j in range(len(distances)):
		nn = [vectors[i,0], vectors[i,1], vectors[i,2] ]
		d = distances[j]
		volume_centroid = calculate_volume_and_centroid(cube,nn,d)
		if np.any(volume_centroid) != None:
			centroids = np.append(centroids,volume_centroid[0:3],axis=0)
			volume.append(volume_centroid[3])

centroids = np.resize(centroids, (len(volume),3))
volume = np.array(volume)

mayavi.mlab.points3d(centroids[0],centroids[1],centroids[2])