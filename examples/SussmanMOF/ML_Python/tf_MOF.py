import math

from IPython import display
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
import tensorflow as tf
from tensorflow.python.data import Dataset

tf.logging.set_verbosity(tf.logging.ERROR)
pd.options.display.max_rows = 10
pd.options.display.float_format = '{:.1f}'.format

# Load Training Date
f = np.loadtxt('MOF_Training.dat',unpack='true')
Norm_x = f[0,:]
Norm_y = f[1,:]
Intersecpt = f[2,:]
Area = f[3,:]
Centroid_x = f[4,:]
Centroid_y = f[5,:]

MOF_Training = pd.DataFrame({ 'Norm X': Norm_x, 'Norm Y': Norm_y, 'Intersecpt': Intersecpt, 'Area': Area, 'Centroid X': Centroid_x, 'Centroid Y': Centroid_y})

N_OUTPUTS = 3
N_INPUTS = 3

N_HIDDEN_UNITS =  
N_EPOCHS =  

input = tf.placeholder(tf.float32, shape=[None, N_INPUTS], name='input')  # input here

outputs = tf.placeholder(tf.float32, shape=[None, N_OUTPUTS], name='output')  # one sample is something like[Ax,Ay,Az]

# one hidden layer with 3 outputs
W = {
    'hidden': tf.Variable(tf.random_normal([N_INPUTS, N_HIDDEN_UNITS])),
    'output': tf.Variable(tf.random_normal([N_HIDDEN_UNITS, N_OUTPUTS]))
}
biases = {
    'hidden': tf.Variable(tf.random_normal([N_HIDDEN_UNITS], mean=1.0)),
    'output': tf.Variable(tf.random_normal([N_OUTPUTS]), mean=1.0)
}