import numpy as np
import tensorflow as tf
from keras.models import Sequential
from keras.layers import Dense

data_dir = 'uniform'

exact_f = tf.convert_to_tensor(np.load(data_dir+'/exact_f.npy'))
exact_angle = tf.convert_to_tensor(np.load(data_dir+'/exact_angle.npy'))
exact_centroid = tf.convert_to_tensor(np.load(data_dir+'/exact_centroid.npy'))
delta_angle = tf.convert_to_tensor(np.load(data_dir+'/delta_angle.npy'))
initial_angle = tf.convert_to_tensor(np.load(data_dir+'/initial_angle.npy'))

# Define Sequential model with 3 layers
model = Sequential()
model.add(Dense(units=100,
                input_dim=4,
                activation='relu'))

# Call model on a test input