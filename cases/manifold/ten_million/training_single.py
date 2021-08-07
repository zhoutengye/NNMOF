import keras
import tensorflow as tf
from numpy.random import seed
seed(5)
import tensorflow
tensorflow.random.set_seed(2)
import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import explained_variance_score
from sklearn.model_selection import GridSearchCV
from keras import backend as K
from keras.wrappers.scikit_learn import KerasRegressor
import json


def coeff_determination(y_true, y_pred):
    from keras import backend as K
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

def model_test():
    model = Sequential()
    model.add(Dense(units=25,
                activation=type_activation))
    model.add(Dense(units=20,
                activation=type_activation))
    model.add(Dense(units=15,
                activation=type_activation))
    model.add(Dense(units=10,
                activation=type_activation))
    model.add(Dense(units=5,
                activation=type_activation))
#    model.add(Dense(units=25,
#                activation=type_activation))
#    model.add(Dense(units=25,
#                activation=type_activation))
#    model.add(Dense(units=20,
#                activation=type_activation))
#    model.add(Dense(units=15,
#                activation=type_activation))
#    model.add(Dense(units=10,
#                activation=type_activation))
#    model.add(Dense(units=5,
#                activation=type_activation))
#    model.add(Dense(units=20,
#                activation=type_activation))
#    model.add(Dense(units=15,
#                activation=type_activation))
#    model.add(Dense(units=10,
#                activation=type_activation))
#    model.add(Dense(units=5,
#                activation=type_activation))
#    model.add(Dense(units=10,
#                activation=type_activation))
#    model.add(Dense(units=5,
#                activation=type_activation))
#    model.add(Dense(units=50,
#                activation=type_activation))
#    model.add(Dense(units=45,
#                activation=type_activation))
#    model.add(Dense(units=40,
#                activation=type_activation))
#    model.add(Dense(units=35,
#                activation=type_activation))
#    model.add(Dense(units=30,
#                activation=type_activation))
#    model.add(Dense(units=25,
#                activation=type_activation))
#    model.add(Dense(units=20,
#                activation=type_activation))
#    model.add(Dense(units=15,
#                activation=type_activation))
#    model.add(Dense(units=10,
#                activation=type_activation))
    model.add(Dense(units=2,
                activation='linear'))
    adam = keras.optimizers.Adam(lr=0.005, beta_1=0.85, beta_2=0.99,
            epsilon=None, decay=0.0, amsgrad=False)
    model.compile(loss='mean_squared_error', optimizer=adam,  metrics=[coeff_determination])
    return model

data_dir = '.'

exact_f = np.load(data_dir+'/exact_f.npy')
exact_angle = np.load(data_dir+'/exact_angle.npy')
exact_centroid = np.load(data_dir+'/exact_centroid.npy')
delta_angle = np.load(data_dir+'/delta_angle.npy')
initial_angle = np.load(data_dir+'/initial_angle.npy')

exact_centroid = exact_centroid - 0.5

exact_f = exact_f.reshape([len(exact_f),1])
# inputs = np.hstack((exact_centroid,exact_f))
inputs = np.hstack((initial_angle,exact_f))
outputs = delta_angle.copy()

# parameters
num_epoch = 200
n_batch_size = 10000
type_activation = 'relu'
type_optimizer = 'adam'

model = model_test()
# history = model.fit(inputs,outputs,validation_split=0.1,epochs=num_epoch,batch_size=n_batch_size,shuffle=False)
history = model.fit(inputs,outputs,validation_split=0.1,epochs=num_epoch,batch_size=n_batch_size)
model.save('test.h5')
with open('history.json', 'w') as f:
    json.dump(history.history,f)

print(model.predict(inputs[:10,:]))
print(outputs[:10,:])
