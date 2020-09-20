import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import explained_variance_score
from keras import backend as K

def coeff_determination(y_true, y_pred):
    from keras import backend as K
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

data_dir = 'uniform'

exact_f = np.load(data_dir+'/exact_f.npy')
exact_angle = np.load(data_dir+'/exact_angle.npy')
exact_centroid = np.load(data_dir+'/exact_centroid.npy')
delta_angle = np.load(data_dir+'/delta_angle.npy')
initial_angle = np.load(data_dir+'/initial_angle.npy')

exact_centroid = exact_centroid - 0.5

exact_f = exact_f.reshape([len(exact_f),1])
inputs = np.hstack((exact_centroid,exact_f))
outputs = delta_angle.copy()

# parameters
num_epoch = 20
n_batch_size = 1000
type_activation = 'relu'
type_optimizer = 'adam'

# Nueral network with keras
nnkeras = Sequential()
nnkeras.add(Dense(units=20,
                input_dim=4,
                activation=type_activation))
nnkeras.add(Dense(units=10,
                activation=type_activation))
nnkeras.add(Dense(units=2,
                activation='linear'))

nnkeras.compile(loss='mean_squared_error', optimizer=type_optimizer,  metrics=[coeff_determination])

# sklearn nueral network, random forest and decision tree
nn = MLPRegressor(hidden_layer_sizes=(20,),verbose=True)
rf = RandomForestRegressor(n_estimators=10,verbose=True,max_depth=20)
dt = DecisionTreeRegressor(max_depth=20)

sk = nn

nnkeras.fit(inputs,outputs,epochs=num_epoch,batch_size=n_batch_size)
sk.fit(inputs,outputs)

# nn_keras_pred = nnkeras.predict(inputs)
# nn_sklearn_pred = sk.predict(inputs)
print(sk.score(inputs,outputs))


vali = slice(100,110)
input = inputs[vali,:]
output = outputs[vali,:]

print(output)
print(sk.predict(input))
print(nnkeras.predict(input))