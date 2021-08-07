h5_file = 'test.h5'
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
from sklearn.externals import joblib
from keras.models import load_model

def keras2f(weights_c,activation_type):
    n_layers = int(len(weights_c)/2)
    intercepts = []
    weights = []
    layer_sizes = []
    layer_sizes.append(np.array(weights_c[0]).shape[0])
    for i in range(n_layers):
        layer_sizes.append(np.array(weights_c[int(i*2)]).shape[1])
        intercepts.append(weights_c[int(i*2+1)])
        weights.append(weights_c[int(i*2)])
    file_name = 'nn_coef.dat'
    nn_output  = open(file_name, 'w')
    nn_output.write('! n_layers\n%d\n'%(n_layers+1))
    nn_output.write('! layer_sizes\n')
    for i in range(n_layers+1):
        nn_output.write('%d '%layer_sizes[i])
    nn_output.write('\n')
    nn_output.write('! intercepts\n')
    for i in range(len(intercepts)):
        for j in range(len(intercepts[i])):
            nn_output.write('%f '%intercepts[i][j])
        nn_output.write('\n')

    for i in range(len(weights)):
        nn_output.write('!! layer%d\n'%i)
        coef = np.transpose(weights[i])
        # coef = weights[i]
        for j in range(len(coef)):
            for k in range(len(coef[j])):
                nn_output.write('%f '%coef[j][k])
            nn_output.write('\n')
    nn_output.write('! activations\n')
    nn_output.write(activation_type)
    nn_output.write('\n')
    nn_output.write('! out_activations\n')
    nn_output.write('identity')
    nn_output.write('\n')

def test_predict(val,weights):
    v = val
    lay = int(len(weights)/2)
    for i in range(lay):
        coef = np.array(weights[int(i*2)])
        bias = np.array(weights[int(i*2)+1])
        v = np.matmul(v,coef) + bias
        if i != lay-1:
            v[v<0] = 0
    print(v)

def coeff_determination(y_true, y_pred):
    from keras import backend as K
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

data_dir = 'uniform'
activation_type = 'relu'

# exact_f = np.load(data_dir+'/exact_f.npy')
# exact_angle = np.load(data_dir+'/exact_angle.npy')
# exact_centroid = np.load(data_dir+'/exact_centroid.npy')
# delta_angle = np.load(data_dir+'/delta_angle.npy')
# initial_angle = np.load(data_dir+'/initial_angle.npy')

# exact_centroid = exact_centroid - 0.5

# exact_f = exact_f.reshape([len(exact_f),1])
# inputs = np.hstack((exact_centroid,exact_f))
# outputs = delta_angle.copy()

model = load_model(h5_file, custom_objects={'coeff_determination':coeff_determination})

weights = model.get_weights()
# test_predict(inputs[:10,:],weights)

# print(model.predict(inputs[:10,:]))
# print(outputs[:10,:])

keras2f(weights,activation_type)
