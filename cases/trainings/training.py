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

def coeff_determination(y_true, y_pred):
    from keras import backend as K
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

def create_model(neurons=1):
    model = Sequential()
    model.add(Dense(units=neurons,
                input_dim=4,
                activation=type_activation))
    model.add(Dense(units=10,
                activation=type_activation))
    model.add(Dense(units=2,
                activation='linear'))
    model.compile(loss='mean_squared_error', optimizer=type_optimizer,  metrics=[coeff_determination])
    return model

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

# Grid search with keras
num_batch_size = [1000,5000,10000,20000,30000,40000,5000]
num_epochs = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
neurons = [5,10,15,20,25,30,35,40,45,50]
model = KerasRegressor(build_fn=create_model)

# # model.fit(inputs,outputs,epochs=num_epoch,batch_size=n_batch_size)

param_grid = dict(batch_size=num_batch_size, 
                    neurons = neurons,
                    epochs=num_epochs)
grid = GridSearchCV(estimator=model, param_grid=param_grid,verbose=1,scoring='r2')
print(grid)
grid_result = grid.fit(inputs, outputs)

result = open('gridsearch.txt','w')
f1 = open('training.py')
print('Best: {} using {}'.format(grid_result.best_score_, grid_result.best_params_))
result.write('Best: {} using {}'.format(grid_result.best_score_, grid_result.best_params_))
means = grid_result.cv_results_['mean_test_score']
stds = grid_result.cv_results_['std_test_score']
params = grid_result.cv_results_['params']

result.write('\n')
for mean, std, param in zip(means, stds, params):
    print("%f (%f) with: %r\n" % (mean, std, param))

joblib.dump(grid_result.best_estimator_, 'gridsearch_best.pkl',compress=1)
joblib.dump(grid_result, 'gridsearch.pkl')