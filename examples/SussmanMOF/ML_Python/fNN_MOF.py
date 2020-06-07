import math

import time
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np
import sklearn
from sklearn import metrics
from sklearn.neural_network import MLPRegressor
from sklearn import preprocessing
from sklearn.utils.extmath import safe_sparse_dot


# Load Training Date and test data
f = np.loadtxt('DATA/MOF_Training4.dat',unpack='true')

X = np.transpose(f[0:3,:])
y = np.transpose(f[4:5,:])
# to change the range of test data, sinply change the subscripts or import other data.
n_test = 1
L = np.zeros((n_test,7))
X2 = np.transpose(f[0:3,1:2])
Y3 = np.transpose(f[4:5,1:2])
n_test = X2.shape[0]

regr = MLPRegressor(hidden_layer_sizes=(1), activation='logistic', max_iter=30, alpha=0.0,\
                     verbose=False, tol=1e-4, random_state=0,\
                     learning_rate_init=.1)

regr.fit(X,y)
y2 = regr.predict(X2)

print(y2)
print(Y3)

# # for i in range(regr.n_layers_-1):
# np.savetxt('DATA/out.dat',regr.coefs_[0])

print(regr.coefs_)
print(regr.intercepts_)

ww = regr.coefs_

activations = [X2]
activations.extend(np.empty((1,3)))
activations.extend(np.empty((1,1)))

activations[1] = safe_sparse_dot(activations[0],ww[0]) + regr.intercepts_[0]
activations[1] = 1.0 / (1.0+np.exp(-activations[1]))
activations[2] = safe_sparse_dot(activations[1],ww[1]) + regr.intercepts_[1]
print(activations[-1])

fl = np.matmul(X2,ww[0]) + regr.intercepts_[0][:]

print(X2.shape)

n1 = fl
n1 = 1.0/(1.0+np.exp(-fl))
# n1 = np.tanh(fl)

print(n1)
print(fl)

fl = np.matmul(n1,ww[1]) + regr.intercepts_[1][:]

n2 = fl

print(fl)
print(-np.log(1.0/fl-1))
print((fl-regr.intercepts_[1][:])/ww[1])

print(regr.n_layers_)
