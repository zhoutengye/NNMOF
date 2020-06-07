import math

import time
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np
import sklearn
from sklearn import metrics
from sklearn.neural_network import MLPRegressor


# Load Training Date and test data
f = np.loadtxt('DATA/Linear_Training_2D.dat',unpack='true')
# f = np.loadtxt('DATA/MOF_Training.dat',unpack='true')
f2 = np.loadtxt('DATA/Linear_Test_2D.dat',unpack='true')

# X: training input data    y: training output data
# X2: test inut data 
# (It is weird that although I put 0:3 and 3:6, it is actually actually 0:2, 3:5 columns.
#    Not sure if there is somthing wrong on my laptop) 
X = np.transpose(f[0:3,:])
y = np.transpose(f[4:7,:])
# to change the range of test data, sinply change the subscripts or import other data.
n_test = 10
L = np.zeros((n_test,7))
X2 = np.transpose(f2[0:3,:])
Y3 = np.transpose(f2[4:7,:])
n_test = X2.shape[0]


# y = np.zeros((len(f[0]),2))
# X = np.transpose(f[3:6,:])
# y[:,0] =  np.transpose(f[0,:])
# y[:,1] =  np.transpose(f[2,:])
# # to change the range of test data, sinply change the subscripts or import other data.
# n_test = 10
# L = np.zeros((n_test,6))
# for i in range(0,n_test):
# 	L[i,:] = np.transpose(f2[:,np.random.randint(low = 1,high = len(f2[0]))])
# 	X2 = L[:,3:6]
# 	Y3 = L[:,0:3]

# sklean stuff
# regr = MLPRegressor(hidden_layer_sizes=(50,), max_iter=10, alpha=1e-4,\
#                      solver='sgd', verbose=10, tol=1e-4, random_state=1,\
#                      learning_rate_init=.1)

# regr = MLPRegressor(activation='logistic',hidden_layer_sizes=(5,5), max_iter=1000, alpha=1e-3,\
#                      solver='sgd', verbose=10, tol=1e-4,  random_state=1,\
#                      learning_rate_init=.1,learning_rate='adaptive')

# f = open("DATA/OPT_NN.namelist", "w")
regr = MLPRegressor(hidden_layer_sizes=(131,),activation='relu',\
			solver='adam', alpha=0.0001, batch_size='auto', \
			learning_rate='constant', learning_rate_init=0.019, \
			power_t=0.5, max_iter=1000, shuffle=True, random_state=None, \
			tol=0.0001, verbose=False, warm_start=False, momentum=0.9, \
			nesterovs_momentum=True, early_stopping=False, validation_fraction=0.1,\
			 beta_1=0.9, beta_2=0.999, epsilon=1e-08)

# timer1 = time.clock()
# # Linear fit
regr.fit(X,y)

# timer2 = time.clock()

# # predict the results with test data input
y2 = regr.predict(X2)

# timer3 = time.clock()

L1_error_xnorm_dt  = sum(abs(y2[:,0]-Y3[:,0])) / n_test
L1_error_ynorm_dt  = sum(abs(y2[:,1]-Y3[:,1])) / n_test
L1_error_area_dt  = sum(abs(y2[:,2]-Y3[:,2]))  / n_test

print(L1_error_xnorm_dt)
print(L1_error_ynorm_dt)
print(L1_error_area_dt)



 
# print(' Input X Centroid, Y centroid, Area')
# print(X2)
# print(' ')

# # print norm X norm, Y norm, Intersecpt predicted by Linear Regression
# print(' predicted norm X norm, Y norm, Intersecpt')
# print(y2) 
# print(' ')

# # print exact norm X norm, Y norm, Intersecpt
# print(' exact norm X norm, Y norm, Intersecpt')
# print(Y3)
# print(' ')

# print(' CPU_time for fit')
# print(timer2-timer1)

# print(' CPU_time for predict')
# print(timer3-timer2)

# # for i in range(regr.n_layers_-1):
# np.savetxt('DATA/out.dat',regr.coefs_[0])

# print(regr.coefs_)
