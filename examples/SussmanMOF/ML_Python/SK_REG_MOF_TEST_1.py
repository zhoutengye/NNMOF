import math

import time
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np
import sklearn
from sklearn import metrics
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neural_network import MLPRegressor


# Load Training Date and test data
f = np.loadtxt('DATA/Linear_Training_2D.dat',unpack='true')
f2 = np.loadtxt('DATA/Linear_Test_2D.dat',unpack='true')
# f = np.loadtxt('DATA/MOF_Training.dat',unpack='true')
# f2 = np.loadtxt('DATA/MOF_test.dat',unpack='true')

# X: training input data    y: training output data
# X2: test inut data 
# (It is weird that although I put 0:3 and 3:6, it is actually actually 0:2, 3:5 columns.
#    Not sure if there is somthing wrong on my laptop) 
kkk = 3

X = np.transpose(f[0:3,:])
y = np.transpose(f[kkk+3,:])

# to change the range of test data, sinply change the subscripts or import other data.
n_test = 10000
# L = np.zeros((n_test,7))
# for i in range(0,n_test):
# 	L[i,:] = np.transpose(f2[:,np.random.randint(low = 1,high = len(f2[0]))])
# 	X2 = L[:,0:3]
# 	Y3 = L[:,kkk+3]
X2 = np.transpose(f2[0:3,:])
Y3 = np.transpose(f2[kkk+3,:])
n_test = X2.shape[0]
print(n_test)

# X = np.transpose(f[3:6,:])
# y = np.transpose(f[0:3,:])
# # to change the range of test data, sinply change the subscripts or import other data.
# n_test = 10000
# L = np.zeros((n_test,6))
# for i in range(0,n_test):
# 	L[i,:] = np.transpose(f2[:,np.random.randint(low = 1,high = len(f2[0]))])
# 	X2 = L[:,3:6]
# 	Y3 = L[:,0:3]
	

# y = np.zeros((len(f[0]),2)
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

#### Linear Model ####
lm = linear_model.LinearRegression()

timer1_lm = time.clock()

lm.fit(X,y)

timer2_lm = time.clock()

y2_lm = lm.predict(X2)

timer3_lm = time.clock()


#### Decision Tree Model####
dt = DecisionTreeRegressor()

timer1_dt = time.clock()

dt.fit(X,y)

timer2_dt = time.clock()

y2_dt = dt.predict(X2)

timer3_dt = time.clock()


#### Nueron Network Model####
nn = MLPRegressor()

timer1_nn = time.clock()

nn.fit(X,y)

timer2_nn = time.clock()

y2_nn = nn.predict(X2)

timer3_nn = time.clock()


#### K Nearest Neighbour Model####
knn = KNeighborsRegressor(n_neighbors=5)

timer1_knn = time.clock()

knn.fit(X,y)

timer2_knn = time.clock()

y2_knn = knn.predict(X2)

timer3_knn = time.clock()


#### Random_Forest Model####
rf = RandomForestRegressor()

timer1_rf = time.clock()

rf.fit(X,y)

timer2_rf = time.clock()

y2_rf = rf.predict(X2)

timer3_rf = time.clock()

#### L_1 error estimation ####
L1_error_lm  = sum(abs(y2_lm[:]-Y3[:])) / n_test
L1_error_nn  = sum(abs(y2_nn[:]-Y3[:])) / n_test
L1_error_dt  = sum(abs(y2_dt[:]-Y3[:])) / n_test
L1_error_knn = sum(abs(y2_knn[:]-Y3[:]))/ n_test
L1_error_rf  = sum(abs(y2_rf[:]-Y3[:])) / n_test




print('Training CPU_time ')
print('LM_Trainging   ',timer2_lm-timer1_lm)
print('NN_Trainging   ',timer2_nn-timer1_nn)
print('DT_Trainging   ',timer2_dt-timer1_dt)
print('KNN_Trainging  ',timer2_knn-timer1_knn)
print('RF_Trainging   ',timer2_rf-timer1_rf)
print(' ')
print('Fitting CPU_time ')
print('LM_Fitting   ',timer3_lm-timer2_lm)
print('NN_Fitting   ',timer3_nn-timer2_nn)
print('DT_Fitting   ',timer3_dt-timer2_dt)
print('KNN_Fitting  ',timer3_knn-timer2_knn)
print('RF_Fitting   ',timer3_rf-timer2_rf)

print(L1_error_lm)
print(L1_error_nn)
print(L1_error_dt)
print(L1_error_knn)
print(L1_error_rf)
print('')

print(dt.n_outputs_)
print(dt.n_features_)
print(dt.tree_)