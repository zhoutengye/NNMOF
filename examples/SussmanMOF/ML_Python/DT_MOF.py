import math

import time
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np
import sklearn
from sklearn.tree import DecisionTreeRegressor



# Load Training Date and test data
f = np.loadtxt('DATA/MOF_Training1.dat',unpack='true')
f[2,:] = -f[2,:]

# X: training input data    y: training output data
# X2: test inut data 
# (It is weird that although I put 0:3 and 3:6, it is actually actually 0:2, 3:5 columns.
#    Not sure if there is somthing wron, on my laptop) 
X = np.transpose(f[3:6,:])
y = np.transpose(f[0:3,:])

# sklean stuff
regr =DecisionTreeRegressor(criterion='mse', splitter='best', max_depth=10, \
				min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0, \
				max_features='auto', random_state=None, max_leaf_nodes=None, \
				min_impurity_decrease=0.0, min_impurity_split=None, presort=False)

# Linear fit
regr.fit(X,y)


f = open("DATA/tree_parameters.namelist", "w")

print(regr.tree_.n_outputs)
f.write("&tree_Parameters\n")
f.write('    node_count = %d\n' %regr.tree_.node_count)
f.write('    n_features = %d\n' %regr.tree_.n_features)
f.write('    n_outputs = %d\n' %regr.tree_.n_outputs)
f.write('    max_depth = %d\n' %regr.tree_.max_depth)
f.write("/\n")
f.write('\n')
f.close()

f = open('DATA/tree_left.dat', "w")
for i in range(0,regr.tree_.node_count):
	f.write('%d\n' %regr.tree_.children_left[i])
f.close()
f = open('DATA/tree_right.dat', "w")
for i in range(0,regr.tree_.node_count):
	f.write('%d\n' %regr.tree_.children_right[i])
f.close()
f = open('DATA/tree_feature.dat', "w")
for i in range(0,regr.tree_.node_count):
	f.write('%d\n' %regr.tree_.feature[i])
f.close()
f = open('DATA/tree_threshold.dat', "w")
for i in range(0,regr.tree_.node_count):
	f.write('%12.8f\n' %regr.tree_.threshold[i])
f.close()
for j in range(0,regr.tree_.n_outputs):
	jj = j+1
	f = open('DATA/tree_value'+str(jj)+'.dat', "w")
	for i in range(0,regr.tree_.node_count):
		f.write('%12.8f\n' %regr.tree_.value[i,j])

# for i in range(0,regr.tree_.node_count):
# 	ii = i+1
# 	f = open('DATA/tree/node_'+str(ii)+'.namelist', "w")
# 	f.write('&node0\n')
# 	f.write('    children_left0 = %d\n' %regr.tree_.children_left[i])
# 	f.write('    children_right0 = %d\n' %regr.tree_.children_right[i])
# 	f.write('    features0 = %d\n' %regr.tree_.feature[i])
# 	f.write('    threshold0 = %.12f\n' %regr.tree_.threshold[i])
# 	f.write('    value0 = ' %regr.tree_.value[i])
# 	for j in range(0,regr.tree_.n_outputs):
# 		f.write('%.12f, ' %regr.tree_.value[i,j])
# 	f.write("\n")
# 	f.write("/\n")
# 	f.write('\n')
# 	f.close()

# Tnode = 0

# s = [[0.68641093,  0.03268939, -0.15239356]]

# for i in range(0,regr.tree_.max_depth):

	# print( Tnode, regr.tree_.feature[Tnode], regr.tree_.threshold[Tnode], regr.tree_.children_left[Tnode], regr.tree_.children_right[Tnode])

	# if s[0][regr.tree_.feature[Tnode]] < regr.tree_.threshold[Tnode]:
		# Tnode = regr.tree_.children_left[Tnode]
	# else:
		# Tnode = regr.tree_.children_right[Tnode]

	
# print(regr.tree_.value[Tnode])

# print(regr.predict(s))
