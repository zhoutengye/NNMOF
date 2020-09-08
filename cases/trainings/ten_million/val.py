import numpy as np
from keras.models import load_model
from backwardc import floodsz_backwardc, angle2norm

def coeff_determination(y_true, y_pred):
    from keras import backend as K
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

model = load_model('test.h5', custom_objects={'coeff_determination': coeff_determination})
# model = load_model('very_deep.h5', custom_objects={'coeff_determination': coeff_determination})

data_dir = '.'

exact_f = np.load(data_dir+'/exact_f.npy')
exact_angle = np.load(data_dir+'/exact_angle.npy')
exact_centroid = np.load(data_dir+'/exact_centroid.npy')
delta_angle = np.load(data_dir+'/delta_angle.npy')
initial_angle = np.load(data_dir+'/initial_angle.npy')

exact_centroid = exact_centroid - 0.5

exact_f = exact_f.reshape([len(exact_f),1])
inputs = np.hstack((initial_angle,exact_f))
outputs = delta_angle.copy()

tn = 1

pred = model.predict(inputs[tn-1:tn,:]) + initial_angle[tn-1:tn,:]

exact_norm = angle2norm(exact_angle[tn-1:tn,:])
initial_norm = angle2norm(initial_angle[tn-1:tn,:])
nn_norm = angle2norm(pred)

print(exact_f[tn-1][0])
ce = floodsz_backwardc(exact_norm,exact_f[tn-1]) -0.5
ce1 = floodsz_backwardc(initial_norm,exact_f[tn-1]) - 0.5
ce2 = floodsz_backwardc(nn_norm,exact_f[tn-1])-  0.5
dce1 = ce1 - ce
dce2 = ce2 - ce
print(np.linalg.norm(dce1,2))
print(np.linalg.norm(dce2,2))

print(exact_centroid[tn-1,:])
print(ce)
print(ce1)
print(ce2)
