import numpy as np
from keras.models import load_model

def coeff_determination(y_true, y_pred):
    from keras import backend as K
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

# model = load_model('test.h5', custom_objects={'coeff_determination': coeff_determination})
model = load_model('very_deep.h5', custom_objects={'coeff_determination': coeff_determination})

data_dir = '.'

exact_f = np.load(data_dir+'/exact_f.npy')
exact_angle = np.load(data_dir+'/exact_angle.npy')
exact_centroid = np.load(data_dir+'/exact_centroid.npy')
delta_angle = np.load(data_dir+'/delta_angle.npy')
initial_angle = np.load(data_dir+'/initial_angle.npy')

exact_centroid = exact_centroid - 0.5

exact_f = exact_f.reshape([len(exact_f),1])
inputs = np.hstack((exact_centroid,exact_f))
outputs = delta_angle.copy()

print(outputs[:10,:]-model.predict(inputs[:10,:]))
print(delta_angle[:10,:])
