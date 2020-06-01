import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, GaussianNoise, Dropout, AlphaDropout
import numpy as np
from tensorflow.keras.models import Model
import tensorflow.keras.backend as K
from tensorflow.keras.callbacks import EarlyStopping
from tempfile import TemporaryFile
from tensorflow.keras.models import load_model


#tailles des representations encodees
encoding_dim1 = 32
encoding_dim2 = 64
encoding_dim3 = 64



#fonction perceptron

def neuralnet(Data, pop, neuralnet=None): #Data -> input, pop -> annotations
  
  es = EarlyStopping(monitor='val_loss', mode='min', restore_best_weights=True, patience=20) #termine training avant overfitting
  
  if neuralnet is None:
    input=Input(shape=(Data.shape[1],)) #data d entree sous forme de tenseur
   
    #input_noise=GaussianNoise(.2)(input) #ajout de bruit
    
    
    #reduction des dimensions
    #dropout lors entrainement
    '''
    encoded1 = Dense(encoding_dim1, activation='relu')(input)
    drop_out1 = Dropout(0.1)(encoded1)
    encoded2 = Dense(encoding_dim2, activation='relu')(drop_out1) #prend la premier encoded en entree
    drop_out2 = Dropout(0.5)(encoded2)
    encoded3 = Dense(encoding_dim3, activation='relu')(drop_out2) #prend la premier encoded en entree
    
    drop_out3 = Dropout(0.5)(encoded3)
    '''
    encoded1 = Dense(encoding_dim1, activation='selu')(input)
    drop_out1 = AlphaDropout(0.1)(encoded1)
    encoded2 = Dense(encoding_dim2, activation='selu')(drop_out1) #prend la premier encoded en entree
    drop_out2 = AlphaDropout(0.3)(encoded2)
    encoded3 = Dense(encoding_dim3, activation='selu')(drop_out2) #prend la premier encoded en entree
    
    drop_out3 = AlphaDropout(0.2)(encoded3)
   
    auxiliary_output1 = Dense(pop.shape[1], activation='softmax', name='aux_output1')(drop_out3) #pour mieux visualiser les data : data doivent etre separees lineairement
    #auxiliary_output2 = Dense(pop.shape[1], activation='softmax', name='aux_output2')(encoded2)
    #auxiliary_output3 = Dense(pop.shape[1], activation='softmax', name='aux_output3')(encoded3)
    
    
    neuralnet = Model(input, outputs=[auxiliary_output1]) # auxiliary_output2, auxiliary_output3]) #formation du modele -> en entree : input et en sortie : decoded3 et auxiliary_output
    
     
    neuralnet.compile(optimizer='adam', 
    loss=['categorical_crossentropy'],  
    metrics={'aux_output1':'accuracy'}) 
    
    neuralnet.summary()
    
  
  entrainement = neuralnet.fit(Data, pop, epochs=200, batch_size=256, callbacks=[es], shuffle=True, validation_split=0.1,verbose = 1) #fit de autoencoder
  

  
  return neuralnet.predict(Data), neuralnet
  
  

