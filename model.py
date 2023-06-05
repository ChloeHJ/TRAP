import numpy as np
import pandas as pd
import tensorflow as tf
import sys
import csv
import random
import os as os
from io import StringIO
import keras
from tensorflow.keras import layers
from keras.layers import Input,Dense,concatenate,Dropout,AveragePooling2D
from keras.models import Model,load_model                                                      
from keras import backend as K
from tensorflow.keras.optimizers import Adam
import matplotlib.pyplot as plt
from sklearn.metrics import auc,precision_recall_curve,roc_curve,confusion_matrix
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Activation, Flatten, Embedding, Input, LSTM, RNN, SimpleRNN, Bidirectional, Concatenate, GlobalMaxPool1D, Conv1D, MaxPooling1D, GlobalMaxPooling1D, Conv2D
from keras.utils.vis_utils import plot_model
import transformers
from transformers import T5Tokenizer, T5Model, BertModel, BertForMaskedLM, BertTokenizer, pipeline, BertConfig, EncoderDecoderConfig, EncoderDecoderModel
from transformers import TFBertModel, TFXLNetModel, XLNetTokenizer,XLNetConfig, TFT5EncoderModel
import re
import pickle
import seaborn as sb
import matplotlib.pyplot as plt
from keras.preprocessing.sequence import pad_sequences

print('transformers version: %s' %transformers.__version__)
print('tensorflow version: %s' %tf.__version__)
print('keras version: %s' %keras.__version__)


# plotting scripts 
from tensorflow.keras.optimizers import Adam

def evaluate_roc_pr_mlp_10cv(test_model, X_train, X_train_mlp, y_train):
    # let's do a train/validation split
    bucket_roc = []
    bucket_pr = []
    for i in range(10):
        array = np.arange(len(X_train))
        train_index = np.random.choice(array,int(len(X_train)*0.9),replace=False)
        valid_index = [item for item in array if item not in train_index]

        input1_train = X_train[train_index]
        input1_valid = X_train[valid_index]
        input1_mlp_train = X_train_mlp.iloc[train_index]
        input1_mlp_valid = X_train_mlp.iloc[valid_index]
        label_train = y_train[train_index]
        label_valid = y_train[valid_index]

        model = test_model(X_train)
        callback_val = keras.callbacks.EarlyStopping(monitor='val_loss', patience=15,restore_best_weights=False)
        callback_train = keras.callbacks.EarlyStopping(monitor='loss',patience=2,restore_best_weights=False)
        history = model.fit(
            x=[input1_train, input1_mlp_train],   # feed a list into
            y=label_train,
            validation_data = ([input1_valid, input1_mlp_valid],label_valid),
            batch_size=128,
            epochs=200,
            callbacks = [callback_val,callback_train])

        y_true = label_valid
        y_pred = model.predict([input1_valid, input1_mlp_valid])

        fpr,tpr,_ = roc_curve(y_true,y_pred)
        area = auc(fpr,tpr)
        bucket_roc.append((fpr,tpr,_,area))

        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        area = auc(recall, precision)
        bucket_pr.append((precision, recall, _, area))

    draw_roc(bucket_roc)
    draw_pr(bucket_pr)


def add_space_to_pep(peptides):
    peptide_space = [] 
    for ele in peptides:
        temp = [[]]
        for char in ele:
            temp.append([])
            temp[-1].append(char) 
        peptide_space.append(' '.join(''.join(ele) for ele in temp))
    peptide_space = [re.sub(r"[UZOB]", "X", sequence.lstrip()) for sequence in peptide_space]

    return peptide_space


# load and process data 
sdata = pd.read_csv('pathogenic_db.csv')
sdata['Immunogenicity'] = sdata['Immunogenicity'].replace(['Positive', 'Positive-Low', 'Positive-Intermediate', 'Positive-High', 'Negative'], [1, 1,1,1,0])
y_train  = sdata['Immunogenicity'].values
peptides = sdata.ContactPosition.values
input_peptides = add_space_to_pep(peptides)

# load embedding model 
tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False )
model = TFT5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50", from_pt=True)

# embed peptides 
tokenized_texts = [tokenizer.tokenize(sent) for sent in input_peptides]
input_ids = pad_sequences([tokenizer.convert_tokens_to_ids(txt) for txt in tokenized_texts], padding="pre")

attention_masks = []
for seq in input_ids:
  seq_mask = [float(i>0) for i in seq]
  attention_masks.append(seq_mask)

embedding = model(input_ids=input_ids)[0]
embedding = np.asarray(embedding)
pickle.dump( embedding, open( "protT5_xl_peptides.pkl", "wb" ) )

# load embedded sequence 
protT5_xl_peptides_1 = pickle.load( open( "protT5_xl_peptides.pkl", "rb" ) )
print(protT5_xl_peptides_1.shape)

# cnn + hydrophobicity + mhc binding model 
X_train_mlp = sdata[['hydrophobicity', 'nlog2Rank']]
bert_embed_matrix = protT5_xl_peptides_1 
def cnn1d_1_MH(bert_embed_matrix):
    
    s = keras.Input(shape=(bert_embed_matrix.shape[1],  bert_embed_matrix.shape[2]))
    emb = bert_embed_matrix.shape[2]

    mlp_input = keras.Input(shape=(2,))

    pep_conv1 = Conv1D(emb, 1, padding='same', activation='relu', kernel_initializer='glorot_normal', name = 'kernel_1')(s)
    pep_pool1 = GlobalMaxPooling1D()(pep_conv1)
    pep_conv3 = Conv1D(emb, 3, padding='same', activation='relu', kernel_initializer='glorot_normal', name = 'kernel_3')(s)
    pep_pool3 = GlobalMaxPooling1D()(pep_conv3)
    pep_conv5 = Conv1D(emb, 5, padding='same', activation='relu', kernel_initializer='glorot_normal', name = 'kernel_5')(s)
    pep_pool5 = GlobalMaxPooling1D()(pep_conv5)
    pep_conv7 = Conv1D(emb, 7, padding='same', activation='relu', kernel_initializer='glorot_normal', name = 'kernel_7')(s)
    pep_pool7 = GlobalMaxPooling1D()(pep_conv7)
    pep_cat = concatenate([pep_pool1, pep_pool3, pep_pool5, pep_pool7])

    mlp = Dense(2000, activation='relu')(mlp_input)
    merge = concatenate([pep_cat, mlp])

    dense = Dense(256, activation='relu')(merge) 
    dense = Dropout(0.1)(dense)
    out = Dense(1, activation='sigmoid')(dense)
    model = Model(inputs=[s, mlp_input], outputs=[out])

    model.compile(loss='binary_crossentropy',
                  optimizer= Adam(learning_rate = 0.00001, decay = 1e-6 ),
                  metrics='accuracy')
    
    return model

model = cnn1d_1_MH(bert_embed_matrix)
model.summary()


X_train = protT5_xl_peptides_1
evaluate_roc_pr_mlp_10cv(cnn1d_1_MH, X_train, X_train_mlp, y_train)
