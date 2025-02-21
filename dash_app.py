# !pip install jupyter-dash
# based on trap_app_v4.2

# !pip install jupyter-dash

import plotly.express as px
from jupyter_dash import JupyterDash
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import base64
import datetime
import os
import io
import plotly.graph_objs as go
import dash
import dash_table
import plotly.express as px
import torch
import numpy as np
import pandas as pd
import tensorflow as tf
import csv
import os
import re
import keras
from tensorflow.keras import layers
from keras.models import Model,load_model
from tensorflow.keras.optimizers import Adam
from keras.models import Sequential, Model
from keras_preprocessing.sequence import pad_sequences
import transformers
from transformers import T5Tokenizer
from transformers import TFT5EncoderModel
import sentencepiece
from sklearn.tree import DecisionTreeClassifier


#cache = Cache(app.server, config={
#    'CACHE_TYPE' : 'filesystem',
#    'CACHE_DIR' : '/trap/cache'
#})
#os.chdir(r"/trap")


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = JupyterDash(__name__, external_stylesheets=external_stylesheets)

def blank_fig():
    fig = go.Figure(go.Scatter(x=[], y = []))
    fig.update_layout(template = None)
    fig.update_xaxes(showgrid = False, showticklabels = False, zeroline=False)
    fig.update_yaxes(showgrid = False, showticklabels = False, zeroline=False)
    
    return fig

colors = {
    "graphBackground": "#F5F5F5",
    "background": "#ffffff",
    "text": "#000000"
}

app.layout = html.Div([
    html.H2("TRAP: Deep learning platform for CD8+ T-cell epitope prediction",
           style = {'color': 'black', 'fontSize': 30, 'font-family' : 'sans-serif',
                   'textAlign': 'center'}),

    dcc.Dropdown(
        options=[
            {"label": "Pathogenic model", "value": "pathogenic"},
            {"label": "Self-antigen model", "value": "self"},
        ],
        value='pathogenic',
        placeholder="Select a model",
        id='dropdown',
        style = {
            'textAlign': 'center'}
    ),
    
    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    ),
    dcc.Graph(id='bargraph',
              figure = blank_fig(),
              style = {
                 "margin-left": "auto",
                 "margin-right": "auto",
             }),
    html.Div(id='output-data-upload')
    ],
    style={"max-width": "900px"},
)


# load models 
os.environ["CURL_CA_BUNDLE"]=""
os.chdir('/trap')
path_trap = tf.keras.models.load_model('model/pathogenic_trap_model')
path_trap_softmax = tf.keras.models.load_model('model/pathogenic_trap_softmax_model')

self_trap = tf.keras.models.load_model('model/self_antigen_trap_model')
self_trap_softmax = tf.keras.models.load_model('model/self_antigen_trap_softmax_model')


files = [f for f in os.listdir('model') if re.match(r'pathogenic_trap_softmaxbucket_model', f)]
files = ['model/' + s  for s in files]
path_trap_softmaxbucket = []
for f in files:
    print(f)
    model = tf.keras.models.load_model(f)
    path_trap_softmaxbucket.append(model)
    
files = [f for f in os.listdir('model') if re.match(r'self_antigen_trap_softmaxbucket_model', f)]
files = ['model/' + s  for s in files]
self_trap_softmaxbucket = []
for f in files:
    print(f)
    model = tf.keras.models.load_model(f)
    self_trap_softmaxbucket.append(model)
    
    
##### OOD classifiers #########
# pathogenic 
path_ood_data = pd.read_csv('data/cal_ood_data_pathogenic.csv')
p_X_train = path_ood_data[['max', 'ensemble_mean_maxprob']]
p_y_train = path_ood_data.iloc[:,-1]
path_trap_ood = DecisionTreeClassifier().fit(p_X_train, p_y_train)

# self-antigen
self_ood_data = pd.read_csv('data/cal_ood_data_selfantigen.csv')
s_X_train = self_ood_data[['max', 'ensemble_mean_maxprob']]
s_y_train = self_ood_data.iloc[:,-1]
self_trap_ood = DecisionTreeClassifier().fit(s_X_train, s_y_train)

######### tokenizer and embedding ###########
tokenizer = T5Tokenizer.from_pretrained("rostlab/prot_t5_xl_uniref50", do_lower_case=False )
TFT5EncoderModel = TFT5EncoderModel.from_pretrained("rostlab/prot_t5_xl_uniref50", from_pt=True)



def parse_data(contents, filename):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV or TXT file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
        elif 'txt' or 'tsv' in filename:
            # Assume that the user upl, delimiter = r'\s+'oaded an excel file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')), delimiter = r'\s+')
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return df


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

def preprocess_test_peptides(test_data):

    # add contact positions 
    cont_peptides = []
    hydrophobicity = []
    for pep in test_data.Peptide:
      cont_pep = pep[2:-1]
      cont_peptides.append(cont_pep)
      hyd_counts = pep.count('A') + pep.count('V') + pep.count('L') + pep.count('M') + pep.count('W')
      length = len(pep)
      hydrophobicity.append(hyd_counts/length)

    cont_peptides = pd.DataFrame(cont_peptides); cont_peptides.columns = ['ContactPosition']
    hydrophobicity = pd.DataFrame(hydrophobicity); hydrophobicity.columns = ['Hydrophobicity']
    test_data = pd.concat([test_data, cont_peptides, hydrophobicity], axis=1)

    return test_data

def embed_test_peptides(test_data, tokenizer, TFT5EncoderModel):
    # embed test peptides 
    test_peptides = test_data.ContactPosition.values
    input_peptides = add_space_to_pep(test_peptides)
      
    tokenized_texts = [tokenizer.tokenize(sent) for sent in input_peptides]
    input_ids = pad_sequences([tokenizer.convert_tokens_to_ids(txt) for txt in tokenized_texts], padding="pre")

    attention_masks = []
    for seq in input_ids:
      seq_mask = [float(i>0) for i in seq]
    attention_masks.append(seq_mask)

    embedding = model(input_ids=input_ids)[0]
    embedding = np.asarray(embedding)

    return embedding

def predict_trap(emebedding, test_data, trap, trap_softmax, trap_softmaxbucket, trap_ood):
    
  X_test = emebedding
  X_test_mlp = test_data[['Hydrophobicity', 'nlog2Rank']]

  # trap 
  y_pred = trap.predict([X_test, X_test_mlp])

  # on softmax 
  y_pred_softmax = trap_softmax.predict([X_test, X_test_mlp])
  softmax = np.max(y_pred_softmax,axis = 1)

  # ensemble on softmax
  bucket_softmax = []
  for i in range(10):
    random_model = trap_softmaxbucket[i]
    random_model_score = Model(random_model.input, random_model.get_layer('logits').output)
    X_random_test_logits = random_model_score.predict(x = [X_test, X_test_mlp])
    smax = np.max(X_random_test_logits,axis = 1)
    bucket_softmax.append(smax)

  c_softmax = np.vstack((bucket_softmax)).T
  ensemble = np.mean(c_softmax, axis = 1)

  # process ood output
  softmax_dt = pd.DataFrame(softmax); softmax_dt.columns = ['max']
  ensemble_dt = pd.DataFrame(ensemble); ensemble_dt.columns = ['ensemble_mean_maxprob']
  ood_test_data = pd.concat([softmax_dt, ensemble_dt], axis=1)

  # ood classifier
  X_test = ood_test_data
  y_pred_ood = trap_ood.predict(X_test)

  # process outputs
  trap_pred = pd.DataFrame(y_pred); trap_pred.columns = ['TRAP']
  ood_test_data.columns = ['MaxProb', 'Ensemble']

  ood_dt = pd.DataFrame(y_pred_ood); ood_dt.columns = ['Confidence']
  ood_dt['Confidence'] = np.where(ood_dt['Confidence']==1, 'High', 'Low')
  prediction_dt = pd.concat([test_data[['Peptide', 'ContactPosition', 'nlog2Rank']], trap_pred, ood_test_data, ood_dt], axis=1)

  return prediction_dt


def blank_fig():
    fig = go.Figure(go.Scatter(x=[], y = []))
    fig.update_layout(template = None)
    fig.update_xaxes(showgrid = False, showticklabels = False, zeroline=False)
    fig.update_yaxes(showgrid = False, showticklabels = False, zeroline=False)
    
    return fig



@app.callback(Output('bargraph', 'figure'),
            [
                Input('upload-data', 'contents'),
                Input('upload-data', 'filename'),
                Input('dropdown', 'value'),
            ])
def update_figure(contents, filename, dropdown):
    
    if contents:
        print('Making predictions...')
        
        contents = contents[0]
        filename = filename[0]
        df = parse_data(contents, filename)
        
        test_data = df[['Peptide', 'nlog2Rank']]
        test_data_p = preprocess_test_peptides(test_data)
        embedding = embed_test_peptides(test_data_p, tokenizer, TFT5EncoderModel)
        if dropdown == 'self':
            df = predict_trap(embedding, test_data_p, self_trap, self_trap_softmax, 
                              self_trap_softmaxbucket, self_trap_ood)
        elif dropdown == 'pathogenic':
            df = predict_trap(embedding, test_data_p, path_trap, path_trap_softmax, 
                              path_trap_softmaxbucket, path_trap_ood)
        
        # plot 
        prediction = df.drop_duplicates()
        prediction_dt = prediction.nlargest(30,'TRAP')
        fig_output = px.bar(prediction_dt, 
                      x="Peptide", 
                      y="TRAP", 
                      color="Confidence", 
                      width=800, height=400, template="simple_white",
                      color_discrete_map={ # replaces default color mapping by value
                            "High": "darkcyan", "Low": "darkgoldenrod"
                      }).add_shape( # add a horizontal "target" line
                 type="line", line_color="salmon", line_width=3, opacity=1, line_dash="dot",
                 x0=0, x1=1, xref="paper", y0=0.5, y1=0.5, yref="y")
        fig_output.update_layout(transition_duration=500)
    
    else: fig_output = blank_fig()
        
    style = {"margin-left": "auto",
             "margin-right": "auto"}
    
    return fig_output



@app.callback(Output('output-data-upload', 'children'),
            [
                Input('upload-data', 'contents'),
                Input('upload-data', 'filename'),
                Input('dropdown', 'value'),
            ])
def update_table(contents, filename, dropdown):
    table = html.Div()

    if contents and dropdown is not None:
        contents = contents[0]
        filename = filename[0]
        df = parse_data(contents, filename)
        
        test_data = df[['Peptide', 'nlog2Rank']]
        test_data_p = preprocess_test_peptides(test_data)
        embedding = embed_test_peptides(test_data_p, tokenizer, TFT5EncoderModel)
        
        if dropdown == 'self':
            df = predict_trap(embedding, test_data_p, self_trap, self_trap_softmax, 
                              self_trap_softmaxbucket, self_trap_ood)
        elif dropdown == 'pathogenic':
            df = predict_trap(embedding, test_data_p, path_trap, path_trap_softmax, 
                              path_trap_softmaxbucket, path_trap_ood)
            
            
        df = df.drop_duplicates().sort_values(by=['TRAP'], ascending=False)
        df['nlog2Rank'] = df['nlog2Rank'].apply(lambda x: float("{:.3f}".format(x)))
        df['TRAP'] = df['TRAP'].apply(lambda x: float("{:.5f}".format(x)))
        df['MaxProb'] = df['MaxProb'].apply(lambda x: float("{:.3f}".format(x)))
        df['Ensemble'] = df['Ensemble'].apply(lambda x: float("{:.3f}".format(x)))

        table = html.Div([
            html.H5(filename, style = {'fontSize': 18, 'font-family' : 'sans-serif', 
                                      'color' : 'darkgrey'}),
            dash_table.DataTable(
                data=df.to_dict('rows'),
                columns=[{'name': i, 'id': i} for i in df.columns],
                style_cell={'textAlign': 'center', 'font-family' : 'sans-serif'},
                style_header={
                    #'backgroundColor': 'white',
                    'fontWeight': 'bold', 'font-family' : 'sans-serif'},
                export_format="csv"
            ),
            #html.Hr(),
            #html.Div('Raw Content'),
            #html.Pre(contents, style={
            #    'whiteSpace': 'pre-wrap',
            #    'wordBreak': 'break-all'
            #})
        ])

    return table




if __name__ == '__main__':
    app.run_server(debug=True) # or whatever you choose



