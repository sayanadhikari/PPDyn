import ini
import time
import os.path
from os.path import join as pjoin
import numpy as np
import pandas as pd

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import plotly.graph_objects as go
import webbrowser
# from multiprocessing import Process

def dash_energy(input):
    path ='data'
    params  = ini.parse(open(input).read())
    N       = int(params['particles']['N'])    # Number of particles
    tmax    = float(params['time']['tmax'])
    realTime    = bool(params['diagnostics']['realTime'])

    pd.options.plotting.backend = "plotly"
    countdown = 20

    if os.path.exists(pjoin(path,'energy.txt')):
        time,energy = np.loadtxt(pjoin(path,'energy.txt'),unpack=True)
        data = np.stack((time, energy), axis=1)
        df = pd.DataFrame(data, columns=['timestep', 'Energy'])
        fig = df.plot(template = 'plotly_dark')
    else:
        fig = go.Figure(data=[go.Scatter(x=[], y=[])])
        fig.layout.template = 'plotly_dark'

    app = dash.Dash(__name__, update_title=None)  # remove "Updating..." from title
    app.layout = html.Div([
        html.H1("PPDyn: Energy"),
                dcc.Interval(
                id='interval-component',
                interval=1*1000, # in milliseconds
                n_intervals=0
            ),
        dcc.Graph(id='graph'),
    ])

    # Define callback to update graph
    @app.callback(
        Output('graph', 'figure'),
        [Input('interval-component', "n_intervals")]
    )
    def streamFig(value):
        global df
        if os.path.exists(pjoin(path,'energy.txt')):
            time,energy = np.loadtxt(pjoin(path,'energy.txt'),unpack=True)
            data = np.stack((time, energy), axis=1)
            df1 = pd.DataFrame(data, columns=['timestep', 'Energy'])
            fig = df1.plot(x= 'timestep', y='Energy',template = 'plotly_dark')
        else:
            fig = go.Figure(data=[go.Scatter(x=[], y=[])])
            fig.layout.template = 'plotly_dark'
        #fig.show()
        return(fig)

    # def run():
    #     app.scripts.config.serve_locally = True
    #     app.run_server(port = 8069, dev_tools_ui=True, debug=False,
    #                       dev_tools_hot_reload =True, threaded=False)
    # s = Process(target=run)
    # s.start()

    webbrowser.open('http://127.0.0.1:8069/')
    app.scripts.config.serve_locally = True
    app.run_server(port = 8069, dev_tools_ui=True, debug=False,
                      dev_tools_hot_reload =True, threaded=False)
