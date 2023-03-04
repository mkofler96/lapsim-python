import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from racetrack import racetrack
from lapsimple import lapsimple, car
import pandas as pd
import numpy as np
import streamlit as st
import os


parent_dir = os.path.abspath('.')  # get absolute path of parent directory
available_tracks = os.listdir(os.path.join(parent_dir, "tracks"))
st.set_page_config(page_title="Lap Time Simulation", layout="wide")


df = pd.DataFrame(
    [
       {"Name": "Formula 1 Car", "Power": 600, "Mass": 1000, "cLA": 4, "cDA": 1, "mu_x": 1, "mu_y": 1},
   ]
)

col1, col2 = st.columns(2)

with col1:
    track_choice = st.selectbox('Select Track', [os.path.splitext(file)[0] for file in available_tracks])

    #mu_x = st.number_input('mu x', value=1)
    #mu_y = st.number_input('mu y', value=1)
    #mass = st.number_input('mass', value=1000)
    #power = st.number_input('power', value=300)
    #cLA = st.number_input('cLA', value=2)
    #cDA = st.number_input('cDA', value=1)
    edited_df = st.experimental_data_editor(df, num_rows="dynamic")

    fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for index, row in edited_df.iterrows():
    my_car = car(mu_x = row["mu_x"], mu_y =row["mu_y"], mass = row["Mass"], power = row["Power"], cLA = row["cLA"], cDA = row["cDA"], rho = 1.2)

    rt = racetrack(os.path.join(parent_dir, "tracks", track_choice)+".csv")

    sim = lapsimple(rt, my_car)
    v_fwd = sim.intTrack(dir="fwd")
    v_bwd = sim.intTrack(dir="bwd")
    v = np.minimum(v_fwd, v_bwd)

    ax1.plot(rt.dist, v*3.6, label=row["Name"])
    ax1.set_ylim(bottom=0)
    ax1.set_xlim(rt.dist.min(), rt.dist.max())
    ax1.legend()
#ax.set_title(f"Laptime: {sum(rt.u/v)}")

ax2.plot(rt.rawx, rt.rawy)
ax2.set_aspect('equal')
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

with col2:
    st.pyplot(fig1)
    st.pyplot(fig2)

