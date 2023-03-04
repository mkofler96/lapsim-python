import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from racetrack import racetrack
from lapsimple import lapsimple, car
import pandas as pd
import numpy as np
import streamlit as st
import os

# todo: put this function into lapsim class together with calculate laptime
def seconds_to_time(seconds):
    minutes, seconds = divmod(seconds, 60)
    seconds, hundreds = divmod(seconds, 1)
    hundreds = int(hundreds * 100)
    return "{:02d}:{:02d}:{:03d}".format(int(minutes), int(seconds), hundreds)

parent_dir = os.path.abspath('.')  # get absolute path of parent directory
available_tracks = os.listdir(os.path.join(parent_dir, "tracks"))
st.set_page_config(page_title="Lap Time Simulation", layout="wide")


df = pd.DataFrame(
    [
       {"Name": "F1 Car", "Power": 1000, "Mass": 800, "cLA": 4.0, "cDA": 2.0, "mu_x": 3.0, "mu_y": 3.0},
       {"Name": "F2 Car", "Power": 620, "Mass": 750, "cLA": 4.0, "cDA": 2.0, "mu_x": 3.0, "mu_y": 3.0}
   ]
)

col1, col2 = st.columns(2)
available_tracks.sort()
track_list = [os.path.splitext(file)[0] for file in available_tracks]

with col1:
    track_choice = st.selectbox('Select Track', track_list)

    #mu_x = st.number_input('mu x', value=1)
    #mu_y = st.number_input('mu y', value=1)
    #mass = st.number_input('mass', value=1000)
    #power = st.number_input('power', value=300)
    #cLA = st.number_input('cLA', value=2)
    #cDA = st.number_input('cDA', value=1)
    edited_df = st.experimental_data_editor(df, num_rows="dynamic")

    fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
maxv = 0
for index, row in edited_df.iterrows():
    my_car = car(mu_x = row["mu_x"], mu_y =row["mu_y"], mass = row["Mass"], power = row["Power"], cLA = row["cLA"], cDA = row["cDA"], rho = 1.2)

    rt = racetrack(os.path.join(parent_dir, "tracks", track_choice)+".csv")

    sim = lapsimple(rt, my_car)
    v_fwd = sim.intTrack(dir="fwd")
    v_bwd = sim.intTrack(dir="bwd")
    v = np.minimum(v_fwd, v_bwd)
    if max(v)>maxv:
        maxv = max(v)
    laptime = round(sum(rt.u/v),3)

    ax1.plot(rt.dist, v*3.6, label=row["Name"] + f": {seconds_to_time(laptime)}")
    
    ax1.set_xlim(rt.dist.min(), rt.dist.max())
    ax1.legend()
#ax.set_title(f"Laptime: {sum(rt.u/v)}")

ax1.set_ylim(0,maxv*1.2*3.6)
ax1.set_ylabel("Velocity (km/h)")
ax1.set_xlabel("Distance (m)")
ax2.plot(rt.rawx, rt.rawy)
ax2.set_aspect('equal')
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

with col2:
    st.pyplot(fig1)
    st.pyplot(fig2)

