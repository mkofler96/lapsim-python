import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy import optimize
from TMSimple import TMSimple
from vehicle_model import vehicle_model, state
import streamlit as st
import pandas as pd
vehicle = vehicle_model("../parameters/test_model.ini")

vx = st.slider('Select Velocity in km/h', 50, 90, 65)
vx = vx/3.6
max_d = 10
max_b = 3
SR = 0

style_choice = st.selectbox(
    'How would you like to view the data?',
    ('classic', 'seaborn-whitegrid', 'dark_background', 'ggplot'))

plt.style.use(style_choice)
all_states = []
fig, ax = plt.subplots()

for beta in np.linspace(-3,3,11):
    states = []
    for delta_ in np.linspace(-max_d,max_d, 31):
        current_state = state(vehicle, vx = vx, beta=beta, SR=SR, delta=delta_)
        current_state.solve()
        current_state = current_state.vehicle.get_residual_forces(current_state.result[0], current_state.result[1], current_state, get_state=True)
        states.append(current_state)
        all_states.append(current_state)
    ax.plot(*zip(*([(state_n.ay, state_n.yaw_moment) for state_n in states])), label=f"beta: {beta}", color="blue")

for delta_ in np.linspace(-max_d,max_d,11):
    states = []
    for beta in np.linspace(-3,3,31):
        current_state = state(vehicle, vx = vx, beta=beta, SR=SR, delta=delta_)
        current_state.solve()
        current_state = current_state.vehicle.get_residual_forces(current_state.result[0], current_state.result[1], current_state, get_state=True)
        states.append(current_state)
        all_states.append(current_state)
    ax.plot(*zip(*([(state_n.ay, state_n.yaw_moment) for state_n in states])), label=f"delta: {delta_}", color="red")
#ax.legend()
#ax.axhline(y=0)
#ax.axhline(x=0)
plt.show()
ax.set_xlabel("Lateral Acceleration [m/s²]")
ax.set_ylabel("Yaw Moment [Nm]")

st.pyplot(fig)

states_df = pd.DataFrame([{"ay": s.ay, "yaw moment": s.yaw_moment, "delta": s.delta, "beta": s.beta, "info": s.warnings} for s in all_states])

st.dataframe(data=states_df)