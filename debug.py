import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy import optimize
from TMSimple import TMSimple
from vehicle_model import vehicle_model, state
from joblib import Parallel, delayed
import pandas as pd
vehicle = vehicle_model("parameters/test_model.ini")

vx = 65/3.6
max_d = 10
max_b = 3
SR = 0



fig, ax = plt.subplots()
states = []
for beta in np.linspace(-3,3,11):
    states = []
    for delta_ in np.linspace(-max_d,max_d, 31):
        current_state = state(vehicle, vx = vx, beta=beta, SR=SR, delta=delta_)
        current_state.solve()
        current_state = current_state.vehicle.get_residual_forces(current_state.result[0], current_state.result[1], current_state, get_state=True)
        states.append(current_state)

    ax.plot(*zip(*([(state_n.ay, state_n.yaw_moment) for state_n in states])), label=f"beta: {beta}", color="blue")

for delta_ in np.linspace(-max_d,max_d,11):
    states = []
    for beta in np.linspace(-3,3,31):
        current_state = state(vehicle, vx = vx, beta=beta, SR=SR, delta=delta_)
        current_state.solve()
        current_state = current_state.vehicle.get_residual_forces(current_state.result[0], current_state.result[1], current_state, get_state=True)
        states.append(current_state)

    ax.plot(*zip(*([(state_n.ay, state_n.yaw_moment) for state_n in states])), label=f"delta: {delta_}", color="red")
#ax.legend()
ax.axhline(y=0)
ax.axhline(x=0)
plt.show()
