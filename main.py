import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy import optimize
from TMSimple import TMSimple
from vehicle_model import vehicle_model, state
from joblib import Parallel, delayed
import pandas as pd
vehicle = vehicle_model("parameters/test_model.ini")

states = []
results = []
processes = []
n = 0
for vx in np.linspace(25/3.6, 100/3.6, 10):
    for beta in np.linspace(-3, 3, 10):
        for SR in np.linspace(-100,100, 10):
            for delta in np.linspace(-30, 30, 10):
                current_state = state(vehicle, vx = vx, beta=beta, SR=SR, delta=delta)
                states.append(current_state)
                

# 113 sekunden für de variante
#for state1 in states:
#    res = optimize.fsolve(lambda x: vehicle.get_residual_forces(x[0], x[1], state1), [0,0])
#    res = np.append(res, state1.vx)
#    results.append(res)
# 107 sekunden ohne multiprocessing
 
#states[0].solve()

#print(states[0])

#for state1 in states:
#    state1.solve()

with Parallel(n_jobs=6, verbose=5) as parallel:
    delayed_funcs = [delayed(lambda x:x.solve())(run) for run in states]
    results = parallel(delayed_funcs)

states_df = pd.DataFrame([vars(state_n) for state_n in results])
states_df.to_csv("statees.csv", encoding='utf-8', index=False)
#fig, ax = plt.subplots()
#results = np.array(results)
#ax.scatter(results[:,0], results[:,1], results[:,2])
#plt.show()