import matplotlib.pyplot as plt
from racetrack import racetrack
from lapsimple import lapsimple, car
import numpy as np

my_car = car(mu_x = 1, mu_y = 1, mass = 1000, power = 600, cLA = 6, cDA = 2, rho = 1.2)

rt = racetrack("tracks/Spielberg.csv")

sim = lapsimple(rt, my_car)

v_fwd = sim.intTrack(dir="fwd")
v_bwd = sim.intTrack(dir="bwd")
v = np.minimum(v_fwd, v_bwd)

plt.plot(rt.dist, v*3.6)
plt.gca().set_xlim(rt.dist.min(), rt.dist.max())
plt.title(f"Laptime: {sum(rt.u/v)}")
plt.show()