import sys
import os
import datetime
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from glob import glob


topo = xr.open_dataset('../topo/topo.nc')

x, y = np.meshgrid(topo.lon, topo.lat)

# Initialize particles where the depth is greater than 40m
depmask = (topo.bathymetry > 60.) * (x > 23.7)

lon_init = x[depmask].flatten()
lat_init = y[depmask].flatten()

nparticles = len(lon_init)

print("Number of particles: {}".format(nparticles))

with open('particles.in', 'w') as f:
    f.write(f"{nparticles}\n")
    for i in range(nparticles):
        f.write(f"{lon_init[i]} {lat_init[i]} 0.0 1.0 259200.0 1300.0 2.e-5\n")

plt.pcolormesh(topo.lon, topo.lat, topo.bathymetry[:-1, :-1], shading="flat")
plt.scatter(lon_init, lat_init, s=1, c='r')
plt.show()
