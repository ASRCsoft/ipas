# ipas
A Python implementation of the [Ice Particle Aggregate Simulator](http://www.carlgschmitt.com/Microphysics.html)

## Installation

```
pip install git+https://github.com/ASRCsoft/ipas_python.git
```

## Manipulating crystals (monomers) and clusters (aggregates)

```python
import numpy as np
import ipas.crystals as crys

# create a hexagonal crystal centered at (1,0,0)
crystal = crys.IceCrystal(length=4, width=6, center=[1, 0, 0])
crystal.points # get a numpy array containing the crystal vertices
# rotate the crystal 45 degrees around the y-axis and 90 degrees
# around the z-axis
crystal.rotate_to([0,np.pi/4,np.pi/2])
# return a shapely MultiLineString representing the crystal edges,
# which plots automatically in a jupyter notebook
crystal.plot()
# project the crystal onto the xy plane, returning a shapely Polygon
poly = crystal.projectxy()
poly.area # area of the polygon

# use the crystal to start a cluster
cluster = crys.IceCluster(crystal)
# add a new crystal to the cluster
crystal2 = crys.IceCrystal(length=4, width=6)
cluster.add_crystal_from_above(crystal2)
cluster.plot() # same as IceCrystal.plot()
# plot the 2D projection of the cluster along with an ellipse fit to
# the area of the projection
cluster.plot_ellipse()
```

## Simulating clusters

```python
import ipas.lab as lab
# create 100 clusters of 2 crystals each with a simulation similar to
# the IDL IPAS code
batch = lab.sim_clusters(length=4, width=6, nclusters=100, ncrystals=2)
# get more information about the available options for simulating ice
# clusters
help(lab.sim_clusters)
# fit ellipses to the x-, y-, and z-axis views of the cluster and use
# them to calculate aspect ratios
batch.calc_aspect_ratios()
batch.ratios # return a list of the aspect ratios
# plot a histogram of the aspect ratios, passing arguments to
# matplotlib.pyplot.hist()
batch.plot_aspect_ratios(bins=20)
# view the first cluster in the batch
cluster = batch.clusters[0]
cluster.plot()
```
