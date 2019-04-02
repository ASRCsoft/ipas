# ipas
A Python implementation of the [Ice Particle Aggregate Simulator](http://www.carlgschmitt.com/Microphysics.html)

## Installation

```shell
pip install git+https://github.com/ASRCsoft/ipas.git
```

## Crystals (monomers)
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
```
![crystal](https://user-images.githubusercontent.com/4205859/27136311-01852f9a-50e9-11e7-8f10-db348cdddd3a.png)
```python
# project the crystal onto the xy plane, returning a shapely Polygon
crystal.projectxy()
```
![crystal_projection](https://user-images.githubusercontent.com/4205859/27136458-5f9d07ba-50e9-11e7-8665-f230dc932c6a.png)

## Clusters (aggregates)
```python
# use the crystal to start a cluster
cluster = crys.IceCluster(crystal)
# add a new crystal to the cluster
crystal2 = crys.IceCrystal(length=4, width=6)
cluster.add_crystal_from_above(crystal2)
cluster.plot() # same as IceCrystal.plot()
```
![cluster](https://user-images.githubusercontent.com/4205859/27136603-bc31c48e-50e9-11e7-88b1-afe0ba2e5790.png)
```python
# plot the 2D projection of the cluster along with an ellipse fit to
# the area of the projection
cluster.plot_ellipse()
```
![cluster_ellipse](https://user-images.githubusercontent.com/4205859/27136608-bfec5d6e-50e9-11e7-8889-6784203a1937.png)


## Simulations
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
```
![aspect_ratios](https://user-images.githubusercontent.com/4205859/27137095-029e7920-50eb-11e7-989e-f71037e10afb.png)
```python
# view the first cluster in the batch
cluster = batch.clusters[0]
cluster.plot()
```
![cluster0](https://user-images.githubusercontent.com/4205859/27137197-49cb2492-50eb-11e7-9e6c-9b14be9315a5.png)
