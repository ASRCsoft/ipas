"""Utilities for running ice particle simulations."""

import random
import numpy as np
import scipy.optimize as opt
import ipas.crystals as crys

def sim_clusters(length, width, nclusters, ncrystals, reorient='random', rotations=50, speedy=False, lodge=0):
    """Simulate crystal aggregates.

    Args:
        length (float): The column length of the crystals.
        width (float): The width of the hexagonal faces of the crystals.
        nclusters (int): The number of clusters to simulate.
        ncrystals (int): The number of crystals in each cluster.
        reorient (str): The method to use for reorienting crystals and clusters.
            Either 'random' or 'IDL'. Default is 'random'.
        rotations (int): The number of rotations to use to reorient crystals and
            clusters. Default is 50.
        speedy (bool): If true, choose an optimal rotation for single crystals
            instead of reorienting them. Default is false.
        lodge (float): A distance. Used to match the output of IPAS IDL code,
            which uses 0.5. Default is zero.

    Returns:
        An IceClusterBatch object containing the simulated clusters.
    """
    if speedy:
        # get optimal y rotation for single crystals
        f = lambda x: -crys.IceCrystal(length=length, width=width, rotation=[0,x,0]).projectxy().area
        yrot = opt.minimize_scalar(f, bounds=(0, np.pi/2), method='Bounded').x
    clusters = []
    for n in range(nclusters):
        if speedy:
            rotation = [0, yrot, random.uniform(0, 2 * np.pi)]
            seedcrystal = crys.IceCrystal(length=length, width=width, rotation=rotation)
            # how this rotation works:
        
            # x) Zero is always one of the best choices for the x
            # rotation (along with pi/3, 2pi/3, ... ). Since the
            # crystal is symmetric it won't make any difference if we
            # rotate it to the other good values. Might as well just
            # choose zero every time.

            # y) We found the optimal y rotation above.

            # z) Rotation around the z-axis has no affect on the
            # projected area. So just choose randomly, no point in
            # simulating it.

            # Since we have values for all 3 rotations, no need to
            # test multiple rotations. These values also apply to the
            # other new crystals, since they have the same dimensions
            # as the seed crystal. So we don't need to run random
            # rotations for them either.
        else:
            # make the seed crystal, orient it
            seedcrystal = crys.IceCrystal(length=length, width=width)
            seedcrystal.reorient(method=reorient, rotations=rotations)
        # create cluster
        cluster = crys.IceCluster(seedcrystal, size=ncrystals)
        # add new crystals
        nmisses = 0
        while cluster.ncrystals < ncrystals:
            # get the cluster's boundaries
            xmax = cluster.max('x')
            ymax = cluster.max('y')
            xmin = cluster.min('x')
            ymin = cluster.min('y')
            random_loc = [random.uniform(xmin, xmax), random.uniform(ymin, ymax), 0]
            if speedy:
                rotation = [0, yrot, random.uniform(0, 2 * np.pi)]
                seedcrystal = crys.IceCrystal(length=length, width=width, center=random_loc, rotation=rotation)
            else:
                # make a new crystal, orient it
                new_crystal = crys.IceCrystal(length=length, width=width)
                new_crystal.reorient(method=reorient, rotations=rotations)
                new_crystal.move(random_loc)
            # add to the cluster
            crystal_hit = cluster.add_crystal_from_above(new_crystal, lodge=lodge) # returns false if the crystal misses
            # the 'lodge' value is just useful for testing against old IPAS code
            if crystal_hit:
                # recenter the cluster around the center of mass and reorient it
                cluster.recenter()
                cluster.reorient(method=reorient, rotations=rotations)
        clusters.append(cluster)
    plates = width > length
    return IceClusterBatch(clusters, plates)

class IceClusterBatch:
    """A collection of IceCluster objects."""
    
    def __init__(self, clusters, plates=None):
        self.clusters = clusters
        self.plates = plates # are they plates or columns?
        self.ratios = None

    def calc_aspect_ratios(self):
        """Calculate the aspect ratios of the clusters using ellipses fitted
        to the 2D cluster projections from the x-, y-, and z-axis
        perspectives.

        """
        if self.plates is None:
            # inform the user that they need to specify whether these
            # clusters are made of columns or plates
            pass
        if self.plates:
            # if the crystals are plates do the plate version
            ratios = [ cluster.aspect_ratio(method='plate') for cluster in self.clusters ]
        else:
            ratios = [ cluster.aspect_ratio(method='column') for cluster in self.clusters ]
        self.ratios = ratios
        return ratios

    def plot_aspect_ratios(self, **kwargs):
        """Plot a histogram of cluster aspect ratios, sending extra arguments
        to `matplotlib.pyplot.hist`.
        """
        import matplotlib.pyplot as plt
        if self.ratios is None:
            self.calc_aspect_ratios()
        plt.hist(self.ratios, **kwargs)
        plt.ylabel('Frequency')
        if self.plates:
            xlabel = 'Aspect Ratio (Plates)'
        else:
            xlabel = 'Aspect Ratio (Columns)'
        plt.xlabel(xlabel)
