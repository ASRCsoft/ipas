import random
import numpy as np
import ipas.crystals as crys

def sim_clusters(length, width, nclusters, ncrystals, random_rotations=50, max_misses=100, lodge=0):
    clusters = []
    for n in range(nclusters):
        # make the seed crystal, orient it
        seedcrystal = crys.IceCrystal(length=length, width=width)
        seedcrystal.reorient()
        # create cluster
        cluster = crys.IceCluster([seedcrystal])
        # add new crystals
        nmisses = 0
        while len(cluster.crystals) < ncrystals:
            # get the cluster's boundaries
            xmax = max([ crystal.points['x'].max() for crystal in cluster.crystals ])
            ymax = max([ crystal.points['y'].max() for crystal in cluster.crystals ])
            xmin = min([ crystal.points['x'].min() for crystal in cluster.crystals ])
            ymin = min([ crystal.points['y'].min() for crystal in cluster.crystals ])
            # make a new crystal, orient it
            rand_center = [random.uniform(xmin, xmax), random.uniform(ymin, ymax), 0]
            new_crystal = crys.IceCrystal(length=length, width=width, center=rand_center)
            new_crystal.reorient()
            # add to the cluster
            crystal_hit = cluster.add_crystal_from_above(new_crystal, lodge=lodge) # returns false if the crystal misses
            # the 'lodge' value is just useful for testing against old IPAS code
            if crystal_hit:
                # recenter the cluster around the center of mass and reorient it
                cluster.move([ -x for x in cluster.center_of_mass() ])
                #cluster.reorient()
            else: # just in case something goes terribly wrong
                nmisses += 1
                if nmisses > max_misses:
                    break
        clusters.append(cluster)
    plates = width > length
    return Batch(clusters, plates)

class Batch:
    def __init__(self, clusters, plates=None):
        self.clusters = clusters
        self.plates = plates # are they plates or columns?
        self.ratios = None

    def calc_aspect_ratios(self):
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
        import matplotlib.pyplot as plt
        if self.ratios is None:
            self.calc_aspect_ratios()
        plt.hist(self.ratios, *kwargs)
        plt.ylabel('Frequency')
        plt.xlabel('Aspect Ratio')
