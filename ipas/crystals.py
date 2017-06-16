"""Classes representing ice crystals (monomers) and ice clusters
(aggregates).

"""

import numpy as np
import scipy.optimize as opt
import shapely.geometry as geom
import shapely.ops as shops
import shapely.affinity as sha

class IceCrystal:
    """A hexagonal prism representing a single ice crystal."""
    
    def __init__(self, length, width, center=[0, 0, 0], rotation=[0, 0, 0]):
        """Create an ice crystal.


        """
        # put together the hexagonal prism
        ca = length # c axis length
        mf = width # maximum face dimension
        f = np.sqrt(3) / 4 # convenient number for hexagons
        x1 = ca / 2
        self.points = np.array([(x1, -mf / 4, mf * f), (x1, mf / 4, mf * f),
                                (x1, mf / 2, 0), (x1, mf / 4, -mf * f),
                                (x1, -mf / 4, -mf * f), (x1, -mf/2, 0),
                                (-x1, -mf / 4, mf * f), (-x1, mf / 4, mf * f),
                                (-x1, mf / 2, 0), (-x1, mf / 4, -mf * f),
                                (-x1, -mf / 4, -mf * f), (-x1, -mf/2, 0)],
                               dtype=[('x', float), ('y', float), ('z', float)])
        # old IDL code
        # #          1       2      3      4       5      6       7       8       9       10      11      12
        # crystal=[[ca/2.  ,ca/2.  ,ca/2. ,ca/2. ,ca/2.  ,ca/2.  ,-ca/2. ,-ca/2. ,-ca/2. ,-ca/2. ,-ca/2. ,-ca/2.],$
        #          [-mf/4. ,mf/4.  ,mf/2. ,mf/4. ,-mf/4. ,-mf/2. ,-mf/4. ,mf/4.  ,mf/2.  ,mf/4.  ,-mf/4. ,-mf/2.],$
        #          [mf*f   ,mf*f   ,0.    ,-mf*f ,-mf*f  ,0.     ,mf*f   ,mf*f   ,0.     ,-mf*f  ,-mf*f  ,0.]]

        self.center = [0, 0, 0] # start the crystal at the origin
        self._rotate(rotation) # rotate the crystal
        self.rotation = rotation
        self.move(center) # move the crystal
        self.maxz = self.points['z'].max()
        self.minz = self.points['z'].min()
        self.tol = 10 ** -11 # used for some calculations

    def move(self, xyz):
        self.points['x'] += xyz[0]
        self.points['y'] += xyz[1]
        self.points['z'] += xyz[2]
        # update the crystal's center:
        for n in range(3):
            self.center[n] += xyz[n]
        # update max and min
        self.maxz += xyz[2]
        self.minz += xyz[2]

    def _rotate(self, angles):
        # do the new rotation:
        [x, y, z] = [self.points['x'], self.points['y'], self.points['z']]
        [y, z] = [y * np.cos(angles[0]) - z * np.sin(angles[0]), y * np.sin(angles[0]) + z * np.cos(angles[0])]
        [x, z] = [x * np.cos(angles[1]) + z * np.sin(angles[1]), -x * np.sin(angles[1]) + z * np.cos(angles[1])]
        [x, y] = [x * np.cos(angles[2]) - y * np.sin(angles[2]), x * np.sin(angles[2]) + y * np.cos(angles[2])]
        # update the crystal's points:
        self.points['x'] = x
        self.points['y'] = y
        self.points['z'] = z
        # update the crystal's rotation:
        #self.rotation = angles
        # update max and min
        self.maxz = z.max()
        self.minz = z.min()
        # update the crystal's center:
        self.center[0] = (x.max() + x.min()) / 2
        self.center[1] = (y.max() + y.min()) / 2
        self.center[2] = (self.maxz + self.minz) / 2
        # old IDL code:
        # pt1r1=[point(0),point(1)*cos(angle(0))-point(2)*sin(angle(0)),point(1)*sin(angle(0))+point(2)*cos(angle(0))]
        # pt1r2=[pt1r1(0)*cos(angle(1))+pt1r1(2)*sin(angle(1)),pt1r1(1),-pt1r1(0)*sin(angle(1))+pt1r1(2)*cos(angle(1))]
        # pt1r3=[pt1r2(0)*cos(angle(2))-pt1r2(1)*sin(angle(2)),pt1r2(0)*sin(angle(2))+pt1r2(1)*cos(angle(2)),pt1r2(2)]

    def _rev_rotate(self, angles):
        # undo a rotation (same as above but in negative and in reverse order):
        angles = [-x for x in angles ]
        [x, y, z] = [self.points['x'], self.points['y'], self.points['z']]
        [x, y] = [x * np.cos(angles[2]) - y * np.sin(angles[2]), x * np.sin(angles[2]) + y * np.cos(angles[2])]
        [x, z] = [x * np.cos(angles[1]) + z * np.sin(angles[1]), -x * np.sin(angles[1]) + z * np.cos(angles[1])]
        [y, z] = [y * np.cos(angles[0]) - z * np.sin(angles[0]), y * np.sin(angles[0]) + z * np.cos(angles[0])]

        # update the crystal's points:
        self.points['x'] = x
        self.points['y'] = y
        self.points['z'] = z
        # update the crystal's rotation:
        #self.rotation = angles
        # update max and min
        self.maxz = z.max()
        self.minz = z.min()
        # update the crystal's center:
        self.center[0] = (x.max() + x.min()) / 2
        self.center[1] = (y.max() + y.min()) / 2
        self.center[2] = (self.maxz + self.minz) / 2

    def rotate_to(self, angles):
        # rotate to the orientation given by the 3 angles

        # first undo the current rotation
        self._rev_rotate(self.rotation)

        # now do the new rotation
        self._rotate(angles)

        # save the new rotation
        self.rotation = angles

    def reorient(self, method='random', rotations=50):
        if method == 'IDL':
            # based on max_area2.pro from IPAS
            # max_area = self.projectxy().area
            # max_rot1 = None
            max_area = 0
            for i in range(rotations):
                [a, b, c] = [np.random.uniform(high=np.pi), np.random.uniform(high=np.pi), np.random.uniform(high=np.pi)]
                # for mysterious reasons we are going to rotate this 3 times
                rot1 = [a, b, c]
                rot2 = [b * np.pi, c * np.pi, a * np.pi]
                rot3 = [c * np.pi * 2, a * np.pi * 2, b * np.pi * 2]
                self._rotate(rot1)
                self._rotate(rot2)
                self._rotate(rot3)
                new_area = self.projectxy().area
                if new_area > max_area:
                    max_area = new_area
                    max_area = new_area
                    max_rot1 = rot1
                    max_rot2 = rot2
                    max_rot3 = rot3
                # now rotate back -- this is fun!
                self._rev_rotate(rot3)
                self._rev_rotate(rot2)
                self._rev_rotate(rot1)
            # rotate new crystal to the area-maximizing rotation
            self._rotate(max_rot1)
            self._rotate(max_rot2)
            self._rotate(max_rot3)
            self.rotation = [0, 0, 0] # set this new rotation as the default
        elif method == 'random':
            # same as schmitt but only rotating one time, with a real
            # random rotation
            # max_rot = None
            # max_area = self.projectxy().area
            max_area = 0
            for i in range(rotations):
                yrot = np.arccos(np.random.uniform(-1, 1)) - np.pi / 2
                rot = [np.random.uniform(high=2 * np.pi), yrot, np.random.uniform(high=2 * np.pi)]
                self._rotate(rot)
                new_area = self.projectxy().area
                if new_area > max_area:
                    max_area = new_area
                    max_rot = rot
                self._rev_rotate(rot)
            # rotate new crystal to the area-maximizing rotation
            self._rotate(max_rot)
            self.rotation = [0, 0, 0]

    def plot(self):
        # return a multiline object representing the edges of the prism
        lines = []
        hex1 = self.points[0:6]
        hex2 = self.points[6:12]
        # make the lines representing each hexagon
        for hex0 in [hex1, hex2]:
            lines.append(geom.LinearRing(list(hex0)))
        # make the lines connecting the two hexagons
        for n in range(6):
            lines.append(geom.LineString([hex1[n], hex2[n]]))

        return geom.MultiLineString(lines)

    def projectxy(self):
        # project the points onto a 2d surface, return a polygon

        # try using points instead of rotation in the if statements
        p0 = self.points[0]
        p6 = self.points[6]
        if p0['x'] == p6['x'] and p0['y'] == p6['y']:
            # It's vertical, so just return one of the hexagons
            points2d = self.points[0:6]
        else:
            # prism is lying flat or tilted sideways
            midx = self.center[0]
            midy = self.center[1]
            # midx = np.mean(self.points['x'])
            # midy = np.mean(self.points['y'])

            if p0['z'] == p6['z']:
                # It's lying flat (looks like a rectangle from above)
                if len(np.unique(self.points['z'])) == 4:
                    # It's rotated so that there's a ridge on the top,
                    # and the sides are vertical. Ignore the two
                    # highest points (the top ridge), get the next set
                    # of 4 points.
                    points2d = self.points[np.argsort(-self.points['z'])[2:6]]
                else:
                    # find the 4 points farthest from the center
                    distances = np.sqrt((self.points['x'] - midx) ** 2 + (self.points['y'] - midy) ** 2)
                    points2d = self.points[np.argsort(-distances)[0:4]]
            else:
                # prism is tilted. Remove the two highest points from
                # the upper hexagon and the two lowest points from the
                # lower hexagon-- they'll always be inside the
                # resulting 2D octagon!
                hex1 = self.points[0:6]
                hex2 = self.points[6:12]
                if hex1['z'].max() > hex2['z'].max():
                    upperhex = hex1
                    lowerhex = hex2
                else:
                    upperhex = hex2
                    lowerhex = hex1
                upperhex = upperhex[np.argsort(upperhex['z'])[0:4]]
                lowerhex = lowerhex[np.argsort(lowerhex['z'])[2:6]]
                points2d = np.concatenate([upperhex, lowerhex])

            # Get the angle of the line connecting the midpoint (which
            # is always inside the 2d projection) to each point, then
            # sort counterclockwise. This ensures that the points are
            # connected in the right order.
            angles = np.arctan2(points2d['y'] - midy, points2d['x'] - midx)
            points2d = points2d[np.argsort(angles)]
            
        # take away the z-values-- now it's 2D! return polygon
        return geom.Polygon(list(points2d[['x', 'y']]))

    def bottom(self):
        # just return everything for now until I get better code
        # written

        # getting the same points regardless of the orientation
        points = [ geom.Point(list(x)) for x in self.points ]
        lines = []
        faces = []
        
        p0 = self.points[0]
        p6 = self.points[6]
        if abs(p0['x'] - p6['x']) < self.tol and abs(p0['y'] - p6['y']) < self.tol:
            # if it's vertical, only return the hexagon faces
            # (for now)
            for hexagon in range(2):
                n0 = hexagon * 6
                for i in range(5):
                    n = n0 + i
                    lines.append(geom.LineString([self.points[n], self.points[n + 1]]))
                lines.append(geom.LineString([self.points[n0 + 5], self.points[n0]]))
            # get the hexagons only-- no rectangles
            for n in range(2):
                i = n * 6
                faces.append(geom.Polygon(list(self.points[i:(i + 6)])))
        elif abs(p0['z'] - p6['z']) < self.tol:
            # lying flat on its side-- not returning hexagon faces
            if len(np.unique(self.points['z'])) == 4:
                # It's rotated so that there's a ridge on the top, and
                # the sides are vertical. Don't return any vertical
                # rectangular sides
                for n in range(5):
                    p1 = self.points[n]
                    p2 = self.points[n + 1]
                    # is it a non-vertical rectangle?
                    if abs(p1['x'] - p2['x']) >= self.tol and abs(p1['y'] - p2['y']) >= self.tol:
                        faces.append(geom.Polygon([self.points[n], self.points[n + 1],
                                                   self.points[n + 7], self.points[n + 6]]))
                # get that last rectangle I missed
                p1 = self.points[5]
                p2 = self.points[0]
                if abs(p1['x'] - p2['x']) >= self.tol and abs(p1['y'] - p2['y']) >= self.tol:
                    faces.append(geom.Polygon([self.points[5], self.points[0],
                                               self.points[6], self.points[11]]))
                # get the lines around the hexagons
                for hexagon in range(2):
                    n0 = hexagon * 6
                    for i in range(5):
                        n = n0 + i
                        p1 = self.points[n]
                        p2 = self.points[n + 1]
                        if abs(p1['x'] - p2['x']) >= self.tol and abs(p1['y'] - p2['y']) >= self.tol:
                            lines.append(geom.LineString([self.points[n], self.points[n + 1]]))
                    p1 = self.points[n0 + 5]
                    p2 = self.points[n0]
                    if abs(p1['x'] - p2['x']) >= self.tol and abs(p1['y'] - p2['y']) >= self.tol:
                        lines.append(geom.LineString([self.points[n0 + 5], self.points[n0]]))
                # get the between-hexagon lines
                for n in range(6):
                    lines.append(geom.LineString([self.points[n], self.points[n + 6]]))
                
                
            # returning only rectangles
            pass
        else:
            # return all the faces

            # get the lines around the hexagons
            for hexagon in range(2):
                n0 = hexagon * 6
                for i in range(5):
                    n = n0 + i
                    lines.append(geom.LineString([self.points[n], self.points[n + 1]]))
                lines.append(geom.LineString([self.points[n0 + 5], self.points[n0]]))
            # get the between-hexagon lines
            for n in range(6):
                lines.append(geom.LineString([self.points[n], self.points[n + 6]]))
            # get the hexagons
            for n in range(2):
                i = n * 6
                faces.append(geom.Polygon(list(self.points[i:(i + 6)])))
            # get the rectangles
            for n in range(5):
                faces.append(geom.Polygon([self.points[n], self.points[n + 1],
                                           self.points[n + 7], self.points[n + 6]]))
            # get that last rectangle I missed
            faces.append(geom.Polygon([self.points[5], self.points[0],
                                       self.points[6], self.points[11]]))
        
        # return the geometry representing the bottom side of the prism

        # # similar to projectxy
        # if self.rotation[1] == math.pi / 2:
        #     # it's vertical, so just return one of the hexagons
        #     points = self.points[0:6]

        # first find top and bottom hexagon

        # remove the top two points

        # make the lines

        # make the faces

        return {'lines': lines, 'points': points, 'faces': faces}

    def top(self):
        # return the geometry representing the top side of the prism

        # first find top and bottom hexagon

        # remove the bottom two points

        # make the lines

        # make the faces

        #return {'lines': lines, 'points': points, 'faces': faces}

        # temporary, until I fix these functions
        top = self.bottom()
        # # unless it's vertical, that's easy
        # if self.rotation[1] / (np.pi / 2) % 4 == 1:
        #     top['points'] = [ geom.Point(list(x)) for x in self.points[0:6] ]
        #     top['lines'] = []
        #     for i in range(5): # get the points around each hexagon
        #         top['lines'].append(geom.LineString([self.points[i], self.points[i + 1]]))
        #     top['lines'].append(geom.LineString([self.points[5], self.points[0]]))
        # elif self.rotation[1] / (np.pi / 2) % 4 == 3:
        #     top['points'] = [ geom.Point(list(x)) for x in self.points[6:12] ]
        #     top['lines'] = []
        #     for i in range(5): # get the points around each hexagon
        #         top['lines'].append(geom.LineString([self.points[i + 6], self.points[i + 7]]))
        #         top['lines'].append(geom.LineString([self.points[11], self.points[6]]))
        
        return top

    def min_vert_dist(self, crystal2):
        # find the minimum directed distance to crystal2 traveling straight downward

        rel_area = self.projectxy().intersection(crystal2.projectxy())
        if not isinstance(rel_area, geom.Polygon):
            return None
        c1_bottom = self.bottom()
        c2_top = crystal2.top()
        mindiffz = self.maxz - crystal2.minz

        # 1) lines and lines
        # all the intersections are calculated in 2d so no need to
        # convert these 3d objects!
        c1_lines = [ l for l in c1_bottom['lines'] if l.intersects(rel_area) ]
        c2_lines = [ l for l in c2_top['lines'] if l.intersects(rel_area) ]
        for line1 in c1_lines:
            for line2 in c2_lines:
                if line1.intersects(line2):
                    # get (2D) point of intersection
                    xy = line1.intersection(line2)
                    if not isinstance(xy, geom.point.Point):
                        # parallel lines don't count
                        continue
                    # get z difference
                    # make sure the damn lines aren't vertical
                    xrange1 = line1.xy[0][1] - line1.xy[0][0]
                    xrange2 = line2.xy[0][1] - line2.xy[0][0]
                    if xrange1 != 0:
                        # interpolate using x value
                        z1 = line1.interpolate((xy.x - line1.xy[0][0]) / (xrange1), normalized=True).z
                    else:
                        # interpolate using y value
                        z1 = line1.interpolate((xy.y - line1.xy[1][0]) / (line1.xy[1][1] - line1.xy[1][0]), normalized=True).z
                    if xrange2 != 0:
                        z2 = line2.interpolate((xy.x - line2.xy[0][0]) / (xrange2), normalized=True).z
                    else:
                        z2 = line2.interpolate((xy.y - line2.xy[1][0]) / (line2.xy[1][1] - line2.xy[1][0]), normalized=True).z
                    diffz = z1 - z2
                    if diffz < mindiffz:
                        mindiffz = diffz
        
        # 2) points and surfaces
        c1_points = [ p for p in c1_bottom['points'] if p.intersects(rel_area) ]
        c2_faces = [ f for f in c2_top['faces'] if f.intersects(rel_area) ]
        for point in c1_points:
            for face in c2_faces:
                if point.intersects(face):
                    # get z difference
                    z1 = point.z
                    # find the equation of the polygon's plane, plug in xy
                    a = np.array(face.exterior.coords[0])
                    AB = np.array(face.exterior.coords[1]) - a
                    AC = np.array(face.exterior.coords[2]) - a
                    normal_vec = np.cross(AB, AC)
                    # find constant value
                    d = -np.dot(normal_vec, a)
                    z2 = -(point.x * normal_vec[0] + point.y * normal_vec[1] + d) / normal_vec[2]
                    diffz = z1 - z2
                    if diffz < mindiffz:
                        mindiffz = diffz
                    # the point can only intersect one face, so we're
                    # done with this one
                    #break
                    # ^ I should be able to do that but I have to fix my 'bottom' function first!
        
        # 3) surfaces and points
        c1_faces = [ f for f in c1_bottom['faces'] if f.intersects(rel_area) ]
        c2_points = [ p for p in c2_top['points'] if p.intersects(rel_area) ]
        for point in c2_points:
            for face in c1_faces:
                if point.intersects(face):
                    # get z difference
                    z2 = point.z # z2 this time!!!
                    # find the equation of the polygon's plane, plug in xy
                    a = np.array(face.exterior.coords[0])
                    AB = np.array(face.exterior.coords[1]) - a
                    AC = np.array(face.exterior.coords[2]) - a
                    normal_vec = np.cross(AB, AC)
                    # find constant value
                    d = -np.dot(normal_vec, a)
                    z1 = -(point.x * normal_vec[0] + point.y * normal_vec[1] + d) / normal_vec[2]
                    diffz = z1 - z2
                    if diffz < mindiffz:
                        mindiffz = diffz
                        # the point can only intersect one face, so we're
                        # done with this one
                    #break

        return mindiffz

    def write_obj(self, filename):
        f = open(filename, 'w')
        # write the vertices
        for n in range(12):
            f.write('v ' + ' '.join(map(str, self.points[n])) + '\n')
        # write the hexagons
        for n in range(2):
            f.write('f ' + ' '.join(map(str, range(n * 6 + 1, (n + 1) * 6 + 1))) + '\n')
        for n in range(5):
            f.write('f ' + ' '.join(map(str, [n + 1, n + 2, n + 8, n + 7])) + '\n')
        f.write('f ' + ' '.join(map(str, [6, 1, 7, 12])) + '\n')
        f.close()

        
class IceCluster:
    def __init__(self, crystal, size=1):
        # needed for bookkeeping:
        self.rotation = [0, 0, 0]
        self.points = np.full((size, 12), np.nan,
                              dtype=[('x', float), ('y', float), ('z', float)])
        self.points[0] = crystal.points
        self.size = size
        self.ncrystals = 1
        # used for some calculations involving shapely objects
        self.tol = 10 ** -11
        # Used for the fit_ellipse function. I really do not like that
        # I have to set this so high, arrr.
        self.tol_ellipse = 10 ** -4.5
        # to store axis lengths
        self.major_axis = {}
        self.minor_axis = {}

    def crystals(self, i=None):
        # return a crystal with the same points and attributes as the
        # nth crystal in the cluster
        if i is None:
            crystals = []
            for n in range(self.ncrystals):
                cr = IceCrystal(1, 1)
                cr.points = self.points[n]
                cr.rotation = self.rotation
                cx = cr.points['x'].mean()
                cy = cr.points['y'].mean()
                cz = cr.points['z'].mean()
                cr.center = [cx, cy, cz]
                cr.maxz = cr.points['z'].max()
                cr.minz = cr.points['z'].min()
                crystals.append(cr)
            return crystals
        else:
            cr = IceCrystal(1, 1)
            cr.points = self.points[i]
            cr.rotation = self.rotation
            cx = cr.points['x'].mean()
            cy = cr.points['y'].mean()
            cz = cr.points['z'].mean()
            cr.center = [cx, cy, cz]
            cr.maxz = cr.points['z'].max()
            cr.minz = cr.points['z'].min()
            return cr

    def _add_crystal(self, crystal):
        n = self.ncrystals
        if n < self.size:
            self.points[n] = crystal.points
        else:
            self.points = np.append(self.points, [crystal.points], axis=0)
        self.ncrystals += 1

    def move(self, xyz):
        # move the entire cluster
        self.points['x'][:self.ncrystals] += xyz[0]
        self.points['y'][:self.ncrystals] += xyz[1]
        self.points['z'][:self.ncrystals] += xyz[2]

    def max(self, dim):
        return self.points[dim][:self.ncrystals].max()

    def min(self, dim):
        return self.points[dim][:self.ncrystals].min()

    def _rotate(self, angles):
        [x, y, z] = [self.points[:self.ncrystals]['x'], self.points[:self.ncrystals]['y'], self.points[:self.ncrystals]['z']]
        [y, z] = [y * np.cos(angles[0]) - z * np.sin(angles[0]), y * np.sin(angles[0]) + z * np.cos(angles[0])]
        [x, z] = [x * np.cos(angles[1]) + z * np.sin(angles[1]), -x * np.sin(angles[1]) + z * np.cos(angles[1])]
        [x, y] = [x * np.cos(angles[2]) - y * np.sin(angles[2]), x * np.sin(angles[2]) + y * np.cos(angles[2])]
        # update the crystal's points:
        self.points['x'][:self.ncrystals] = x
        self.points['y'][:self.ncrystals] = y
        self.points['z'][:self.ncrystals] = z

    def _rev_rotate(self, angles):
        angles = [-x for x in angles ]
        [x, y, z] = [self.points[:self.ncrystals]['x'], self.points[:self.ncrystals]['y'], self.points[:self.ncrystals]['z']]
        [x, y] = [x * np.cos(angles[2]) - y * np.sin(angles[2]), x * np.sin(angles[2]) + y * np.cos(angles[2])]
        [x, z] = [x * np.cos(angles[1]) + z * np.sin(angles[1]), -x * np.sin(angles[1]) + z * np.cos(angles[1])]
        [y, z] = [y * np.cos(angles[0]) - z * np.sin(angles[0]), y * np.sin(angles[0]) + z * np.cos(angles[0])]
        # update the crystal's points:
        self.points['x'][:self.ncrystals] = x
        self.points['y'][:self.ncrystals] = y
        self.points['z'][:self.ncrystals] = z

    def rotate_to(self, angles):
        # rotate the entire cluster

        # first get back to the original rotation
        if any(np.array(self.rotation) != 0):
            self._rev_rotate(self.rotation)

        # now add the new rotation
        self._rotate(angles)

        # save the new rotation
        self.rotation = angles

    def center_of_mass(self):
        x = np.mean(self.points[:self.ncrystals]['x'])
        y = np.mean(self.points[:self.ncrystals]['y'])
        z = np.mean(self.points[:self.ncrystals]['z'])
        return [x, y, z]

    def recenter(self):
        self.move([ -x for x in self.center_of_mass() ])

    def plot(self):
        return geom.MultiLineString([ lines for crystal in self.crystals() for lines in crystal.plot() ])

    def _crystal_projectxy(self, n):
        # try using points instead of rotation in the if statements
        points = self.points[n]
        p0 = points[0]
        p6 = points[6]
        # If differences are less than 'tol' we just consider them
        # equal. I don't know why but these small values trip up
        # shapely.
        if abs(p0['x'] - p6['x']) < self.tol and abs(p0['y'] - p6['y']) < self.tol:
            # It's vertical, so just return one of the hexagons
            points2d = points[0:6]
        else:
            # prism is lying flat or tilted sideways
            midx = points['x'].mean()
            midy = points['y'].mean()
            # midx = np.mean(self.points['x'])
            # midy = np.mean(self.points['y'])

            if abs(p0['z'] - p6['z']) < self.tol:
                # It's lying flat (looks like a rectangle from above)
                if len(np.unique(points['z'])) == 4:
                    # It's rotated so that there's a ridge on the top,
                    # and the sides are vertical. Ignore the two
                    # highest points (the top ridge), get the next set
                    # of 4 points.
                    points2d = points[np.argsort(-points['z'])[2:6]]
                else:
                    # find the 4 points farthest from the center
                    distances = np.sqrt((points['x'] - midx) ** 2 + (points['y'] - midy) ** 2)
                    points2d = points[np.argsort(-distances)[0:4]]
            else:
                # prism is tilted. Remove the two highest points from
                # the upper hexagon and the two lowest points from the
                # lower hexagon-- they'll always be inside the
                # resulting 2D octagon!
                hex1 = points[0:6]
                hex2 = points[6:12]
                if hex1['z'].max() > hex2['z'].max():
                    upperhex = hex1
                    lowerhex = hex2
                else:
                    upperhex = hex2
                    lowerhex = hex1
                upperhex = upperhex[np.argsort(upperhex['z'])[0:4]]
                lowerhex = lowerhex[np.argsort(lowerhex['z'])[2:6]]
                points2d = np.concatenate([upperhex, lowerhex])

            # Get the angle of the line connecting the midpoint (which
            # is always inside the 2d projection) to each point, then
            # sort counterclockwise. This ensures that the points are
            # connected in the right order.
            angles = np.arctan2(points2d['y'] - midy, points2d['x'] - midx)
            points2d = points2d[np.argsort(angles)]
            
        # take away the z-values-- now it's 2D! return polygon
        return geom.Polygon(list(points2d[['x', 'y']]))

    def projectxy(self):
        polygons = [ self._crystal_projectxy(n) for n in range(self.ncrystals) ]
        return shops.cascaded_union(polygons)

    def add_crystal_from_above(self, crystal, lodge=0):
        # drop a new crystal onto the cluster

        # first see which other crystals intersect with the new one in
        # the xy plane
        newpoly = crystal.projectxy()
        # use the bounding box to determine which crystals to get
        xmax = max(newpoly.exterior.coords.xy[0])
        ymax = max(newpoly.exterior.coords.xy[1])
        xmin = min(newpoly.exterior.coords.xy[0])
        ymin = min(newpoly.exterior.coords.xy[1])
        close = np.all([self.points['x'][:self.ncrystals].max(axis=1) > xmin,
                        self.points['x'][:self.ncrystals].min(axis=1) < xmax,
                        self.points['y'][:self.ncrystals].max(axis=1) > ymin,
                        self.points['y'][:self.ncrystals].min(axis=1) < ymax], axis=0)
        which_close = np.where(close)
        close_crystals = [ self.crystals(n) for n in which_close[0] ]
        # see which crystals could actually intersect with the new crystal
        close_crystals = [ x for x in close_crystals if x.projectxy().intersects(newpoly) ]
        # close_crystals = [ x for x in self.crystals() if x.projectxy().intersects(newpoly) ]
        if len(close_crystals) == 0:
            return False # the crystal missed!
        
        # we know highest hit is >= max(minzs), therefore the first
        # hit can't be below (max(minzs) - height(crystal))
        minzs = [ crystal2.minz for crystal2 in close_crystals ]
        first_hit_lower_bound = max(minzs) - (crystal.maxz - crystal.minz)
        # remove the low crystals, sort from highest to lowest
        close_crystals = [ x for x in close_crystals if x.maxz > first_hit_lower_bound ]
        close_crystals.sort(key=lambda x: x.maxz, reverse=True)

        # look to see where the new crystal hits the old ones
        mindiffz = crystal.minz - first_hit_lower_bound # the largest it can possibly be
        for crystal2 in close_crystals:
            if first_hit_lower_bound > crystal2.maxz:
                break # stop looping if the rest of the crystals are too low
            diffz = crystal.min_vert_dist(crystal2)
            #return diffz
            # update if needed
            if diffz < mindiffz:
                mindiffz = diffz
                first_hit_lower_bound = crystal.minz - mindiffz
                # take the highest hit, move the crystal to that level
        crystal.move([0, 0, -mindiffz - lodge])

        # append new crystal to list of crystals
        # self.crystals.append(crystal)
        self._add_crystal(crystal)
        
        # fin.
        return True

    def reorient(self, method='random', rotations=50):
        if method == 'IDL':
            # based on max_agg3.pro from IPAS
            max_area = 0
            for i in range(rotations):
                [a, b, c] = [np.random.uniform(high=np.pi / 4), np.random.uniform(high=np.pi / 4), np.random.uniform(high=np.pi / 4)]
                # for mysterious reasons we are going to rotate this 3 times
                rot1 = [a, b, c]
                rot2 = [b * np.pi, c * np.pi, a * np.pi]
                rot3 = [c * np.pi * 2, a * np.pi * 2, b * np.pi * 2]
                self._rotate(rot1)
                self._rotate(rot2)
                self._rotate(rot3)
                new_area = self.projectxy().area
                if new_area > max_area:
                    max_area = new_area
                    max_rot1 = rot1
                    max_rot2 = rot2
                    max_rot3 = rot3
                # now rotate back -- this is fun!
                self._rev_rotate(rot3)
                self._rev_rotate(rot2)
                self._rev_rotate(rot1)
            # rotate new crystal to the area-maximizing rotation(s)
            self._rotate(max_rot1)
            self._rotate(max_rot2)
            self._rotate(max_rot3)
            self.rotation = [0, 0, 0]
            
        elif method == 'random':
            # same as schmitt but only rotating one time, with a real
            # random rotation
            max_area = 0
            for i in range(rotations):
                yrot = np.arccos(np.random.uniform(-1, 1)) - np.pi/2
                rot = [np.random.uniform(high=2 * np.pi), yrot, np.random.uniform(high=2 * np.pi)]
                self._rotate(rot)
                new_area = self.projectxy().area
                if new_area > max_area:
                    max_area = new_area
                    max_rot = rot
                self._rev_rotate(rot)
            # rotate new crystal to the area-maximizing rotation
            self._rotate(max_rot)
            self.rotation = [0, 0, 0]
        
        # elif method == 'bh':
        #     # use a basin-hopping algorithm to look for the optimal rotation
        #     def f(x):
        #         # yrot = np.arccos(x[1]) - np.pi/2
        #         # self.rotate_to([x[0], yrot, 0])
        #         self.rotate_to([x[0], x[1], 0])
        #         return -self.projectxy().area
        #     # lbfgsb_opt = {'ftol': 1, 'maxiter': 5}
        #     # min_kwargs = {'bounds': [(0, np.pi), (0, np.pi)], 'options': lbfgsb_opt}
        #     # # min_kwargs = {'bounds': [(0, np.pi), (0, np.pi)]}
        #     # opt_rot = opt.basinhopping(f, x0=[np.pi/2, np.pi/2], niter=15, stepsize=np.pi / 7,
        #     #                            interval=5, minimizer_kwargs=min_kwargs)
        #     lbfgsb_opt = {'ftol': 1, 'maxiter': 0, 'maxfun': 4}
        #     min_kwargs = {'bounds': [(0, np.pi), (-0.99, 0.99)], 'options': lbfgsb_opt}
        #     # min_kwargs = {'bounds': [(0, np.pi), 0, np.pi)]}
        #     opt_rot = opt.basinhopping(f, x0=[np.pi/2, np.pi/2], niter=30, stepsize=np.pi / 4,
        #                                interval=10, minimizer_kwargs=min_kwargs)
        #     # xrot = opt_rot.x[0]
        #     # yrot = np.arccos(opt_rot.x[1]) - np.pi / 2
        #     [xrot, yrot] = opt_rot.x
        #     # area at rotation + pi is the same, so randomly choose to
        #     # add those
        #     # if np.random.uniform() > .5:
        #     #     xrot += np.pi
        #     # if np.random.uniform() > .5:
        #     #     yrot += np.pi
        #     zrot = np.random.uniform(high=2 * np.pi) # randomly choose z rotation
        #     self.rotate_to([xrot, yrot, zrot])
        #     self.rotation = [0, 0, 0]
        #     return opt_rot

        # elif method == 'diff_ev':
        #     def f(x):
        #         # yrot = np.arccos(x[1]) - np.pi/2
        #         # self.rotate_to([x[0], yrot, 0])
        #         self.rotate_to([x[0], x[1], 0])
        #         return -self.projectxy().area
        #     opt_rot = opt.differential_evolution(f, [(0, np.pi), (-1, 1)],
        #                                          maxiter=10, popsize=15)
        #     # xrot = opt_rot.x[0]
        #     # yrot = np.arccos(opt_rot.x[1]) - np.pi / 2
        #     [xrot, yrot] = opt_rot.x
        #     zrot = np.random.uniform(high=2 * np.pi) # randomly choose z rotation
        #     self.rotate_to([xrot, yrot, zrot])
        #     self.rotation = [0, 0, 0]
        #     return opt_rot

    def _get_moments(self, poly):
        # get 'mass moments' for this cluster's 2D polygon using a
        # variation of the shoelace algorithm
        xys = poly.exterior.coords.xy
        npoints = len(xys[0])
        # values for the three points-- point[n], point[n+1], and
        # (0,0)-- making up triangular slices from the origin to the
        # edges of the polygon
        xmat = np.array([xys[0][0:-1], xys[0][1:], np.zeros(npoints - 1)]).transpose()
        ymat = np.array([xys[1][0:-1], xys[1][1:], np.zeros(npoints - 1)]).transpose()
        # arrange the points in left-center-right order
        x_order = np.argsort(xmat, axis=1)
        ordered_xmat = xmat[np.array([range(npoints - 1)]).transpose(), x_order]
        ordered_ymat = ymat[np.array([range(npoints - 1)]).transpose(), x_order]
        xl = ordered_xmat[:, 0]
        xm = ordered_xmat[:, 1]
        xr = ordered_xmat[:, 2]
        yl = ordered_ymat[:, 0]
        ym = ordered_ymat[:, 1]
        yr = ordered_ymat[:, 2]
        # which slices have areas on the left and right sides of the
        # middle point? Ignore values smaller than 'tol' so we don't
        # run into terrible problems with division.
        left = xm - xl > self.tol_ellipse
        right = xr - xm > self.tol_ellipse
        # slope and intercept of line connecting left and right points
        has_area = xr != xl
        m3 = np.zeros(npoints - 1)
        m3[has_area] = (yr[has_area] - yl[has_area]) / (xr[has_area] - xl[has_area])
        b3 = -xl * m3 + yl
        # the y coordinate of the line connecting the left and right
        # points at the x position of the middle point
        m3_mid = yl + m3 * (xm - xl)
        # is the midpoint above or below that line?
        mid_below = ym < m3_mid
        # line connecting left and middle point (where applicable)
        m1 = (ym[left] - yl[left]) / (xm[left] - xl[left])
        b1 = -xl[left] * m1 + yl[left]
        # line connecting middle and right point (where applicable)
        m2 = (yr[right] - ym[right]) / (xr[right] - xm[right])
        b2 = -xr[right] * m2 + yr[right]
        # now that we have the points in a nice format + helpful
        # information we can calculate the integrals of the slices
        xx = np.zeros(npoints - 1)
        xy = np.zeros(npoints - 1)
        yy = np.zeros(npoints - 1)
        dxl = (xm[left] - xl[left])
        dx2l = (xm[left] ** 2 - xl[left] ** 2)
        dx3l = (xm[left] ** 3 - xl[left] ** 3)
        dx4l = (xm[left] ** 4 - xl[left] ** 4)
        dxr = (xr[right] - xm[right])
        dx2r = (xr[right] ** 2 - xm[right] ** 2)
        dx3r = (xr[right] ** 3 - xm[right] ** 3)
        dx4r = (xr[right] ** 4 - xm[right] ** 4)
        # x^2
        xx[left] = dx4l * (m1 - m3[left]) / 4 +\
                   dx3l * (b1 - b3[left]) / 3
        xx[right] += dx4r * (m2 - m3[right]) / 4 +\
                     dx3r * (b2 - b3[right]) / 3
        # x*y
        xy[left] = dx4l * (m1 ** 2 - m3[left] ** 2) / 8 +\
                   dx3l * (b1 * m1 - b3[left] * m3[left]) / 3 +\
                   dx2l * (b1 ** 2 - b3[left] ** 2) / 4
        xy[right] += dx4r * (m2 ** 2 - m3[right] ** 2) / 8 +\
                     dx3r * (b2 * m2 - b3[right] * m3[right]) / 3 +\
                     dx2r * (b2 ** 2 - b3[right] ** 2) / 4
        # y^2
        yy[left] = dx4l * (m1 ** 3 - m3[left] ** 3) / 12 +\
                   dx3l * (b1 * m1 ** 2 - b3[left] * m3[left] ** 2) / 3 +\
                   dx2l * (b1 ** 2 * m1 - b3[left] ** 2 * m3[left]) / 2 +\
                   dxl * (b1 ** 3 - b3[left] ** 3) / 3
        yy[right] += dx4r * (m2 ** 3 - m3[right] ** 3) / 12 +\
                     dx3r * (b2 * m2 ** 2- b3[right] * m3[right] ** 2) / 3 +\
                     dx2r * (b2 ** 2 * m2 - b3[right] ** 2 * m3[right]) / 2 +\
                     dxr * (b2 ** 3 - b3[right] ** 3) / 3
        # if the middle point was below the other points, multiply by
        # minus 1
        xx[mid_below] *= -1
        xy[mid_below] *= -1
        yy[mid_below] *= -1
        # find out which slices were going clockwise, and make those
        # negative        
        points = np.array([xys[0], xys[1]]).transpose()
        cross_prods = np.cross(points[:-1], points[1:])
        clockwise = cross_prods < 0
        xx[clockwise] *= -1
        xy[clockwise] *= -1
        yy[clockwise] *= -1
        # add up the totals across the entire polygon
        xxtotal = np.sum(xx)
        yytotal = np.sum(yy)
        xytotal = np.sum(xy)
        # and if the points were in clockwise order, flip the sign
        if np.sum(cross_prods) < 0:
            xxtotal *= -1
            yytotal *= -1
            xytotal *= -1
        # also need to account for the holes, if they exist
        for linestring in list(poly.interiors):
            hole = geom.Polygon(linestring)
            hole_moments = self._get_moments(hole)
            xxtotal -= hole_moments[0]
            yytotal -= hole_moments[1]
            xytotal -= hole_moments[2]
        return [xxtotal, yytotal, xytotal]

    def fit_ellipse(self):
        # Emulating this function, but for polygons in continuous
        # space rather than blobs in discrete space:
        # http://www.idlcoyote.com/ip_tips/fit_ellipse.html
        poly = self.projectxy()
        xy_area = poly.area
        
        # center the polygon around the centroid
        centroid = poly.centroid
        poly = sha.translate(poly, -centroid.x, -centroid.y)

        # occasionally we get multipolygons
        if isinstance(poly, geom.MultiPolygon):
            xx = 0
            yy = 0
            xy = 0
            for poly2 in poly:
                moments = self._get_moments(poly2)
                xx += moments[0] / xy_area
                yy += moments[1] / xy_area
                xy -= moments[2] / xy_area
        else:
            moments = self._get_moments(poly)
            xx = moments[0] / xy_area
            yy = moments[1] / xy_area
            xy = -moments[2] / xy_area

        # do magic
        m = np.matrix([[yy, xy], [xy, xx]])
        evals, evecs = np.linalg.eigh(m)
        semimajor = np.sqrt(evals[0]) * 2
        semiminor = np.sqrt(evals[1]) * 2
        major = semimajor * 2
        minor = semiminor * 2
        # semiAxes = [semimajor, semiminor]
        # axes = [major, minor]
        evec = np.squeeze(np.asarray(evecs[0]))
        orientation = np.arctan2(evec[1], evec[0]) * 180 / np.pi

        ellipse = {'xy': [centroid.x, centroid.y], 'width': minor,
                   'height': major, 'angle': orientation}
        return ellipse

    def plot_ellipse(self):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse
        # get the ellipse
        params = self.fit_ellipse()
        ellipse = Ellipse(**params)
        # get the polygon
        poly = self.projectxy()
        fig = plt.figure(0)
        ax = fig.add_subplot(111)
        ax.add_artist(ellipse)
        ellipse.set_alpha(.5)
        if isinstance(poly, geom.multipolygon.MultiPolygon):
            for poly2 in poly:
                x, y = poly2.exterior.xy
                ax.plot(x, y)
        else:
            x, y = poly.exterior.xy
            ax.plot(x, y)
        #maxdim = max([params['width'], params['height']])# / 2
        #ax.set_xlim([-maxdim + params['xy'][0], maxdim + params['xy'][0]])
        #ax.set_ylim([-maxdim + params['xy'][1], maxdim + params['xy'][1]])
        ax.set_aspect('equal', 'datalim')

    def write_obj(self, filename):
        f = open(filename, 'w')

        faces = []
        for i, crystal in enumerate(self.crystals()):
            nc = i * 12
            # write the vertices
            for n in range(12):
                f.write('v ' + ' '.join(map(str, crystal.points[n])) + '\n')
            # write the hexagons
            for n in range(2):
                coords = range(n * 6 + 1 + nc, (n + 1) * 6 + 1 + nc)
                faces.append('f ' + ' '.join(map(str, coords)))
            # write the rectangles
            for n in range(5):
                coords = [n + 1 + nc, n + 2 + nc, n + 8 + nc, n + 7 + nc]
                faces.append('f ' + ' '.join(map(str, coords)))
            # write the last rectangle I missed
            coords = [nc + 6, nc + 1, nc + 7, nc + 12]
            faces.append('f ' + ' '.join(map(str, coords)))
        f.write('\n'.join(faces))
        f.close()

    def aspect_ratio(self, method):
        rotation = self.rotation
        
        # getting ellipse axes from 3 perspectives
        ellipse = {}
        self.rotate_to([0, 0, 0])
        ellipse['z'] = self.fit_ellipse()
        self.rotate_to([np.pi / 2, 0, 0])
        ellipse['y'] = self.fit_ellipse()
        self.rotate_to([np.pi / 2, np.pi / 2, 0])
        ellipse['x'] = self.fit_ellipse()
        
        # put the cluster back
        self.rotate_to(rotation)

        for dim in ellipse.keys():
            self.major_axis[dim] = max(ellipse[dim]['height'], ellipse[dim]['width'])
            self.minor_axis[dim] = min(ellipse[dim]['height'], ellipse[dim]['width'])

        if method == 1:
            return max(self.major_axis.values()) / max(self.minor_axis.values())
        elif method == 'plate':
            return self.minor_axis['y'] / self.major_axis['z']
        elif method == 'column':
            return self.major_axis['z'] / self.minor_axis['y']
