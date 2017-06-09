# do crystal things!

import random
import numpy as np
import shapely.geometry as geom
import shapely.ops as shops
import shapely.affinity as sha

class IceCrystal:
    def __init__(self, length, width, center=[0, 0, 0], rotation=[0, 0, 0]):
        self.length = length
        self.width = width
        self.center = [0, 0, 0] # start the crystal at the origin

        # put together the hexagonal prism
        ca = self.length # c axis length
        mf = self.width # maximum face dimension
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

        self._rotate(rotation) # rotate the crystal
        self.rotation = rotation
        self.move(center) # move the crystal
        self.maxz = self.points['z'].max()
        self.minz = self.points['z'].min()

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

    def reorient(self, method='schmitt', random_rotations=50):
        if method == 'schmitt':
            max_area = self.projectxy().area
            max_rot1 = None
            for i in range(random_rotations):
                [a, b, c] = [random.uniform(0, np.pi), random.uniform(0, np.pi), random.uniform(0, np.pi)]
                # for mysterious reasons we are going to rotate this 3 times
                rot1 = [a, b, c]
                rot2 = [b * np.pi, c * np.pi, a * np.pi]
                rot3 = [c * np.pi / 2, a * np.pi / 2, b * np.pi / 2]
                self.rotate_to(rot1)
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
            # rotate new crystal to the area-maximizing rotation(s)
            if max_rot1 is not None:
                self.rotate_to(max_rot1)
                self._rotate(max_rot2)
                self._rotate(max_rot3)
                # if none of the new rotations was better we can leave
                # it at the original rotation
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
        if self.points['z'][0] == self.points['z'][6]:
        # if (self.rotation[1] / (np.pi / 2)) % 2 == 1:
            # (odd multiple of pi/2) it's vertical, so just return one
            # of the hexagons
            points2d = self.points[0:6]
        else:
            # prism is lying flat or tilted sideways
            midx = self.center[0]
            midy = self.center[1]
            # midx = np.mean(self.points['x'])
            # midy = np.mean(self.points['y'])

            if self.points['y'][0] == self.points['y'][3]:
            # if self.rotation[1] == 0:
                # it's lying flat (looks like a rectangle from above)
                if len(np.unique(self.points['z'])) == 3:
                # if (self.rotation[0] / (np.pi / 6)) % 2 == 1:
                    # (odd multiple of pi/6) It's rotated so that
                    # there's a ridge on the top, and the sides are
                    # vertical. Ignore the two highest points (the top
                    # ridge), get the next set of 4 points.
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

        points = []
        lines = []
        faces = []
        # # unless it's vertical, that's easy
        # if self.points['x'][0] == self.points['x'][6] and self.points['y'][0] == self.points['y'][6]:
        # # if (self.rotation[1] / (np.pi / 2)) % 4 == 1:
        #     points = [ geom.Point(list(x)) for x in self.points[6:12] ]
        #     for i in range(5): # get the points around each hexagon
        #         lines.append(geom.LineString([self.points[i + 6], self.points[i + 7]]))
        #         lines.append(geom.LineString([self.points[11], self.points[6]]))
        # elif (self.rotation[1] / (np.pi / 2)) % 4 == 3:
        #     points = [ geom.Point(list(x)) for x in self.points[0:6] ]
        #     for i in range(5): # get the points around each hexagon
        #         lines.append(geom.LineString([self.points[i], self.points[i + 1]]))
        #         lines.append(geom.LineString([self.points[5], self.points[0]]))
        # else:
        
        # get the points
        for n in range(12):
            points.append(geom.Point(list(self.points[n])))
        # get the hexagon lines
        for hexn in range(2): # cycle through 2 hexagons
            for n in range(5): # get the points around each hexagon
                i = n + hexn * 6
                lines.append(geom.LineString([self.points[i], self.points[i + 1]]))
            # get that last line I missed
            i = hexn * 6
            lines.append(geom.LineString([self.points[i + 5], self.points[i]]))
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


# an ice cluster!!
class IceCluster:
    def __init__(self, crystals, size=1):
        self.crystals = crystals
        self.ncrystals = len(self.crystals)
        # needed for bookkeeping:
        self.rotation = [0, 0, 0]
        self.points = np.full((3, 12, size), np.nan, float)
        for i, crystal in enumerate(crystals):
            self.points[:, :, i] = crystal.points

    def move(self, xyz):
        # move the entire cluster
        # self.points[0] += xyz[0]
        # self.points[1] += xyz[1]
        # self.points[2] += xyz[2]
        for crystal in self.crystals:
            crystal.move(xyz)

    def _rotate(self, angles):
        for crystal in self.crystals:
            crystal._rotate(angles)

    def _rev_rotate(self, angles):
        # get back to the original rotation
        for crystal in self.crystals:
            crystal._rev_rotate(angles)

    def rotate_to(self, angles):
        # rotate the entire cluster

        # first get back to the original rotation
        self._rev_rotate(self.rotation)

        # now add the new rotation
        self._rotate(angles)

        # save the new rotation
        self.rotation = angles

    def center_of_mass(self):
        x = np.mean([ x.center[0] for x in self.crystals ])
        y = np.mean([ x.center[1] for x in self.crystals ])
        z = np.mean([ x.center[2] for x in self.crystals ])
        return [x, y, z]

    def plot(self):
        return geom.MultiLineString([ lines for crystal in self.crystals for lines in crystal.plot() ])
        # return geom.MultiLineString([ lines for n in range(len(self.crystals)) for lines in self._plot_one(n) ])

    def projectxy(self):
        polygons = [ crystal.projectxy() for crystal in self.crystals ]
        return shops.cascaded_union(polygons)

    def add_crystal_from_above(self, crystal, lodge=0):
        # drop a new crystal onto the cluster

        # first see which other crystals intersect with the new one in
        # the xy plane
        newpoly = crystal.projectxy()
        close_crystals = [ x for x in self.crystals if x.projectxy().intersects(newpoly) ]
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
        mindiffz = crystal.maxz - first_hit_lower_bound
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
        self.crystals.append(crystal)
        
        # fin.
        return True

    def reorient(self, method='schmitt', random_rotations=50):
        if method == 'schmitt':
            max_rot1 = None
            max_area = self.projectxy().area
            for i in range(random_rotations):
                [a, b, c] = [random.uniform(0, np.pi / 4), random.uniform(0, np.pi / 4), random.uniform(0, np.pi / 4)]
                # for mysterious reasons we are going to rotate this 3 times
                rot1 = [a, b, c]
                rot2 = [b * np.pi, c * np.pi, a * np.pi]
                rot3 = [c * np.pi / 2, a * np.pi / 2, b * np.pi / 2]
                self.rotate_to(rot1)
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
            if max_rot1 is not None:
                self.rotate_to(max_rot1)
                self._rotate(max_rot2)
                self._rotate(max_rot3)
                # if none of the new rotations was better we can leave
                # it at the original rotation
            self.rotation = [0, 0, 0]

    def _get_slice_x(self, p1, p2):
        # if the area goes up along the y-axis, we will have to
        # calculate things differently
        x1 = p1[0]
        x2 = p2[0]
        if (x1 / abs(x1)) != (x2 / abs(x2)):
            # first find height of triangle at x=zero
            y1 = p1[1]
            y2 = p2[1]
            m0 = (y1 - y2) / (x1 - x2)
            height = -m0 * x2 + y2
            # then integrate both sides
            m1 = -height / x1
            chunk1 = m1 * x1 ** 4 / 4 + height * x1 ** 3 / 3
            m2 = -height / x2
            chunk2 = -(m2 * x2 ** 4 / 4 + height * x2 ** 3 / 3)
            return abs(chunk1 + chunk2)
        # get the slopes from the origin to the points
        x1 = p1[0]
        x2 = p2[0]
        m1 = p1[1] / p1[0]
        m2 = p2[1] / p2[0]
        # get the integral of the main chunk
        main_chunk = (m2 - m1) * x2 ** 4 / 4
        x_range = x1 - x2
        if x_range != 0:
            x_height = (m2 - m1) * x2
            small_chunk = x_height *\
                          (x1 * (x1 ** 3 - x2 ** 3) / 3 - (x1 ** 4 - x2 ** 4) / 4) /\
                          x_range
        else:
            small_chunk = 0
        return abs(main_chunk + small_chunk)

    def _get_slice_xy(self, p1, p2):
        # figure out which point is in the middle
        points = [p1, p2, [0, 0]]
        xs = [ p[0] for p in points]
        xorder = np.argsort(xs)
        # the coordinates of the left, middle, and right points,
        # respectively
        [xl, yl] = points[xorder[0]]
        [xm, ym] = points[xorder[1]]
        [xr, yr] = points[xorder[2]]
        # slope of the line connecting the left and right points
        m3 = (yr - yl) / (xr - xl)
        b3 = -xl * m3 + yl
        # the y coordinate of the line connecting the left and right
        # points at the x position of the middle point
        m3_mid = yl + m3 * (xm - xl)
        # is the midpoint above or below that line?
        mid_above = ym > m3_mid
        # Writing this integral assuming the midpoint is above the
        # other points. If it's not then I have to flip the sign.
        if xm != xl:
            m1 = (ym - yl) / (xm - xl)
            b1 = -xl * m1 + yl
            leftside1 = (xm ** 4 - xl ** 4) * (m1 ** 2 - m3 ** 2) / 8
            leftside2 = (xm ** 3 - xl ** 3) * (b1 * m1 - b3 * m3) / 3
            leftside3 = (xm ** 2 - xl ** 2) * (b1 ** 2 - b3 ** 2) / 4
            leftside = leftside1 + leftside2 + leftside3
        else:
            leftside = 0
        if xm != xr:
            m2 = (yr - ym) / (xr - xm)
            b2 = -xr * m2 + yr
            rightside1 = (xr ** 4 - xm ** 4) * (m2 ** 2 - m3 ** 2) / 8
            rightside2 = (xr ** 3 - xm ** 3) * (b2 * m2 - b3 * m3) / 3
            rightside3 = (xr ** 2 - xm ** 2) * (b2 ** 2 - b3 ** 2) / 4
            rightside = rightside1 + rightside2 + rightside3
        else:
            rightside = 0

        if mid_above:
            return leftside + rightside
        else:
            return -(leftside + rightside)

    def _get_moments(self, poly):
        # get 'mass' moments for this polygon
        # using a variation of the shoelace algorithm
        points = list(poly.exterior.coords)
        npoints = len(points)
        xtotal = 0
        ytotal = 0
        xytotal = 0
        # the shapely polygons end with a copy of the first point, so
        # this completes the loop even though it only goes to
        # npoints-1
        for n in range(npoints - 1):
            p1 = list(points[n])
            p2 = list(points[n + 1])
            # are we going counterclockwise? if so make it positive
            # (shapely appears to store polygon vertices in
            # counterclockwise order, so this works out)
            counterclockwise = np.cross(p1, p2) > 0
            if counterclockwise:
                c = 1
            else:
                c = -1
            xtotal += c * self._get_slice_x(p1, p2)
            # instead of writing a y-specific function I'm just going
            # to transpose x and y whoa mind blown man
            ytotal += c * self._get_slice_x(list(reversed(p1)), list(reversed(p2)))
            xytotal -= c * self._get_slice_xy(p1, p2)
        # now loop through the holes as needed, subtract areas
        for linestring in list(poly.interiors):
            poly2 = geom.Polygon(linestring)
            # see if the polygon points are listed clockwise or
            # counterclockwise
            areas = []
            points = list(poly2.exterior.coords)
            for n in range(len(points) - 1):
                areas.append(np.cross(points[n], points[n+1]))
            clockwise = np.sum(areas) < 0
            if clockwise:
                c = -1
            else:
                c = 1
            # whoa recursion mind blown man
            minus_moments = self._get_moments(poly2)
            xtotal -= minus_moments[0]
            ytotal -= minus_moments[1]
            xytotal -= c * minus_moments[2]

        return [abs(xtotal), abs(ytotal), xytotal]

    def fit_ellipse(self):
        # Emulating this function, but for polygons in continuous
        # space rather than blobs in discrete space:
        # http://www.idlcoyote.com/ip_tips/fit_ellipse.html
        # Using a variation of the shoelace formula to calculate mass
        # moments for complicated polygons, integrating one triangular
        # slice at a time.
        poly = self.projectxy()
        xy_area = poly.area
        
        # temporarily center the cluster around the centroid
        centroid = poly.centroid
        poly = sha.translate(poly, -centroid.x, -centroid.y)

        # # see if the polygon points are listed clockwise or
        # # counterclockwise
        # areas = []
        # points = list(poly.exterior.coords)
        # for n in range(len(points) - 1):
        #     areas.append(np.cross(points[n], points[n+1]))
        # clockwise = np.sum(areas) < 0
        # if clockwise:
        #     c = 1
        # else:
        #     c = -1

        # if poly is actually a multipolygon, have to do a bit more
        # work
        if isinstance(poly, geom.multipolygon.MultiPolygon):
            x = 0
            y = 0
            xy = 0
            for poly2 in poly:
                # see if the polygon points are listed clockwise or
                # counterclockwise
                areas = []
                points = list(poly2.exterior.coords)
                for n in range(len(points) - 1):
                    areas.append(np.cross(points[n], points[n+1]))
                clockwise = np.sum(areas) < 0
                if clockwise:
                    c = -1
                else:
                    c = 1
                moments = self._get_moments(poly2)
                x += moments[0] / xy_area
                y += moments[1] / xy_area
                xy += c * (moments[2] / xy_area)
        else:
            # see if the polygon points are listed clockwise or
            # counterclockwise
            areas = []
            points = list(poly.exterior.coords)
            for n in range(len(points) - 1):
                areas.append(np.cross(points[n], points[n+1]))
            clockwise = np.sum(areas) < 0
            if clockwise:
                c = -1
            else:
                c = 1
            moments = self._get_moments(poly)
            x = moments[0] / xy_area
            y = moments[1] / xy_area
            xy = c * (moments[2] / xy_area)
            # the xy integral is calculated assuming the points are
            # ordered counterclockwise, if not have to multiply by minus 1

        # do magic
        m = np.matrix([[y, xy], [xy, x]])
        evals, evecs = np.linalg.eigh(m)
        semimajor = np.sqrt(evals[0]) * 2
        semiminor = np.sqrt(evals[1]) * 2
        major = semimajor * 2
        minor = semiminor * 2
        semiAxes = [semimajor, semiminor]
        axes = [major, minor]
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
        maxdim = max([params['width'], params['height']])# / 2
        ax.set_xlim([-maxdim + params['xy'][0], maxdim + params['xy'][0]])
        ax.set_ylim([-maxdim + params['xy'][1], maxdim + params['xy'][1]])
        ax.set_aspect('equal', 'datalim')

    def write_obj(self, filename):
        f = open(filename, 'w')

        faces = []
        for i, crystal in enumerate(self.crystals):
            # write the vertices
            for n in range(12):
                f.write('v ' + ' '.join(map(str, crystal.points[n])) + '\n')
            # write the hexagons
            nc = i * 12
            for n in range(2):
                coords = range(n * 6 + 1 + nc, (n + 1) * 6 + 1 + nc)
                faces.append('f ' + ' '.join(map(str, coords)))
            for n in range(5):
                coords = [n + 1 + nc, n + 2 + nc, n + 8 + nc, n + 7 + nc]
                faces.append('f ' + ' '.join(map(str, coords)))
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

        major_lengths = {}
        minor_lengths = {}
        for dim in ellipse.keys():
            major_lengths[dim] = max(ellipse[dim]['height'], ellipse[dim]['width'])
            minor_lengths[dim] = min(ellipse[dim]['height'], ellipse[dim]['width'])

        if method == 1:
            return max(major_lengths.values()) / max(minor_lengths.values())
        elif method == 'plate':
            return minor_lengths['y'] / major_lengths['z']
        elif method == 'column':
            return major_lengths['z'] / minor_lengths['y']
