import unittest
from unittest import TestCase
import copy as cp
import ipas.crystals as crys

class TestCrystals(TestCase):
    def test_crystal_moves(self):
        crystal = crys.IceCrystal(1, 1)
        x0 = cp.copy(crystal.points['x'])
        y0 = cp.copy(crystal.points['y'])
        z0 = cp.copy(crystal.points['z'])
        crystal.move([1, 1, 1])
        self.assertTrue((crystal.points['x'] == x0 + 1).all())
        self.assertTrue((crystal.points['y'] == y0 + 1).all())
        self.assertTrue((crystal.points['z'] == z0 + 1).all())

if __name__ == '__main__':
    unittest.main()
