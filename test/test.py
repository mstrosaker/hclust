#!/usr/bin/env python

# Copyright (c) 2013 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import sys

try:
    import unittest2 as unittest
except ImportError:
    import unittest

sys.path.append("../hclust")
from hclust import DistanceMatrix, HClust

class TestDistanceMatrixConstructor(unittest.TestCase):

    def test_distance_matrix(self):
        with open('cytochromec.dist', 'rb') as f:
            dmat = DistanceMatrix(f)

        self.assertEqual(dmat.n_nodes, 7)

    def test_full_matrix(self):
        with open('bwa_cfg.dist', 'rb') as f:
            dmat = DistanceMatrix(f, full=True)

        self.assertEqual(dmat.n_nodes, 288)

class TestHierarchicalClustering(unittest.TestCase):

    def test_clustering(self):
        with open('cytochromec.dist', 'rb') as f:
            dmat = DistanceMatrix(f)

        hc = HClust(dmat)

        self.assertEqual(len(hc.clusters), 13)

    def test_cut(self):
        with open('cytochromec.dist', 'rb') as f:
            dmat = DistanceMatrix(f)

        hc = HClust(dmat)

        cluster_list = hc.cut(5)
        self.assertEqual(len(cluster_list), 5)

        cluster_list = hc.cut(10)
        self.assertEqual(len(cluster_list), 3)

        cluster_list = hc.cut(50)
        self.assertEqual(len(cluster_list), 1)

    def test_number_clusters(self):
        with open('bwa_cfg.dist', 'rb') as f:
            dmat = DistanceMatrix(f, full=True)

        hc = HClust(dmat)

        cluster_list = hc.n_clusters(5)
        self.assertEqual(len(cluster_list), 5)

        cluster_list = hc.n_clusters(10)
        self.assertEqual(len(cluster_list), 10)

        cluster_list = hc.n_clusters(1)
        self.assertEqual(len(cluster_list), 1)



if __name__ == '__main__':
    unittest.main()
