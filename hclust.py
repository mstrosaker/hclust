#!/usr/bin/env python

# Copyright (c) 2013 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import sys
from copy import deepcopy


class Node(object):
    """
    Internal representation of a node for hierarchical clustering.  A node
    may be a leaf, or a cluster of multiple objects (leaves or other clusters).
    Instances of this class are hashable.

    The only use for this class outside of this module is when accessing
    HClust.clusters to obtain a list of all the clusters.
    """
    def __init__(self, id, depth, parent=None, children=None):
        self.id = id
        self.parent = parent
        self.children = children
        self.depth = depth

    def __repr__(self):
        return '(%s %f)' % (self.id, self.depth)

    def __key(self):
        return (self.id, self.depth)

    def __eq__(x, y):
        return type(x) == type(y) and x.__key() == y.__key()

    def __hash__(self):
        return hash(self.__key())


def bottom_triangle(stream):
    """
    This is a generator that returns the bottom triangle of a full matrix
    of distance values.  The matrix must be symmetric across the diagonal.
    For example, it converts:

        'title  id1     id2     id3     id4'
        'id1    0.0     2.0     1.5     2.5'
        'id2    2.0     0.0     3.0     1.0'
        'id3    1.5     3.0     0.0     4.5'
        'id4    2.5     1.0     4.5     0.0'

    to:

        '       id1     id2     id3'
        'id2    2.0'
        'id3    1.5     3.0'
    """
    line_num = 0
    for line in stream:
        line_num += 1
        entries = line.rstrip().split('\t')
        if line_num == 1:
            yield '\t' + '\t'.join(entries[1:-1])
        elif line_num > 2:
            yield '\t'.join(entries[0:line_num-1])


class DistanceMatrix:
    def __init__(self, stream, full=False):
        """
        Create a new distance matrix object based on the filelike 'stream'.
        The strings from the 'stream' parameter can be in one of two formats;
        the format is specified with the 'full' parameter (defaults
        to the first format):

        DISTANCE MATRIX:

        '	id1	id2	id3'
        'id2	2.0'
        'id3	1.5	3.0'
        'id4	2.5	1.0	4.5'

        There are a total of 4 observations in this distance matrix.
        Note that there is no row corresponding to the first observation,
        and no column corresponding to the last.  The entries in each line
        must be tab-delimited.

        FULL MATRIX:

        'title	id1	id2	id3	id4'
        'id1	0.0	2.0	1.5	2.5'
        'id2	2.0	0.0	3.0	1.0'
        'id3	1.5	3.0	0.0	4.5'
        'id4	2.5	1.0	4.5	0.0'
	
        There are a total of 4 observations.  The matrix must be symmetric
        across the diagonal.  The entries in each line must be tab-delimited.

        Note that the first distance matrix is the bottom triangle of the
        full matrix.
        """
        if full:
            f = bottom_triangle(stream)
        else:
            f = stream
        self.obs = []
        self.distances = {}
        for row in f:
            row = row.rstrip()
            if row[0] == '\t':
                # header row
                self.obs = row.split('\t')[1:]
            else:
                entries = row.split()
                self.distances[entries[0]] = {}
                i = 0
                for entry in entries[1:]:
                    self.distances[entries[0]][self.obs[i]] = float(entry)
                    i += 1
                if entries[0] not in self.obs:
                    self.obs.append(entries[0])

    def closest(self):
        """
        Return a tuple containing the two entries in the distance matrix
        that are closest together, or None if there are less than two in
        the matrix.
        """
        if len(self.obs) < 2:
            return None

        current_closest = None
        current_a = None
        current_b = None

        for a, bs in self.distances.iteritems():
            for b, dist in bs.iteritems():
                if current_closest is None or \
                            dist < current_closest:
                    current_closest = dist
                    current_a = a
                    current_b = b

        return (current_a, current_b)

    def distance(self, a, b):
        """
        Return the distance between two nodes, or None if one or both of
        nodes are not in the distance matrix.
        """
        if a in self.distances:
            if b in self.distances[a]:
                return self.distances[a][b]
        if b in self.distances:
            if a in self.distances[b]:
                return self.distances[b][a]
        return None

    def to_nodes(self):
        """
        Convert the members of the matrix to Node objects.
        """
        new_distances = {}
        for a in self.distances:
            new_a = Node((a,), 0)
            new_distances[new_a] = {}
            for b in self.distances[a]:
                new_b = Node((b,), 0)
                new_distances[new_a][new_b] = self.distances[a][b]
        self.distances = new_distances

        new_obs = []
        for o in self.obs:
            new_obs.append(Node((o,), 0))
        self.obs = new_obs

    def _get_n_nodes(self):
        """
        Return the number of observations in the distance matrix.
        """
        return len(self.distances) + 1

    n_nodes = property(_get_n_nodes,
                       doc='the number of nodes in the distance matrix')


class HClust:
    def __init__(self, dist_matrix, linkage_criterion='average'):
        """
        Perform agglomerative hierarchical clustering on the DistanceMatrix
        object specified as dist_matrix.  The possible values for
        linkage_criterion are:
          'average': the mean distance between the elements of each
                     cluster, also known as average linkage clustering
                     and used in UPGMA
          'max': the maximum distance between the elements of each cluster,
                 also known as complete-linkage clustering
          'min': the minimum distance between the elements of each cluster,
                 also known as single-linkage clustering
        """
        # this distance matrix will not be modified; it will maintain
        # the original user-specified distances
        self.dist_matrix = dist_matrix

        self.linkage_criterion = linkage_criterion

        # create a working distance matrix whose rows/columns will
        # be merged during clustering
        self.working_matrix = deepcopy(dist_matrix)
        self.working_matrix.to_nodes()

        self.clusters = deepcopy(self.working_matrix.obs)

        while True:
            nodes = self.working_matrix.closest()
            if nodes is None:
                break
            self._merge(nodes[0], nodes[1])

    def _linkage(self, a, b):
        """
        Calculate the new distance between two clusters based on the
        linkage criteria specified when the object was created.
        """
        if self.linkage_criterion == 'average':
            sum = 0.0
            n_obs = 0
            for obs_a in a.id:
                for obs_b in b.id:
                    sum += self.dist_matrix.distance(obs_a, obs_b)
                    n_obs += 1
            return sum / n_obs

        elif self.linkage_criterion == 'max':
            m = 0
            for obs_a in a.id:
                for obs_b in b.id:
                    d = self.dist_matrix.distance(obs_a, obs_b)
                    if d > m:
                        m = d
            return m

        elif self.linkage_criterion == 'min':
            m = sys.maxint
            for obs_a in a.id:
                for obs_b in b.id:
                    d = self.dist_matrix.distance(obs_a, obs_b)
                    if d < m:
                        m = d
            return m

        else:
            raise ValueError('invalid linkage criterion specified: %s' % \
                             self.linkage_criterion)

    def _cluster_parent(self, child, parent):
        """
        Find a node in the list of all clusters and assign a parent to it.
        """
        for n in self.clusters:
            if n == child:
                n.parent = parent

    def _merge(self, a, b):
        """
        Merge two nodes into a cluster.
        """
        c = Node(a.id + b.id, self.working_matrix.distance(a, b) * 0.5,
                 children=(a, b))
        self._cluster_parent(a, c)
        self._cluster_parent(b, c)
        self.working_matrix.obs.remove(a)
        self.working_matrix.obs.remove(b)
        self.clusters.append(c)

        # remove the nodes that were replaced by the new cluster
        if a in self.working_matrix.distances:
            del self.working_matrix.distances[a]
        if b in self.working_matrix.distances:
            del self.working_matrix.distances[b]

        # from remaining nodes, remove the distances to the removed nodes
        for x, ys in self.working_matrix.distances.iteritems():
            if a in ys:
                del ys[a]
            if b in ys:
                del ys[b]

        # for remaining nodes, add distances to the new node
        added = []
        for x in self.working_matrix.distances:
            self.working_matrix.distances[x][c] = self._linkage(x, c)
            added.append(x)

        self.working_matrix.distances[c] = {}

        # for the new node, calculate distances to each existing node
        for x in self.working_matrix.obs:
            if x not in added:
                self.working_matrix.distances[c][x] = self._linkage(x, c)

        if len(self.working_matrix.distances[c]) == 0:
            del self.working_matrix.distances[c]

        self.working_matrix.obs.append(c)

    def _ids(self, l):
        """
        Change a list of Node objects to a list of identifying strings.
        """
        ret = []
        for node in l:
            ret.append(node.id)
        return ret

    def _leaves(self):
        """
        Return a list of the Node objects corresponding to the leaves
        of the dendrogram.
        """
        l = []
        for clust in self.clusters:
            if clust.children is None:
                l.append(clust)
        return l

    def _get_leaves(self):
        """
        Return a list of IDs corresponding to the leaves of the dendrogram.
        """
        return self._ids(self._leaves())

    def _trunk(self):
        """
        Return the Node object corresponding to the trunk (root) of the
        dendrogram.
        """
        return self.clusters[-1]

    def _get_trunk(self):
        """
        Return the ID corresponding to the trunk (root) of the dendrogram.
        """
        return self._trunk.id

    trunk = property(_get_trunk, doc='the trunk (root) of the dendrogram')
    leaves = property(_get_leaves, doc='the leaves of the dendrogram')

    def cut(self, n):
        """
        Return the tuple of clusters obtained when the dendrogram is cut
        at the specified depth n.

        Always returns a tuple, even if there is only one member.
        """
        if n < 0:
            return None

        lvs = self._leaves()
        ret = []

        for l in lvs:
            node = l
            while node.parent is not None and node.parent.depth <= n:
                node = node.parent
            if node not in ret:
                ret.append(node)

        return self._ids(ret)

    def n_clusters(self, n):
        """
        Return the tuple of clusters obtained when the clustering is stopped
        at the specified number of clusters n.

        Always returns a tuple, even if there is only one member.
        """
        if n <= 0:
            return None

        ret = self._leaves()

        while len(ret) > n:
            # find the parent with the smallest depth
            lowest = None
            for node in ret:
                if node.parent.children[0] in ret and \
                                node.parent.children[1] in ret:
                    if lowest is None:
                        lowest = node.parent
                    elif node.parent.depth < lowest.depth:
                        lowest = node.parent
            ret.remove(lowest.children[0])
            ret.remove(lowest.children[1])
            ret.append(lowest)

        return self._ids(ret)

