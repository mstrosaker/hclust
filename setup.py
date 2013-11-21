# Copyright (c) 2013 Michael Strosaker
# MIT License
# http://opensource.org/licenses/MIT

import os, sys
from distutils.core import setup

try:
    with open('README.rst', 'rt') as readme:
        description = readme.read()
except IOError:
    description = ''

setup(
    # metadata
    name='hclust',
    description='Agglomerative Hierarchical Clustering',
    long_description=description,
    license='MIT License',
    version='0.10',
    author='Mike Strosaker',
    maintainer='Mike Strosaker',
    author_email='mstrosaker@gmail.com',
    url='https://github.com/mstrosaker/hclust',
    platforms='Cross Platform',
    classifiers = [
        'Programming Language :: Python :: 2',
        ],

    # All packages and sub-packages must be listed here
    packages=[
        'hclust',
        ],
)

