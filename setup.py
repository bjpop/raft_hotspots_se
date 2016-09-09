#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''The program reads one or more input FASTA files. 
For each file it computes a variety of statistics, and then
prints a summary of the statistics as output.
        
The goal is to provide a solid foundation for new bioinformatics command line tools,
and is an ideal starting place for new projects.'''


setup(
    name='raft_hotspots_se-py',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['raft_hotspots_se'],
    package_dir={'raft_hotspots_se': 'raft_hotspots_se'},
    entry_points={
        'console_scripts': ['raft_hotspots_se-py = raft_hotspots_se.raft_hotspots_se:main']
    },
    url='https://github.com/bjpop/raft_hotspots_se',
    license='LICENSE',
    description=('A prototypical bioinformatics command line tool'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["biopython==1.66"],
)
