#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''
Software tools for analysis data from the RAFT (rapid amplification of forum termini) protocol
to compute sites of DNA double stranded breakpoint hotspots"
'''


setup(
    name='raft_hotspots_se',
    version='0.1.0.0',
    author='Daniel Park',
    author_email='djp@unimelb.edu.au',
    packages=['raft_hotspots_se'],
    package_dir={'raft_hotspots_se': 'raft_hotspots_se'},
    entry_points={
        'console_scripts': [
            'raft_fastq_2sites_parse = raft_hotspots_se.raft_fastq_2sites_parse:main',
            'raft_bed_2sites_parse = raft_hotspots_se.raft_bed_2sites_parse:main']
    },
    url='https://github.com/bjpop/raft_hotspots_se',
    license='LICENSE',
    description=('Compute DNA double stranded breakpoint hotspots from RAFT sequencing data'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["biopython==1.66"],
)
