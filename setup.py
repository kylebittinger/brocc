#!/usr/bin/env python
from distutils.core import setup

setup(
    name='BLAST READ and OTU Consensus Classifier',
    version='1.3.0',
    description='BLAST READ and OTU Consensus Classifier',
    author='Serena Dollive, Kyle Bittinger',
    author_email='kbit@mail.med.upenn.edu',
    maintainer='Kyle Bittinger',
    maintainer_email='kbit@mail.med.upenn.edu',
    url='',
    packages=['brocclib'],
    scripts=[
        'scripts/brocc.py', 'scripts/create_local_taxonomy_db.py',
        'scripts/compare_brocc_assignments.py'],
)
