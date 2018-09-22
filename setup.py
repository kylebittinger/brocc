import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='brocc',
    version='1.4.0',
    author='Serena Dollive, Kyle Bittinger',
    author_email='kylebittinger@gmail.com',
    maintainer='Kyle Bittinger',
    maintainer_email='kylebittinger@gmail.com',
    description='BLAST READ and OTU Consensus Classifier',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/kylebittinger/brocc',
    packages=['brocclib'],
    scripts=[
        'scripts/brocc.py', 'scripts/create_local_taxonomy_db.py',
        'scripts/compare_brocc_assignments.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
