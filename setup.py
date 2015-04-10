"""Distutils file for metagenomix."""

from setuptools import setup, find_packages
import inspect
import shutil
import os

from metagenomix import __version__


requires = [
    'sqlalchemy',
    'lxml',
    'biopython',
    'numpy',
    'pysam',
    'matplotlib'
]

setup(
    name = 'metagenomix',
    version = __version__,
    description = 'Metagenomic analysis pipeline',
    url = 'github/metagenomix',
    author = 'Ana Bulovic',
    author_email = 'bulovic.ana@gmail.com',
    license = 'MIT',
    long_description = open('README.txt').read(),
    packages = find_packages(),
    scripts = ['bin/store-input-to-db', 'bin/store-blast-to-db',
               'bin/analyze-alignment', 'bin/prot-freq-db', 'bin/parse-mgb',
               'bin/greedy-OTU-assign', 'bin/db-analysis', 'bin/extract-subtaxa',
               'bin/lca-read-assign', 'bin/dmp2taxtree', 'bin/generate-config',
               'bin/host-greedy-otu', 'bin/metagenomix-ms-eval', 'bin/create-msim-profile',
               'bin/extract-reads', 'bin/simple-profile', 'bin/detailed-profile',
               'bin/generate-consensus'],
    package_data = {'metagenomix': ['taxid2namerank', 'ncbi_tax_tree'],
                    'metagenomix.visualization': ['*.html', '*.css', '*.js']},
    data_files = [('', ['README.txt'])],
    install_requires=requires,
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: Freeware',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    )
