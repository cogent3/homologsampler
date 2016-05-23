#!/usr/bin/env python
from setuptools import setup
import sys, os, re, subprocess

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2014, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

# Check Python version, no point installing if unsupported version inplace
if sys.version_info < (2, 7):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError("Python-2.7 or greater is required, Python-%s used." % py_version)

short_description = "homologsampler"

# This ends up displayed by the installer
long_description = """homologsampler
Download homolog's from the Ensembl database
"""

setup(
    name="homologsampler",
    version=__version__,
    author="Gavin Huttley",
    author_email="gavin.huttley@anu.edu.au",
    description=short_description,
    long_description=long_description,
    platforms=["any"],
    license=["GPL"],
    keywords=["science", "bioinformatics", "genetics", "evolution"],
    classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Software Development :: Libraries :: Python Modules",
            "Operating System :: OS Independent",
            ],
    packages=['homologsampler'],
    install_requires=[
              'cogent==1.5.3-dev',
              'scitrack',
              'sqlalchemy',
              'PyMySQL',
          ],
    dependency_links=['https://github.com/GavinHuttley/pycogent/archive/master.zip#egg=cogent-1.5.3-dev'],
    entry_points={
            'console_scripts': ['one2one=homologsampler.__init__:main',
                            ],
        }
    url="https://bitbucket.org/gavin.huttley/homologsampler"
    )
