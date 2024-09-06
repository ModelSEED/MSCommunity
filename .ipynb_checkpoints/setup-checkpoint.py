# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

# with open("README.rst") as f:
#     readme = f.read()

# with open("LICENSE") as f:
#     license = f.read()

setup(
    name="MSCommunity",
    version="0.0.1",
    description="Python package for building and analyzing microbial communities using ModelSEED",
    # long_description_content_type="text/x-rst",
    # long_description=readme,
    author="Andrew Freiburger",
    author_email="afreiburger@anl.gov",
    url="https://github.com/ModelSEED/MSCommunity",
    # license=license,
    packages=find_packages(),
    # package_data={
    #     "modelseedpy": ["config.cfg", "community/*.html", "core/*.html", "data/*", "data/categories/*", "data/templates/*"],
    # },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Natural Language :: English",
    ],
    include_package_data =True,
    install_requires=[
        "networkx >= 2.4",
        "cobra >= 0.28.0",
        "scikit-learn >= 1.2.0",
        "scipy >= 1.5.4",
        "chemicals >= 1.0.13",
        "chemw >= 0.3.2",
        "matplotlib >= 3.0.0",
        "icecream",
        "deepdiff",
        "openpyxl",
        "jinja2",
        "multiprocess",
        "h5py",
        "graphviz"
    ],
    tests_require=[
        "pytest",
    ],
    project_urls={
        # "Documentation": "https://modelseedpy.readthedocs.io/en/latest/",
        "Issues": "https://github.com/ModelSEED/MSCommunity/issues",
    },
)
