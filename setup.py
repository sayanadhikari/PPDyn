#!/usr/bin/env python
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
with open('.version') as f:
    version = f.read().strip()

setuptools.setup(
    name="PPDyn", # Replace with your own username
    version=version,
    author="Sayan Adhikari",
    author_email="sayanadhikari207@gmail.com",
    description="A python package to simulate plasma particles using Molecular Dynamics Algorithm.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sayanadhikari/PPDyn",
    project_urls={
        "Bug Tracker": "https://github.com/sayanadhikari/PPDyn/issues",
    },
    packages=['PPDyn'],
    install_requires=['numpy', 'scipy', 'ini-parser', 'numba', 'h5py', 'matplotlib'],
    package_data={'': ['input.ini']},
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # package_dir={"": "PPDyn"},
    # packages=setuptools.find_packages(where="PPDyn"),
    python_requires=">=3.6",
)
