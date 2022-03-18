import setuptools
import sys
import os

#sys.path.insert(0, os.path.abspath("./"))

#from composipy._version import VERSION

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="composipy", # Replace with your own username
    version="0.2.0",
    install_requires=["numpy>=1.0", "scipy>=1.0"],
    author="Rafael Pereira",
    author_email="rafaelpsilva07@gmail.com",
    description="This package intends to perform composite material calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rafaelpsilva07/mimo_composipy.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)