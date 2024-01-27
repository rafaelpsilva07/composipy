import setuptools
import sys

sys.path.insert(0, ("./composipy"))

from version import __version__


with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    requires = f.readlines()
    requires = [r.replace('\n', '') for r in requires]


setuptools.setup(
    name="composipy",
    version=__version__,
    install_requires=requires,
    author="Rafael Pereira",
    author_email="rafaelpsilva07@gmail.com",
    description="This package intends to perform composite material calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rafaelpsilva07/composipy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)