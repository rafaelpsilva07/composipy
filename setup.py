import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mimo_composipy", # Replace with your own username
    version="0.1.1",
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