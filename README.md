[![Downloads](https://static.pepy.tech/badge/composipy)](https://pepy.tech/project/composipy)
[![Deploy](https://github.com/rafaelpsilva07/composipy/actions/workflows/python-publish.yml/badge.svg)](https://github.com/rafaelpsilva07/composipy/actions/workflows/python-publish.yml)
[![Test](https://github.com/rafaelpsilva07/composipy/actions/workflows/pytest_test.yml/badge.svg)](https://github.com/rafaelpsilva07/composipy/actions/workflows/pytest_test.yml)
![PyPI](https://img.shields.io/pypi/v/composipy)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/rafael-pereira-da-silva07/)
[![Youtube](https://img.shields.io/badge/YouTube-FF0000?style=for-the-badge&logo=youtube&logoColor=white)](https://www.youtube.com/watch?v=W4foqYR8IL0&t=622s&ab_channel=EngenhariaRefinada)
[![DOI](https://zenodo.org/badge/332543985.svg)](https://zenodo.org/badge/latestdoi/332543985)


# Composipy

## What it is

Composipy is a Python-based library designed to address the challenges of composite plate analysis and optimization in the aerospace industry, where weight reduction is crucial for efficient and profitable aircraft design. The library offers tools for plate buckling analysis and lamination parameter optimization, empowering engineers to streamline the design process and enhance structural performance. Utilizing object-oriented programming and native Python structures, Composipy ensures a seamless workflow and easy integration into existing engineering practices. Through continuous integration and delivery practices, Composipy maintains reliability and efficiency. 

Composipy is able to perform buckling calculation considering different boundary conditions and in-plane load applications. See [Examples](https://rafaelpsilva07.github.io/composipy/notebooks/Examples_BCs.html).


<img src="https://github.com/rafaelpsilva07/composipy/blob/main/doc/images/load_bcs_examples.PNG" width="700">



## How to install

### Directly From PYPI

```shell
pip install composipy
```

### Built from source

```shell
python setup.py install
```

## Documentation

- [Examples](https://rafaelpsilva07.github.io/composipy/notebooks/index.html)

- [Code References](https://rafaelpsilva07.github.io/composipy/reference/index.html)

- [Optitmization functions](https://rafaelpsilva07.github.io/composipy/reference/optimization_functions.html)

## Quick start

### Create the Material Properties.

```python
from composipy import OrthotropicMaterial
 
E1 = 60800
E2 = 58250
v12 = 0.07
G12 = 4550
t = 0.21

mat_1 = OrthotropicMaterial(E1, E2, v12, G12, t)
```

See [OrthotropicMaterial](https://rafaelpsilva07.github.io/composipy/reference/classes.html) for reference.


### Define the Laminate.

```python
from composipy import LaminateProperty
stacking = [-45, 45, 90, 0, 0, 0, 0, 90, 45, -45]
laminate1 = LaminateProperty(stacking, mat_1)
```

See [LaminateProperty](https://rafaelpsilva07.github.io/composipy/reference/classes.html#laminateproperty) for reference.

### Calculate Stiffnnes Matrix and Lamination Parameters

```python
print(laminate1.ABD) # prints the ABD matrix as a np.ndarray
print(laminate1.xiA) # prints lamination parameters of extension as a np.ndarray
print(laminate1.xiD) # prints lamination parameters of bending as a np.ndarray
```

### Create a Plate Structure.

```python
from composipy import PlateStructure
 
constraints = {    
     'x0' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'],
     'xa' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'],
     'y0' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'],
     'yb' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
}
panel = PlateStructure(laminate1, 360, 360, m=10, n=10, Nxx=-1, constraints=constraints)
```

See [PlateStructure](https://rafaelpsilva07.github.io/composipy/reference/classes.html#platestructure) for reference.


### Calculate Buckling
```python
print(panel.buckling_analysis()) # solve the eigenvalue problem.
```

### Plot Buckling shape mode
```python
print(panel.plot_eigenvalue())
```


## Theoretical References

The implementation of composipy is based on the following reference:

- SILVA, Rafael Pereira da. Composite Plate optimization combining semi-analytical model, Lamination Parameters and a Gradient-Based Optimizer. 2023. 82f. Dissertation of Master of Science – Instituto Tecnológico de Aeronáutica, São José dos Campos.
- rafaelpsilva07. (2024). rafaelpsilva07/rafaelmscdissertation: v1.0.0 (1.0.0). Zenodo. https://doi.org/10.5281/zenodo.10546621
- Application repository: https://github.com/rafaelpsilva07/rafaelmscdissertation