[![Downloads](https://static.pepy.tech/badge/composipy)](https://pepy.tech/project/composipy)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/rafaelpsilva07/composipy/python-publish.yml)
![PyPI](https://img.shields.io/pypi/v/composipy)
[![DOI](https://zenodo.org/badge/332543985.svg)](https://zenodo.org/badge/latestdoi/332543985)


# Composipy

## What it is

**composipy** is a python library to calculate composite plates using the classical laminate theory. This library is designed to be simple, userfriendly and helpfull.

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

- [Composipy documentation](https://rafaelpsilva07.github.io/composipy/#contents)

- [Examples](https://rafaelpsilva07.github.io/composipy/notebooks/index.html)

- [Code References](https://rafaelpsilva07.github.io/composipy/reference/index.html)


## Quick start

### Create the Material Properties.

```python
>>> from composipy import OrthotropicMaterial
>>> 
>>> E1 = 60800
>>> E2 = 58250
>>> v12 = 0.07
>>> G12 = 4550
>>> t = 0.21
>>>
>>> mat_1 = OrthotropicMaterial(E1, E2, v12, G12, t)
```

See [OrthotropicMaterial](https://rafaelpsilva07.github.io/composipy/reference/classes.html) for reference.


### Define the Laminate.

```python
>>> from composipy import LaminateProperty
>>> stacking = [-45, 45, 90, 0, 0, 0, 0, 90, 45, -45]
>>> laminate1 = LaminateProperty(stacking, mat_1)
```

See [LaminateProperty](https://rafaelpsilva07.github.io/composipy/reference/classes.html#laminateproperty) for reference.

### Calculate Stiffnnes Matrix and Lamination Parameters

```python
>>> print(laminate1.ABD) # prints the ABD matrix as a np.ndarray
>>> print(laminate1.xiA) # prints lamination parameters of extension as a np.ndarray
>>> print(laminate1.xiD) # prints lamination parameters of bending as a np.ndarray
```

### Create a Plate Structure.

```python
>>> from composipy import PlateStructure
>>> 
>>> constraints = {    
---     'x0' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'],
---     'xa' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'],
---     'y0' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ'],
---     'yb' : ['TX', 'TY', 'TZ', 'RX', 'RY', 'RZ']
--- }
>>> panel = PlateStructure(laminate1, 360, 360, m=10, n=10, Nxx=-1, constraints=constraints)
```

See [PlateStructure](https://rafaelpsilva07.github.io/composipy/reference/classes.html#platestructure) for reference.


### Calculate Buckling
```python
>>> print(panel.buckling_analysis()) # solve the eigenvalue problem.
```

### Plot Buckling shape mode
```python
>>> print(panel.plot_eigenvalue())
```



## Theoretical References

The implementation of composipy is based on the following references

- [Mechanics Of Composite Materials by Robert M. Jones](https://www.routledge.com/Mechanics-Of-Composite-Materials/Jones/p/book/9781560327127)
- [Castro S.G.P., Donadon M.V. Assembly of semi-analytical models to address linear buckling and vibration of stiffened composite panels with debonding defect. Compos. Struct., 160 (2017), pp. 232-247,](https://www.sciencedirect.com/science/article/abs/pii/S026382231631008X)
