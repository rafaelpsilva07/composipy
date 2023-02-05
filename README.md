
# Composipy

## What it is

**composipy** is a python library to calculate composite plates using the classical laminate theory. This library is designed to be simple, userfriendly and helpfull.



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

[Composipy documentation](https://rafaelpsilva07.github.io/composipy/)

## Quick start

Create the material.

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

Define the laminate.

```python
>>> from composipy import LaminateProperty
>>> stacking = [-45, 45, 90, 0, 0, 0, 0, 90, 45, -45]
>>> laminate1 = LaminateProperty(stacking, mat_1)
>>> print(laminate1.ABD) # prints the ABD matrix as a np.ndarray
>>> print(laminate1.xiA) # prints lamination parameters of extension as a np.ndarray
>>> print(laminate1.xiD) # prints lamination parameters of bending as a np.ndarray
```

Create a plate structure.

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
>>> print(panel.buckling_analysis()) # solve the eigenvalue problem.
```
