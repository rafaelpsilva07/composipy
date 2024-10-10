[![Downloads](https://static.pepy.tech/badge/composipy)](https://pepy.tech/project/composipy)
[![Deploy](https://github.com/rafaelpsilva07/composipy/actions/workflows/python-publish.yml/badge.svg)](https://github.com/rafaelpsilva07/composipy/actions/workflows/python-publish.yml)
[![Test](https://github.com/rafaelpsilva07/composipy/actions/workflows/pytest_test.yml/badge.svg)](https://github.com/rafaelpsilva07/composipy/actions/workflows/pytest_test.yml)
![PyPI](https://img.shields.io/pypi/v/composipy)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/rafael-pereira-da-silva07/)
[![Youtube](https://img.shields.io/badge/YouTube-FF0000?style=for-the-badge&logo=youtube&logoColor=white)](https://www.youtube.com/watch?v=W4foqYR8IL0&t=622s&ab_channel=EngenhariaRefinada)
[![DOI](https://zenodo.org/badge/332543985.svg)](https://zenodo.org/badge/latestdoi/332543985)


# Composipy

## What it is

Composipy is a Python-based library designed to address the challenges of composite plate analysis and optimization in the aerospace industry. The library offers tools for laminate stress-strain, plate buckling analysis and lamination parameter optimization, helping engineers during the design process. Utilizing object-oriented programming and native Python structures, Composipy ensures a intuitive workflow and easy integration into existing engineering practices, such as defining material, defining properties, so on.

The library is built using the most powerful Python numerical libraries, so users can think of using Composipy in the same way as Pandas, NumPy, and SciPy. In fact, Composipy is built using their objects and structures.


<img src="https://github.com/rafaelpsilva07/composipy/blob/main/doc/images/composipy_features.PNG" width="700">


It is especially useful for leveraging the data-driven culture in companies and can be used to build response surfaces, generate samples, and much more. Refer to the [Composipy Examples](https://rafaelpsilva07.github.io/composipy/notebooks/index.html)


## Main features
- ABD Matrix
- Laminate Stress and Strain
- Lamination Parameters
- Plate Buckling
- Plate Optimization


## How to install

### Directly From PYPI

```shell
pip install composipy
```

### Built from source

```shell
python setup.py install
```

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


### Create a Laminate Strength Analysis.
```python
from composipy import LaminateStrength

laminate_strength = LaminateStrength(laminate1, Nxx=100, Mxx=10)
laminate_strength.epsilon0() #strains at the midplane
laminate_strength.calculate_strain() #strain ply by ply
laminate_strength.calculate_strain() #stress ply by ply
```
See [LaminateStrength](https://rafaelpsilva07.github.io/composipy/reference/classes.html#laminatestrength) for reference.

Also, check the [Stress Strain Calculation of a Laminate](https://rafaelpsilva07.github.io/composipy/notebooks/Stress_strain_of_laminate.html) to see a complete example.

### Create a Plate Structure for Buckling Analysis

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


#### Calculate Buckling
```python
print(panel.buckling_analysis()) # solve the eigenvalue problem.
```

#### Plot Buckling shape mode
```python
print(panel.plot_eigenvalue())
```

Composipy is able to perform buckling calculation considering different boundary conditions and in-plane load applications. Check the [Critical Buckling Examples With Varying Boundary Conditions](https://rafaelpsilva07.github.io/composipy/notebooks/Critical_buckling_varying_BCs.html) to see a complete example.


If you are interesting in plate optimization, you may want to check the [Optimization of a Plate Subjected to Buckling Loads](https://rafaelpsilva07.github.io/composipy/notebooks/Optimization_buckling.html)


## Theoretical References

The implementation of composipy is based on the following reference:

- SILVA, Rafael Pereira da. Composite Plate optimization combining semi-analytical model, Lamination Parameters and a Gradient-Based Optimizer. 2023. 82f. Dissertation of Master of Science – Instituto Tecnológico de Aeronáutica, São José dos Campos.
- rafaelpsilva07. (2024). rafaelpsilva07/rafaelmscdissertation: v1.0.0 (1.0.0). Zenodo. https://doi.org/10.5281/zenodo.10546621
- Application repository: https://github.com/rafaelpsilva07/rafaelmscdissertation