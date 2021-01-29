# Overview

**mimo_composipy** is a python library to calculate composite plates using the classical laminate theory. This library is designed to be simple, userfriendly and helpfull.
Now, you can perform a composite plate buckling calculation using 3 lines of python code.


This library is a creation of the Techmimo project.


## Download using PYPI

**pip install mimo-composipy**

Access the PYPI project:

https://pypi.org/project/mimo-composipy/





# Realeases

## v 0.0.1 (2021/02)

Basic release. It contains:

- **Ply** instances to calculate a lamina macromechanical behavior
- **Laminate** instances to perform laminate calculation
- **critical_buckling** function that calculates the critical buckling load of a composite plate


# First steps

## Application example:
In this example we're going to perform a buckling calculation from the scretch. Consider the following composite plate:

### Ply mechanical properties
$$ E_1 = 129500 MPa $$
$$ E_2 = 9370 MPa $$
$$ v12 = 0.38 MPa $$
$$ G12 = 5240 MPa $$
$$ ply_{thickness} = 0.2 mm $$

### Laminate stacking sequence
$$ 90°/90°/0°/90°/90°/0°/90°/90° $$

### Plate dimensions
$$ plate_{width} = 360mm$$
$$ plate_{length} = 360mm$$

## Coding solution
```python
>>> #Importing the objects and functions
>>> from mimo_composipy import Laminate, Ply, critical_buckling

>>> #Creating the variables
>>> E1 = 129500
>>> E2 = 9370
>>> v12 = 0.38
>>> G12 = 5240
>>> t = 0.2

>>> #Creating the Ply instance
>>> ply_1 = Ply(E1,E2,v12,G12,t)

>>> #Creating the layup argument
>>> #The argument is a list filled with tuples containing the information
>>> layup_1 = [(90,ply_1),(90,ply_1),(0,ply_1),(90,ply_1),(90,ply_1),(0,ply_1),(90,ply_1),(90,ply_1)]
>>>
>>> #Creating the Laminate instance
>>> laminate_1 = Laminate(layup_1)

>>> #Accessing properties

>>> #Lamina stiffness matrix
>>> print(ply_1.Q_0)
>>> print('\n')

>>> #Laminate [A][B][D] matrices
>>> print(laminate_1.A)
>>> print('\n')
>>> print(laminate_1.B)
>>> print('\n')
>>> print(laminate_1.D)

>>> #Calculating the critical buckling load and ploting the shape
>>> critical_buckling( plate_width, plate_lenght, laminate_1.D,method = 'polynomial',n=3)
>>> # method may be 'polinomial' or 'sin_series'
>>> # n is the degree of refinement
```
**The outputs are:**

![instances_printing](https://github.com/rafaelpsilva07/mimo_composipy/blob/main/images/instance_properties.PNG)

![buckling_results](https://github.com/rafaelpsilva07/mimo_composipy/blob/main/images/buckling_output.PNG)
