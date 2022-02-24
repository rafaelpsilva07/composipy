Getting started
===============

Installation
------------

The composipy can be installed via pip from PyPI.

>>> pip install composipy


Intro to composipy
------------------
Composipy was developed to perform composites mechanics calculation analytically.

.. note::
	Current version of composipy contains Ply and Laminate classes.
	The Rayleigh Ritz method for buckling calculation has been removed for simplification purposes.
	But buckling functions can be found in previous releases (see v0.1.3 and previous).

Examples
--------
This example presents the creation of a Laminate from the scratch.

Importing Ply and Laminate classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> from composipy import Ply
>>> from composipy import Laminate

Defining Ply instances
~~~~~~~~~~~~~~~~~~~~~~
The following lines presents how an instance of Ply is created.

>>> ply_1 = Ply(129500, 9370, 0.38, 5240, 0.2, name='Tape Carbon Fiber')
>>> ply_2 = Ply(129500, 129500, 0.38, 5240, 0.2, name='Fabric Carbon Fiber')

Builting a Laminate
~~~~~~~~~~~~~~~~~~~
>>> layup_1 = [
    (90, ply_1),
    (90, ply_2),
    (0, ply_1),
    (45, ply_2),
    (-45, ply_2),
    (0, ply_1),
    (90, ply_2),
    (90, ply_1)
    ]
>>> laminate_1 = Laminate(layup_1)

.. note::
    You made it! Now lets explore the options behind your composite material!

A nice representation of your material!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> print(ply_1)
Out:
Ply(Tape Carbon Fiber, E1 = 129500, E2 = 9370,
    v12 = 0.38, G12 = 5240, thickness = 0.2)
>>> print(laminate_1)
Out:
Tape Carbon Fiber   90                  |||||
Fabric Carbon Fiber 90                  |||||
Tape Carbon Fiber   0                   =====
Fabric Carbon Fiber 45                  /////
Fabric Carbon Fiber -45                 /////
Tape Carbon Fiber   0                   =====
Fabric Carbon Fiber 90                  |||||
Tape Carbon Fiber   90                  |||||

Get compliance Matrix Q of a ply
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> ply_1.Q_0
Out:
array([[130867.31382151,   3598.19426713,      0.        ],
       [  3598.19426713,   9468.93228191,      0.        ],
       [     0.        ,      0.        ,   5240.        ]])

Get A, B, D matrices of the laminate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> laminate_1.A
Out:
array([[ 1.60547001e+05,  6.55628266e+04, -1.92885098e-12],
       [ 6.55628266e+04,  1.60547001e+05,  4.90225377e-12],
       [-1.92885098e-12,  4.90225377e-12,  2.50561159e+04]])

.. note::
    Try to reproduce the B and D matrices.

Show the complete result for ABD matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> laminate_1.show_ABD()
[A]:
[[ 1.60547001e+05  6.55628266e+04 -1.92885098e-12]
 [ 6.55628266e+04  1.60547001e+05  4.90225377e-12]
 [-1.92885098e-12  4.90225377e-12  2.50561159e+04]]
[B]:
[[ 2.91038305e-11  1.21076482e-11 -4.73316543e-28]
 [ 1.21076482e-11  2.54658516e-11  1.00974196e-27]
 [-4.73316543e-28  1.00974196e-27  4.09272616e-12]]
[D]:
[[ 2.26765700e+04  7.20162516e+03 -4.61547603e-13]
 [ 7.20162516e+03  4.21003111e+04  1.92842631e-12]
 [-4.61547603e-13  1.92842631e-12  2.01088155e+03]]