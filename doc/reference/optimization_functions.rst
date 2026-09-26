Composipy Optimization Functions
================================
This page contains functions used to optimize the composite plate.
Optimization algorithms are gradient based and they use Scipy as optimization engine.

Since optimization is continuous they result in Lamination Parameters as output. 
Lamination Parameters are then supposed to be converted into stacking sequences by the user.

The two minimize functions address different failure modes independently.
To size a panel, run both and take the result with the larger thickness T — that
is the governing constraint.


Maximize Buckling Load
----------------------
.. autofunction:: composipy.optimize.maximize_buckling_load

Minimize Panel Weight (Buckling)
---------------------------------
.. autofunction:: composipy.optimize.minimize_panel_weight

Minimize Panel Weight (Strain)
--------------------------------
.. autofunction:: composipy.optimize.minimize_panel_weight_strain
