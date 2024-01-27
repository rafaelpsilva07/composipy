Composipy Optimization Functions
================================
This page contains functiosn used to optimize the composite plate.
Optimization algorithms are gradient based and they use Scipy and optimization engine.

Since optimization is in continuous they result in Lamination Parameters as result. 
Then, Lamination Parameters are supposed to be convereted into stacking sequences somehow by the user.


Maximize Buckling Load
----------------------
.. autofunction:: composipy.optimize.maximize_buckling_load

Minimize panel weight
----------------------
.. autofunction:: composipy.optimize.minimize_panel_weight
