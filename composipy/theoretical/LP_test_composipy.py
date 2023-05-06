from composipy import OrthotropicMaterial, LaminateProperty, PlateStructure

a = 200.
b = 200.
T = 1.


E1 = 60000.
E2 = 7000.
v12 = 0.07
G12 = 400
t = 0.21 # Not used

mat1 = OrthotropicMaterial(E1, E2, v12, G12, t)
laminate1 = LaminateProperty(stacking={'xiD': [0., 0., -1., 0.], 'T': 1.0}, plies=mat1)
plate1 = PlateStructure(laminate1, a, b, m=10, n=10, Nxx=-1, constraints="PINNED")


print(plate1.buckling_analysis())
plate1.plot_eigenvalue()



