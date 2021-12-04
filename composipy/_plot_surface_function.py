
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import sympy as sp

def _plot_surface(a,b, C_values, Nw_values ,n_divisions = 25):

    #Data
    a_array = np.linspace(0,a,n_divisions)
    b_array = np.linspace(0,b,n_divisions)
    a_array, b_array = np.meshgrid(a_array, b_array)
    z = np.zeros(n_divisions**2).reshape(n_divisions,n_divisions)
    x = sp.symbols('x')
    y = sp.symbols('y')

    for i in range(n_divisions):
        for j in range(n_divisions): 
            z[i,j] =  np.array(Nw_values.subs(x,a_array[i,j]).subs(y,b_array[i,j])).astype('float')@C_values
    
    #Plotting
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot the surface.
    surf = ax.plot_surface(a_array, b_array, z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)


    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()