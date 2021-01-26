from .ply_class import Ply
from .critical_buckling_function import critical_buckling
from .laminate_class import Laminate


if __name__ =='__main__':
    ply1 = Ply(125500, 9450, 0.32, 4700, 0.35)
    print(ply1.Q_0)