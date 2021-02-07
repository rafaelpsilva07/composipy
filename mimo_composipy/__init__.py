from .ply_class import Ply
from .critical_buckling_function import critical_buckling
from .laminate_class import Laminate
from .buckling_load_function import buckling_load

if __name__ =='__main__':
    print('this script is being executed')
    ply1 = Ply(125500, 9450, 0.32, 4700, 0.35)
    layup_1 = [(45,ply1),(0,ply1),(0,ply1),(45,ply1),(0,ply1),(0,ply1),(45,ply1)]
    plate_1 = Laminate(layup_1)
    critical_buckling(360,360,plate_1.D, method='polynomial')