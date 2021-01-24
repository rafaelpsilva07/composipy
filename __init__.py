from src.ply_class import Ply
from src.critcal_buckling_function import critcal_buckling
from src.laminate_class import Laminate


if __name__ =='__main__':
    ply1 = Ply(125500, 9450, 0.32, 4700, 0.35)
    print(ply1.Q_0)