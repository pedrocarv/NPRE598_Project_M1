import numpy as np

def unit_vec(vec):
    return vec / np.linalg.norm(vec)

def trapping(Bmax, Bmin,Vx, Vy, Vz):
    mirror_ratio = Bmax/Bmin
    alpha = np.arcsin(np.sqrt(1/mirror_ratio))
    #B = [Bx, By, Bz]
    V = [Vx, Vy, Vz]
    Vperp = np.sqrt(Vx**2 + Vy**2)
    Vpar = np.abs(Vz)
    if Vperp**2 >= Vpar**2/mirror_ratio:
        return 1
    else:
        return 0
    # b = unit_vec(B)
    # v = unit_vec(V)
    # angle = np.arccos(np.dot(v, b))
    # print(angle, alpha)
    # if angle < alpha:
    #     return 1
    # else:
    #     return 0