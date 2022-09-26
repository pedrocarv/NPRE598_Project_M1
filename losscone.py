import numpy as np
import matplotlib.pyplot as plt
import mirror

def unit_vec(vec):
    return vec / np.linalg.norm(vec)

def fieldline_tracing(x0, y0, z0, dx, dy, dz, nstep, whole=False):

    x = x0
    y = y0
    z = z0
    Z = [z]
    Y = [y]
    X = [x]
    Bmax = 0.0
    Bmin = 0.0

    for i in range(0, nstep):
        Bx, By, Bz = mirror.mirror(x,y,z)
        Bnorm = np.sqrt(Bx*Bx + By*By + Bz*Bz)

        if i == 0:
            Bmin = Bnorm
        if Bnorm > Bmax:
            Bmax = Bnorm

        x += dx*Bx/Bnorm
        y += dy*By/Bnorm
        z += dz*Bz/Bnorm
        X = np.append(X,x)
        Y = np.append(Y,y)
        Z = np.append(Z,z)

    if whole == True: # This reflects the field line in the other direction (due to symmetry)
        X = np.append(X[::-1], X)
        Y = np.append(Y[::-1], Y)
        Z = np.append(Z[::-1], -Z)

    return X, Y, Z, Bmax, Bmin

def particle_trapping(X, Y, Z, x0, y0, Vx, Vy, Vz, Bmax, Bmin):

    # This assumes that at the mid plane the Vz component is the only parallel component to the magnetic field
    angle = np.arctan2(np.sqrt(Vx**2 + Vy**2),Vz)
    if angle > np.pi/2:
        angle = np.pi - angle

    alpha = np.arctan2(np.sqrt((Bmax-Bmin)), np.sqrt(Bmin)) # loss-cone angle
    
    if angle > alpha:
        return 1, angle, alpha # Particle is trapped
    else:
        return 0, angle, alpha # Particle escaped 

def test_pitch():
    Bx0, By0, Bz0 = mirror.mirror(0,0,0)
    B = [Bx0, By0, Bz0]
    V = [0.1, 0.1, 0.0]
    b = unit_vec(B)
    v = unit_vec(V)
    angle = np.arccos(np.dot(v, b))
    print('test pitch angle =', angle*180/np.pi)
    return 


