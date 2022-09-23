import numpy as np
import matplotlib.pyplot as plt
import mirror

qe = 1.60217662e-19
me = 9.10938356e-31

def frequency_correction(x):
    if x == 0: # To avoid division by zero and taking into account that tan(x) -> 1 as x -> 0
        alpha = 1.0
    else:
        alpha = np.tan(x)/x
    return alpha

def boris_bunemann(time, x0, params,Ymin, Ymax, Zmin, Zmax, dY, dZ, Bx_grid, By_grid, Bz_grid, correction=False):
#def boris_bunemann(time, x0, params, Y, Z, Bx_grid, By_grid, Bz_grid, correction=False):
#def boris_bunemann(time, x0, params, correction=False):

    dt = params[0]
    qmdt2 = params[1]
    N = np.size(time)
    M = np.size(x0)
    X = np.zeros((N,M))
    X[0,:] = x0
    x = X[0,0]
    y = X[0,1]
    z = X[0,2]
    vx = X[0,3]
    vy = X[0,4]
    vz = X[0,5]

    for n in range(0, N): 

        Ex, Ey, Ez = [0.0, 0.0, 0.0]

        Bx, By, Bz = mirror.Bfield_interpolator(y, z, Ymin, Ymax, Zmin, Zmax, dY, dZ, Bx_grid, By_grid, Bz_grid)

        #Bx, By, Bz = mirror.interp2D(Y,Z,Bx_grid,By_grid,Bz_grid,y,z)

        #Bx, By, Bz = mirror.mirror(x, y, z)

        if correction:
            if n == 0:
                alpha = frequency_correction(0.5*qmdt2*np.sqrt(Bx*Bx + By*By + Bz*Bz))
            else:
                alpha = frequency_correction(qmdt2*np.sqrt(Bx*Bx + By*By + Bz*Bz))
        else:
            alpha = 1.

        if n == 0: # Do velocity push back
            # Step 1. Half acceleration (E-field)
            vx -= 0.5 * qmdt2 * Ex * alpha
            vy -= 0.5 * qmdt2 * Ey * alpha
            vz -= 0.5 * qmdt2 * Ez * alpha
            # Step 2. B-field rotation
            tx = -0.5 * qmdt2 * Bx * alpha
            ty = -0.5 * qmdt2 * By * alpha
            tz = -0.5 * qmdt2 * Bz * alpha
        else:
            # Step 1. Half acceleration (E-field)
            vx += qmdt2 * Ex * alpha
            vy += qmdt2 * Ey * alpha
            vz += qmdt2 * Ez * alpha
            # Step 2. B-field rotation
            tx = qmdt2 * Bx * alpha
            ty = qmdt2 * By * alpha
            tz = qmdt2 * Bz * alpha

        tmagsq = tx*tx + ty*ty + tz*tz
        sx = 2.0 * tx/(1 + tmagsq)
        sy = 2.0 * ty/(1 + tmagsq)
        sz = 2.0 * tz/(1 + tmagsq)

        vpx = vx + vy*tz - vz*ty
        vpy = vy + vz*tx - vx*tz
        vpz = vz + vx*ty - vy*tx

        vx += vpy*sz - vpz*sy
        vy += vpz*sx - vpx*sz
        vz += vpx*sy - vpy*sx

        if n == 0:
            # Step 3. Half acceleration (E-field) (half timestep back)
            vx -= 0.5 * qmdt2 * Ex * alpha
            vy -= 0.5 * qmdt2 * Ey * alpha
            vz -= 0.5 * qmdt2 * Ez * alpha
        else:
            # Step 3. Half acceleration (E-field)
            vx += qmdt2 * Ex * alpha
            vy += qmdt2 * Ey * alpha
            vz += qmdt2 * Ez * alpha
            # Step 4. Push position
            x += vx*dt
            y += vy*dt
            z += vz*dt

        X[n,0] = x
        X[n,1] = y
        X[n,2] = z
        X[n,3] = vx
        X[n,4] = vy
        X[n,5] = vz

        if abs(z) > 1.0 or abs(y) > 1.0:
            return X

    return X
    
