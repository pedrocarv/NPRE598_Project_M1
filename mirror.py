import numpy as np
import bfield

def mirror(x, y, z):
    '''
    This function calculates the cartesian components of the magnetic field generated by two current loops.
    
    Input:
    x, y, x --> cartesian coordinates of the point where the field will be evaluated

    Output:
    Bx, By, Bz --> cartesian components of the magnetic field
    
    '''
    Br, Bz = [0.0, 0.0]

    Loops = np.array([[0.25, 5000., 5], 
                      [0.25, 5000., 5]]) # [Ra, I0, Nturns] of the loops\

    Center = np.array([[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]])

    for l in range(2):
        Ra     = Loops[l][0]
        I0     = Loops[l][1]
        Nturns = Loops[l][2]
        z1 = z - Center[l][2]
        r = np.sqrt(x*x + y*y)
        phi = np.arctan2(y,x)
        Br_temp, Bz_temp = bfield.loopbrz(Ra, I0, Nturns, r, z1)
        Br += Br_temp
        Bz += Bz_temp

    Bx = Br*np.cos(phi)
    By = Br*np.sin(phi)
    Bz = Bz

    return Bx, By, Bz

def Bfield_interpolator(x, y, xmin, xmax, ymin, ymax, dx, dy, Bx_grid, By_grid, Bz_grid):
    '''
    This function does 2D linear interpolation of the B-field at any point
    '''
    i = int(np.abs(np.floor((x-xmin)/dx)))
    j = int(np.abs(np.floor((y-ymin)/dy)))
   
    xi = xmin + i * dx
    yj = ymin + j * dy

    A0 = (xi + dx - x)*(yj + dy - y)
    A1 = (x - xi)*(yj + dy - y)
    A2 = (x - xi)*(y - yj)
    A3 = (xi + dx - x)*(y - yj)
    At = dx*dy

    w0 = A0/At
    w1 = A1/At
    w2 = A2/At
    w3 = A3/At

    if xi == xmax:
        if yj == ymax:
            Bx = w0*Bx_grid[i,j] + w1*Bx_grid[i,j] + w2*Bx_grid[i,j] + w3*Bx_grid[i,j]
            By = w0*By_grid[i,j] + w1*By_grid[i,j] + w2*By_grid[i,j] + w3*By_grid[i,j]
            Bz = w0*Bz_grid[i,j] + w1*Bz_grid[i,j] + w2*Bz_grid[i,j] + w3*Bz_grid[i,j]
        elif yj < ymax:
            Bx = w0*Bx_grid[i,j] + w1*Bx_grid[i,j] + w2*Bx_grid[i,j+1] + w3*Bx_grid[i,j+1]
            By = w0*By_grid[i,j] + w1*By_grid[i,j] + w2*By_grid[i,j+1] + w3*By_grid[i,j+1]
            Bz = w0*Bz_grid[i,j] + w1*Bz_grid[i,j] + w2*Bz_grid[i,j+1] + w3*Bz_grid[i,j+1]
    elif yj == ymax and xi < xmax:
            Bx = w0*Bx_grid[i,j] + w1*Bx_grid[i+1,j] + w2*Bx_grid[i+1,j] + w3*Bx_grid[i,j]
            By = w0*By_grid[i,j] + w1*By_grid[i+1,j] + w2*By_grid[i+1,j] + w3*By_grid[i,j]
            Bz = w0*Bz_grid[i,j] + w1*Bz_grid[i+1,j] + w2*Bz_grid[i+1,j] + w3*Bz_grid[i,j]
    elif xi > xmax:
            Bx = w0*Bx_grid[i,j] + w1*Bx_grid[i-1,j] + w2*Bx_grid[i-1,j] + w3*Bx_grid[i,j]
            By = w0*By_grid[i,j] + w1*By_grid[i-1,j] + w2*By_grid[i-1,j] + w3*By_grid[i,j]
            Bz = w0*Bz_grid[i,j] + w1*Bz_grid[i-1,j] + w2*Bz_grid[i-1,j] + w3*Bz_grid[i,j]

    if xi < xmax and yj < ymax:
        Bx = w0*Bx_grid[i,j] + w1*Bx_grid[i+1,j] + w2*Bx_grid[i+1,j+1] + w3*Bx_grid[i,j+1]
        By = w0*By_grid[i,j] + w1*By_grid[i+1,j] + w2*By_grid[i+1,j+1] + w3*By_grid[i,j+1]
        Bz = w0*Bz_grid[i,j] + w1*Bz_grid[i+1,j] + w2*Bz_grid[i+1,j+1] + w3*Bz_grid[i,j+1]
    
    return Bx, By, Bz

