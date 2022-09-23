import numpy as np
from scipy.stats import maxwell

def get_vel(n, T):
    '''
    This function uses inverse transform sampling to get the 3D velocity components
    '''
    k = 1.380649e-23
    m = 1.66e-27
    uniform_rand = np.random.rand(n)

    speeds = maxwell.ppf(uniform_rand, scale=np.sqrt(k*T/m))

    # Generate random angle for vectors in spherical coordinates uniformly
    theta = np.arccos(np.random.uniform(-1,1,n))
    phi = np.random.uniform(0,2*np.pi,n)

    # Converting to cartesian coordinates
    vx = speeds * np.sin(theta) * np.cos(phi) 
    vy = speeds * np.sin(theta) * np.sin(phi)
    vz = speeds * np.cos(theta)

    return speeds, vx, vy, vz