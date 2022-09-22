import numpy as np
from scipy.stats import maxwell

# To generate random velocity components that obey a maxwellian distribution, we need the cumulative distribution function (CDF), which can be found from the PDF:
#                                                        C_v(v) = \int_{-\inf}^{v} f_v(v') \; {\rm d}v'
# Then, to get the velocity components, we need to invert the CDF. We'll do it numerically using an interpolation function.

def get_vel(n):
    '''
    This function uses scipy interpolation to transform the CDF into a function that can be called to generate the 3D velocity components
    '''
    k = 1.380649e-23
    T = 400
    amu = 1.66e-27
    m = 85*amu
    uniform_rand = np.random.uniform()

    speeds = maxwell.ppf(uniform_rand, scale=np.sqrt(k*T/m))

    # Generate random angle for vectors in spherical coordinates uniformly
    theta = np.arccos(np.random.uniform(-1,1,n))
    phi = np.random.uniform(0,2*np.pi,n)

    # Converting to cartesian coordinates
    vx = speeds * np.sin(theta) * np.cos(phi) 
    vy = speeds * np.sin(theta) * np.sin(phi)
    vz = speeds * np.cos(theta)

    return speeds, vx, vy, vz