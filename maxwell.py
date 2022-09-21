import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d

# To generate random velocity components that obey a maxwellian distribution, we need the cumulative distribution function (CDF), which can be found from the PDF:
#                                                        C_v(v) = \int_{-\inf}^{v} f_v(v') \; {\rm d}v'
# Then, to get the velocity components, we need to invert the CDF. We'll do it numerically using an interpolation function.

def MB_speed(v,m,T):
    """ Maxwell-Boltzmann speed distribution for speeds """
    kB = 1.38e-23
    return (m/(2*np.pi*kB*T))**1.5 * 4*np.pi * v**2 * np.exp(-m*v**2/(2*kB*T))

def CDF(v, m, T):
    '''
    This function returns the cumulative distribution of a Maxwell-Boltzmann velocity distribution in 3D
    '''
    k = 1.380649e-23 # J/K 
    ffactor = np.sqrt(k*T/m)
    CDF = erf(v/(np.sqrt(2)*ffactor)) - np.sqrt(2/np.pi) * v * np.exp(-v**2/(2*ffactor**2))/ffactor

    return CDF

def get_vel(CDF, v, n):
    '''
    This function uses scipy interpolation to transform the CDF into a function that can be called to generate the 3D velocity components
    '''
    CDFinv = interp1d(CDF, v)
    speeds = CDFinv(np.random.random(n))

    # Generate random angle for vectors in spherical coordinates uniformly
    theta = np.arccos(np.random.uniform(-1,1,n))
    phi = np.random.uniform(0,2*np.pi,n)

    # Converting to cartesian coordinates
    vx = speeds * np.sin(theta) * np.cos(phi) 
    vy = speeds * np.sin(theta) * np.sin(phi)
    vz = speeds * np.cos(theta)

    return speeds, vx, vy, vz

def main():
    kB = 1.38e-23
    amu = 1.66e-27
    mass = 85*amu
    T = 400
    vs = np.arange(0,1500)

    fv = MB_speed(vs,mass,T)
    cdf = CDF(vs,mass,T)

    spd, vx, vy, vz = get_vel(cdf,vs,100000)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #generate histogram of velocities
    ax.hist(spd,bins=50,density=True, stacked=True,fc='b',alpha=0.4,lw=0.2)

    #compare this histogram to f(v) - this is MB_speed that we wrote earlier

    ax.plot(vs,fv,'k')

    ax.set_xlabel('Speed (m/s)')
    ax.set_ylabel('PDF')

    plt.show()

if __name__ == '__main__':
    main()