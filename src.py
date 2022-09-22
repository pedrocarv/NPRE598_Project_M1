import numpy as np
import matplotlib.pyplot as plt
import boris
import mirror
import maxwell

kB = 1.38e-23
amu = 1.66e-27
m = 85*amu
T = 298.15
qe = 1.60217662e-19 # eV
Q = -qe

def main():
    # R, Z grid
    X = np.linspace( -1, 1, 20 )
    Y = np.linspace( -1, 1, 20 )
    Z = np.linspace( -1, 1, 20 )

    # Cell size
    dY = Y[1] - Y[0]
    dZ = Z[1] - Z[0]

    Ymin = Y[0]
    Ymax = Y[-1]
    Zmin = Z[0]
    Zmax = Z[-1]

    Bx_grid = np.zeros((Y.size, Z.size))
    By_grid = np.zeros((Y.size, Z.size))
    Bz_grid = np.zeros((Y.size, Z.size))
    Bnorm   = np.zeros((Y.size, Z.size))

    for j in range(0,Y.size):
        for k in range(0,Z.size):
                Bx, By, Bz = mirror.mirror(0.0, Y[j], Z[k])
                Bx_grid[j,k] = Bx
                By_grid[j,k] = By
                Bz_grid[j,k] = Bz
                Bnorm[j,k]   = np.sqrt( Bx*Bx + By*By + Bz*Bz )

    Bmin = np.min(Bnorm)
    Bmax = np.max(Bnorm)

    theta = np.arcsin(np.sqrt(Bmin/Bmax))

    # Initial conditions
    v = np.arange(0, 2)
    cdf = maxwell.CDF(v,m,T)
    nparticles = 10
    speeds, vx, vy, vz = maxwell.get_vel(nparticles)

    X0 = 0.0
    Y0 = 0.0
    Z0 = 0.0

    Nperiods = 10
    Nsteps_per_period = 100
    time = np.linspace(0.0, 1.0, Nperiods*Nsteps_per_period)
    dt = time[1] - time[0]
    params = np.array([dt, Q/m*dt/2])

    particleTraj = np.zeros((nparticles, 6, time.size))

    for i in range(nparticles):
        Vx0 = vx[i]
        Vy0 = vy[i]   
        Vz0 = vz[i]
        x0 = np.array([X0, Y0, Z0, Vx0, Vy0, Vz0])
        #Traj = boris.boris_bunemann(time, x0, params, Ymin, Ymax, Zmin, Zmax, dY, dZ, Bx_grid, By_grid, Bz_grid, correction=True)
        Traj = boris.boris_bunemann(time, x0, params, correction=True)
        for t in range(time.size):
            particleTraj[i, 0, t] = Traj[t,0] 
            particleTraj[i, 1, t] = Traj[t,1]
            particleTraj[i, 2, t] = Traj[t,2]
            particleTraj[i, 3, t] = Traj[t,3]
            particleTraj[i, 4, t] = Traj[t,4]
            particleTraj[i, 5, t] = Traj[t,5]

    # plt.figure(1, figsize=(9,9))
    # ax = plt.subplot(211)
    # XX,YY = np.meshgrid(X,Y)
    # plt.contourf(XX,YY,Bnorm[:,:,4],50)
    # plt.colorbar()
    # plt.xlabel('X [m]')
    # plt.ylabel('Y [m]')
    # plt.title('B-field magnitude [T] - Simple Magnetic Mirror')

    # plt.figure(1)
    # ax = plt.subplot(312)
    # XX,YY = np.meshgrid(X,Y)
    # plt.streamplot(XX, YY, Bx_grid[:,:,10], By_grid[:,:,10])
    # plt.xlabel('X [m]')
    # plt.ylabel('Y [m]')
    # #plt.title('B-field components - Simple Magnetic Mirror')

    # plt.figure(1)
    # ax = plt.subplot(212)
    # XX,YY = np.meshgrid(Z, Y)
    # ZY = plt.contourf(XX,YY,Bnorm[25,:,:],50)
    # plt.colorbar()
    # plt.xlabel('Z [m]')
    # plt.ylabel('Y [m]')


    # XX,YY = np.meshgrid(Y, Z)
    # plt.contourf(np.transpose(XX),np.transpose(YY),Bnorm,50)
    # plt.colorbar()
    # plt.xlabel('Y [m]')
    # plt.ylabel('Z [m]')
    # #plt.title('B-field magnitude [T] - Simple Magnetic Mirror')
    

    plt.subplots()
    plt.plot(particleTraj[0,2,:], particleTraj[0,1,:], 'r-')
    plt.plot(particleTraj[1,2,:], particleTraj[1,1,:], 'b-')
    plt.plot(particleTraj[2,2,:], particleTraj[2,1,:], 'g-')
    plt.xlabel('Z [m]')
    plt.ylabel('Y [m]')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.title('Multiple particle test')
    plt.show()

if __name__ == '__main__':
   main()