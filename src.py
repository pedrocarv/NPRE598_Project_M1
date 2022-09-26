import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib.pyplot as plt
import boris
import mirror
import maxwell
import losscone

kB = 1.38e-23
m = 1.66e-27 # kg
T = 5000 # K
qe = 1.60217662e-19 # eV
Q = qe

def main():
    # R, Z grid
    X = np.linspace( -1, 1, 101 )
    Y = np.linspace( -1, 1, 101 )
    Z = np.linspace( -1, 1, 101 )

    # Cell size
    dY = Y[1] - Y[0]
    dZ = Z[1] - Z[0]

    Ymin = Y[0]
    Ymax = Y[-1]
    Zmin = Z[0]
    Zmax = Z[-1]

    # Calculating the Bfield on the grid points

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
    
    # Initial conditions for the velocities
    nparticles = 10000
    speeds, vx, vy, vz = maxwell.get_vel(nparticles, T)

    # Randomizing the initial positions along the Y-Z mid plane
    r_rand = np.ones(nparticles)*0.3
    p_rand = np.random.rand(nparticles)*2*np.pi
    X0 = np.empty(nparticles)
    Y0 = np.empty(nparticles)
    Z0 = 0.0

    for i in range(p_rand.size):
        X0[i] = r_rand[i]*np.cos(p_rand[i])
        Y0[i] = r_rand[i]*np.sin(p_rand[i])

    Nperiods = 100
    Nsteps_per_period = 100
    # Cyclotron frequency [rad/s]
    omega_c = np.abs(Q)*Bnorm[int((Y.size-1)/2), int((Z.size-1)/2)]/m
    # Cyclotron period
    Tc = 2.0*np.pi/omega_c

    time = np.linspace(0.0, Tc*Nperiods, Nperiods*Nsteps_per_period)
    dt = time[1] - time[0]
    params = np.array([dt, Q/m*dt/2])

    particleTraj = np.empty((nparticles, 6, time.size))
    particleTraj.fill(np.nan)

    trapped_particles = 0
    theory_trapped_particles = 0

    # Arrays used for field line tracing

    nstep = 100
    X_line = np.empty((nparticles, 101))
    Y_line = np.empty((nparticles, 101))
    Z_line = np.empty((nparticles, 101))

    # Arrays to store angles

    pitchangle = np.empty(nparticles)
    lossconeAngle = np.empty(nparticles)


    # Reconstruction of the velocity distribution
    # fig, axs = plt.subplots(1)
    # N, bins, patches = plt.hist(speeds, bins=100, density=True, alpha=0.5, histtype='bar', ec='black')
    # plt.title('Reconstruction of the velocity distribution')
    # plt.ylabel('Frequency')
    # plt.xlabel('Velocity [$km/s$]')
    # plt.tight_layout()
    # plt.savefig('Reconstruction.eps')
    # plt.show()

    for i in range(nparticles):
        X_line[i,:], Y_line[i,:], Z_line[i,:], Bmax, Bmin = losscone.fieldline_tracing(X0[i], Y0[i], 0.0, dY, dY, dZ, nstep)
        theoryTrapped, pitchangle[i], lossconeAngle[i] = losscone.particle_trapping(X_line[i,:], Y_line[i,:], Z_line[i,:], X0[i], Y0[i], vx[i], vy[i], vz[i], Bmax, Bmin)
        x0 = np.array([X0[i], Y0[i], Z0, vx[i], vy[i], vz[i]])
        Traj, isTrapped = boris.boris_bunemann(time, x0, params, Ymin, Ymax, Zmin, Zmax, dY, dZ, Bx_grid, By_grid, Bz_grid, correction=True)
        trapped_particles += isTrapped
        theory_trapped_particles += theoryTrapped

        for t in range(time.size):
            particleTraj[i, 0, t] = Traj[t,0] 
            particleTraj[i, 1, t] = Traj[t,1]
            particleTraj[i, 2, t] = Traj[t,2]
            particleTraj[i, 3, t] = Traj[t,3]
            particleTraj[i, 4, t] = Traj[t,4]
            particleTraj[i, 5, t] = Traj[t,5]

    print(trapped_particles, "particles are trapped")
    print(theory_trapped_particles, "particles are trapped (theory)")

    vperp = [0, np.max(np.sqrt(vx**2 + vy**2))]
    vpar = [0, np.max(np.sqrt(vx**2 + vy**2))/np.tan(lossconeAngle[0])]

    # fig, axs = plt.subplots(1)
    # plt.scatter(abs(vz), np.sqrt(vx**2 + vy**2), s = 6)
    # plt.plot(vpar, vperp, color='r', label='Loss cone')
    # plt.title('Loss cone for 10,000 particles')
    # plt.ylabel('$v_{\perp}$')
    # plt.xlabel('$v_{\parallel}$')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig('Losscone.eps')
    # plt.show()

    # plt.figure(2)
    # #plt.plot(np.transpose(Z_line), np.transpose(Y_line), 'r-')
    # plt.axvline(x=-0.5, ymin=0.25, ymax= 0.75, color='k', linewidth = 2, label='Boundary')
    # plt.plot(np.transpose(particleTraj[:,2,:]), np.transpose(particleTraj[:,1,:]), linewidth = 1 )
    # #plt.axvline(x=-0.5, ymin=0.375, ymax= 0.625, color='b', linewidth = 4, label='Current loop')
    # plt.axvline(x=0.5, ymin=0.25, ymax= 0.75, color='k', linewidth = 2)
    # #plt.axvline(x=0.5, ymin=0.375, ymax= 0.625, color='b', linewidth = 4)
    # plt.axhline(y=0.5, xmin=0.25, xmax=0.75,color='k', linewidth = 2)
    # plt.axhline(y=-0.5,  xmin=0.25, xmax=0.75, color='k', linewidth = 2)
    # plt.xlabel('Z [m]')
    # plt.ylabel('Y [m]')
    # plt.title('Boundaries for the loss criterion ($YZ$-plane)')
    # plt.legend()
    # plt.xlim(-1,1)
    # plt.ylim(-1,1)
    # plt.savefig('10kparticles.png', dpi = 1200)
    # plt.show()

    # plt.figure(1)
    # #XX,YY = np.meshgrid(Z,Y)
    # #plt.streamplot(XX, YY, Bz_grid[:,:], By_grid[:,:])
    # #plt.plot(np.transpose(Z_line), np.transpose(Y_line), 'r-')
    # plt.plot(np.transpose(particleTraj[:,2,:]), np.transpose(particleTraj[:,1,:]), linewidth = 1 )
    # plt.axvline(x=-0.5, ymin=0.375, ymax= 0.625, color='k', linewidth = 3, label='Current loop')
    # plt.axvline(x=0.5, ymin=0.375, ymax= 0.625, color='k', linewidth = 3)
    # plt.xlabel('Z [m]')
    # plt.ylabel('Y [m]')
    # plt.title('Trajectory projection on $YZ$-plane')
    # plt.legend()
    # plt.xlim(-1,1)
    # plt.ylim(-1,1)
    # plt.savefig('100_particles.eps')
    # plt.show()

    # ax = plt.subplot()
    # ZZ, YY = np.meshgrid(Z,Y)
    # plt.contourf(ZZ, YY, Bnorm, 50)
    # plt.xlabel('Z [m]')
    # plt.ylabel('Y [m]')
    # plt.title('Magnetic field magnitude [T] on $YZ$-plane')
    # plt.colorbar()
    # plt.xlim(-1,1)
    # plt.ylim(-1,1)
    # plt.savefig("BFieldMag.eps")
    # plt.show()

if __name__ == '__main__':
   main()