# NPRE598_Project_M1

PROJECT: Loss cone of a magnetic mirror for a Maxwell-Boltzmann ion population

• First, prepare a routine returning the 3 components of the B-field of a 2-
loop magnetic mirror, e.g.: 
                 Bx, By, Bz = bfield.mirror(x,y,z)

• Sample the B-field on the R,Z plane of the mirror, and write a function 
returning the 3 interpolated components of the B-field at a generic point 
P(x,y,z). This will greatly improve computational speed w.r.t. a direct B-field 
calculation.  

• Solve the Newton-Lorentz equation of an ion moving in the grid; do not 
solve the B-field at every point, but rather use the interpolator; initialize ions 
using a Maxwell-Boltzmann distribution at the midplane (green line) and 
integrate their trajectory over time. 

• Establish a criterion to decide automatically if the particle is ‘trapped’ within 
the mirror, or ‘lost’ into the loss cone. 

• Repeat the calculation for a large number of particles, O(104), and 
reconstruct the distribution function 

• Compare the loss-cone you find from your calculation against the 
theoretical value.  

• Optional: repeat the same calculation adding an electrostatic perturbation 
near the cyclotron frequency, eφ(z, t) = eφ0 sin(k⊥x − ω0t), ω0 = Ω(±z0), 
where Ω is the cyclotron frequency at the midpoint, and study the onset of 
chaotic motion as a function of the ES wave amplitude. 
