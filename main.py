import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

# Program params
size = 100  # Plot axis size (x and y length of the grid)
dt = 0.01  # Time step
diffusion = 0.00005  # Diffusion (how quickly the fluid density spreads out)
viscosity = 0.00005  # Viscosity (how quickly the velocity field smooths out)
iter = 8  # Iterations for linear solver

# fluid properties
class Fluid:
    def __init__(self, diffusion, viscosity, dt):
        self.size = size
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity
        
        self.s = np.zeros((size, size), dtype=float) # temp array for calculations
        self.density = np.zeros((size, size), dtype=float) # density field
        
        self.Vx = np.zeros((size, size), dtype=float)  # h velocity
        self.Vy = np.zeros((size, size), dtype=float)  # v velocity
        
        # temp arrays for velocity
        self.Vx0 = np.zeros((size, size), dtype=float) 
        self.Vy0 = np.zeros((size, size), dtype=float)
    
    # updates variables and calls all methods to update at each time step
    def step(self):
        # updating vars
        visc = self.visc
        diff = self.diff
        dt = self.dt
        Vx = self.Vx
        Vy = self.Vy
        Vx0 = self.Vx0
        Vy0 = self.Vy0
        s = self.s
        density = self.density

        # projection ensures incompressibility
        # diffuse smooths out velocity and density fields
        # advection moves density and velocity along the velocity field

        diffuse(1, Vx0, Vx, visc, dt)
        diffuse(2, Vy0, Vy, visc, dt)

        projection(Vx0, Vy0, Vx, Vy) 

        advection(1, Vx, Vx0, Vx0, Vy0, dt)
        advection(2, Vy, Vy0, Vx0, Vy0, dt)

        projection(Vx, Vy, Vx0, Vy0)

        diffuse(0, s, density, diff, dt)
        advection(0, density, s, Vx, Vy, dt)
    
    # changing density 
    def addDensity(self, x, y, increment):
        self.density[x, y] += increment
    
    # changing velocity
    def addVelocity(self, x, y, incrementX, incrementY):
        self.Vx[x, y] += incrementX
        self.Vy[x, y] += incrementY
    
    # reducing density over time to emulate fluid dissipation
    def fade(self):
        fade_rate = 0.985  # fade rate to make fluid dissipate gradually
        self.density *= fade_rate
        self.density = np.clip(self.density, 0, 1)

# Fluid simulation functions

# Solves system of equations using Gauss-Seidel method to solve for velocity / density fields during diffusion / projection
def linear_solve(b, x, x0, a, c):
    cReciprocal = 1.0 / c
    # iteratively updates field x based on averages of neighboring cells
    for _ in range(iter):
        x[1:-1, 1:-1] = (x0[1:-1, 1:-1] + a * (x[2:, 1:-1] + x[:-2, 1:-1] + x[1:-1, 2:] + x[1:-1, :-2])) * cReciprocal
    # enforce boundary conditions
    set_bnd(b, x)

# setting plot boundaries for computation
def set_bnd(b, x):
    x[0, :] = x[1, :]
    x[-1, :] = x[-2, :]
    x[:, 0] = x[:, 1]
    x[:, -1] = x[:, -2]

    # handling edge cases (top of plot)
    if b == 1: # h velocity
        x[0, :] *= -1
        x[-1, :] *= -1
    elif b == 2:  # v velocity
        x[:, 0] *= -1
        x[:, -1] *= -1

# simulates diffusion of velocity / density field
# b - boundary condition type : 0 for density, 1 for h velocity, 2 for v velocity
# diff - diffusion rate
def diffuse(b, x, x0, diff, dt):
    # coeff a determines how much field diffuses in a single step
    a = dt * diff * (size - 2) * (size - 2)
    linear_solve(b, x, x0, a, 1 + 6 * a)

# simulates diffusion of velocity / density field
# b - boundary condition type : 0 for density, 1 for h velocity, 2 for v velocity
def advection(b, d, d0, velocX, velocY, dt):
    # calculates where fluid at cell comes from in prev step
    dtx = dt * (size - 2)
    dty = dt * (size - 2)
    
    # bilinear search for cell value is previous cell
    for j in range(1, size - 1):
        for i in range(1, size - 1):
            x = i - dtx * velocX[i, j]
            y = j - dty * velocY[i, j]
            
            x = max(0.5, min(size - 1.5, x))
            y = max(0.5, min(size - 1.5, y))
            
            i0 = int(x)
            i1 = i0 + 1
            j0 = int(y)
            j1 = j0 + 1
            
            s1 = x - i0
            s0 = 1 - s1
            t1 = y - j0
            t0 = 1 - t1
            
            d[i, j] = (
                s0 * (t0 * d0[i0, j0] + t1 * d0[i0, j1]) +
                s1 * (t0 * d0[i1, j0] + t1 * d0[i1, j1])
            )
    # enforce boundary conditions
    set_bnd(b, d)

# ensures velocity field is incompressible ((u * vector field) * u = 0)
# p - temp array to store pressure field
# div - temp array to store divergence of velocity field 
def projection(velocX, velocY, p, div):
    # how much velocity field is spreading out
    div[1:-1, 1:-1] = -0.5 * (velocX[2:, 1:-1] - velocX[:-2, 1:-1] + velocY[1:-1, 2:] - velocY[1:-1, :-2]) / size
    p[1:-1, 1:-1] = 0

    set_bnd(0, div)
    set_bnd(0, p)

    # solves for pressure field p to correct velocity field
    linear_solve(0, p, div, 1, 6)

    # velocty field updated by subtracting gradient of pressure field
    velocX[1:-1, 1:-1] -= 0.5 * (p[2:, 1:-1] - p[:-2, 1:-1]) * size
    velocY[1:-1, 1:-1] -= 0.5 * (p[1:-1, 2:] - p[1:-1, :-2]) * size

    # enforce boundary conditions
    set_bnd(1, velocX)
    set_bnd(2, velocY)

# Initialize fluid simulation
fluid = Fluid(diffusion, viscosity, dt)

# fluid source location
source_x, source_y = size // 2, 10  # Move source to bottom center
source_radius = 5  

# Oscillating angle variation
normal_theta = math.pi / 2  # Straight up
variation_speed = 0.1  # Speed of oscillation
variation_magnitude = math.pi / 12  # Max angle deviation

# Set up plot
fig, ax = plt.subplots()

# replace arg value assigned to "cmap=" to ...
#  'gist_gray' for smoke colors
#  'jet' for traditional spectral colors
# Full list of colormaps at https://matplotlib.org/stable/users/explain/colors/colormaps.html

im = ax.imshow(fluid.density, cmap='gist_gray', origin='lower', vmin=0, vmax=1)
plt.colorbar(im, label='Fluid Density')

# Update func for animation
def update(frame):
    global fluid, normal_theta
    
    # Oscillate the angle over time for smoother changes
    angle_offset = variation_magnitude * math.sin(frame * variation_speed)
    curr_theta = normal_theta + angle_offset
    
    # Velocity components
    velocity_x = math.sin(curr_theta) * 2
    velocity_y = math.cos(curr_theta) * 0.8
    
    for i in range(-source_radius, source_radius + 1):
        for j in range(-source_radius, source_radius + 1):
            if (i * i + j * j) <= source_radius * source_radius:
                fluid.addDensity(source_y + i, source_x + j, 1.0)
                fluid.addVelocity(source_y + i, source_x + j, velocity_x, velocity_y)
    
    fluid.step()
    fluid.fade()
    im.set_array(fluid.density)
    return im,

# Run animation and show plot
ani = FuncAnimation(fig, update, frames = 60, interval = 16, blit = True)
plt.show()