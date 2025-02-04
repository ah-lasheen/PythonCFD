# PythonCFD

Computational fluid dynamics simulation using Eulerian methods, NumPy for computation, and Matplotlib for visual representation. This project implements the Navier–Stokes equations with Lagrangian advection to create a real-time 2D fluid simulation.

# Overview

This project simulates fluid dynamics on a 2D grid by solving the Navier–Stokes equations. It uses a stable fluids algorithm that incorporates:
- **Diffusion:** Smooths out the density and velocity fields.
- **Advection:** Transports these fields along the velocity vectors.
- **Projection:** Enforces incompressibility in the fluid by solving for the pressure field.

A dynamic fluid source injects density and momentum (velocity) into the simulation. The source’s velocity direction oscillates over time, adding a natural variation to the flow.

# Features

- **Real-Time Simulation:** Visualize fluid density evolving over time.
- **Stable Fluids Algorithm:** Implements diffusion, advection, and projection steps to simulate incompressible flow.
- **Customizable Parameters:** Easily tweak grid size, time step, diffusion rate, viscosity, and solver iterations.
- **Visualization:** Uses Matplotlib’s `imshow` and `FuncAnimation` for smooth, animated rendering.
- **Dynamic Fluid Source:** A fluid injection point with an oscillating velocity to mimic realistic fluid behavior.

# Prerequisites

Make sure you have Python installed. This project requires the following Python packages:
- [NumPy](https://numpy.org)
- [Matplotlib](https://matplotlib.org/)

# Installation

1. **Clone the Repository**

   Open your terminal and run:
   ```bash
   git clone https://github.com/your-username/PythonCFD.git
   cd PythonCFD
   
2. **Install Dependencies**

   Install the required packages via pip in your termainal:
   ```bash
   pip install numpy matplotlib
   ```

# Usage 

   To run the fluid simulation, execute the Python script:
   ```bash
   python main.py
   ```

# Customization

## Colormap
The simulation uses Matplotlib’s colormap to visualize fluid density. By default, I have set it to 'gist_gray' as my goal of this project was to represent smoke as it would flow from a lit stick of incense:

   ```
   im = ax.imshow(fluid.density, cmap='gist_gray', origin='lower', vmin=0, vmax=1)
   ```
You can change this to another colormap by modifying the cmap parameter. Full list of color maps can be viewed [Here](https://matplotlib.org/stable/users/explain/colors/colormaps.html)

## Fluid Source
The location and behavior of the fluid source can be adjusted by modifying:

   ```
   source_x, source_y = size // 2, 10  # Source location (bottom center)
   source_radius = 5  # Radius of the injection area
   ```

The oscillation in the source velocity is controlled by:

   ```
   normal_theta = math.pi / 2  # Base angle (upwards)
   variation_speed = 0.1         # Oscillation speed
   variation_magnitude = math.pi / 12  # Maximum angle deviation
   ```

# Acknowledgements

This project draws inspiration from the paper [Stable Fluids](https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf) by Jos Stam, which has been widely used in computer graphics and physical simulations to achieve visually stable and realistic fluid motion. Much of the implementation of the Navier-Stokes Equation in my program mimics the explanation and derivation from this paper.
