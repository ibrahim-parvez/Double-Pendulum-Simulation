# Double-Pendulum-Simulation

This C++ program simulates the motion of a double pendulum using the Runge-Kutta method for solving the system of ordinary differential equations (ODEs) governing the motion of the pendulums. The program allows users to input parameters such as the length of the rods, masses of the pendulums, initial angles, and velocities. It also includes the option to simulate damping (air resistance) for a more realistic scenario.

Features
Double Pendulum Simulation: Simulates the motion of a two-mass, two-rod system.
Runge-Kutta Method: Solves the ODEs governing the motion of the system numerically.
Damping Effect: Includes an option to add air resistance that gradually slows down the motion over time.
Energy Conservation: Models the system's potential energy (PE), kinetic energy (KE), and total energy (TE) to ensure that the simulation conserves momentum.
Equations of Motion
The motion of the double pendulum is governed by a set of second-order ordinary differential equations (ODEs). These can be solved using the Runge-Kutta method. The key equations used are:

Without Damping:

θ1'' and θ2'' represent the angular accelerations of the two masses.
The equations include gravitational forces, interaction forces between the pendulums, and lengths of the rods.
With Damping:

The equations are modified to include a damping term proportional to the angular velocity of the pendulums, simulating air resistance.
Energy Equations:

Potential Energy (PE): The potential energy of the pendulums due to gravity.
Kinetic Energy (KE): The kinetic energy of both masses based on their velocities.
Total Energy (TE): The sum of potential and kinetic energy, used for verification of the simulation.
