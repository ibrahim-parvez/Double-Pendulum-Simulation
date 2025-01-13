#include <cmath>
#include <iostream>
#include "cpgplot.h"
#include <tuple>
#include <unistd.h> // For usleep
// Constants
const int nintervals = 5000; // Number of intervals for integration
const float g = 9.81; // Gravity (m/s^2)
const float dt = 0.01; // Time step for numerical integration (s)
// Function to integrate the double pendulum using the Runge-Kutta method
void simulation(float L1, float L2, float m1, float m2, float b, float theta1_0,
float theta2_0, float omega1_0, float omega2_0) {
float t[nintervals + 1];
float theta1[nintervals + 1], theta2[nintervals + 1];
float omega1[nintervals + 1], omega2[nintervals + 1];
float x1[nintervals + 1], y1[nintervals + 1], x2[nintervals + 1],
y2[nintervals + 1];
float KE[nintervals + 1], PE[nintervals + 1], TE[nintervals + 1];
float min_E = 0.0;
// Initial conditions
t[0] = 0.0;
theta1[0] = theta1_0;
omega1[0] = omega1_0;
theta2[0] = theta2_0;
omega2[0] = omega2_0;
// Initial positions of the pendulums
x1[0] = L1 * sin(theta1_0);
y1[0] = -L1 * cos(theta1_0);
x2[0] = x1[0] + L2 * sin(theta2_0);
y2[0] = y1[0] - L2 * cos(theta2_0);
// Trajectory arrays
float traj_x1[nintervals + 1], traj_y1[nintervals + 1];
float traj_x2[nintervals + 1], traj_y2[nintervals + 1];
traj_x1[0] = x1[0];
traj_y1[0] = y1[0];
traj_x2[0] = x2[0];
traj_y2[0] = y2[0];
KE[0] = 0.5 * m1 * (omega1[0] * omega1[0] * L1 * L1) +
0.5 * m2 * ((omega1[0] * omega1[0] * L1 * L1) + (omega2[0] *
omega2[0] * L2 * L2) +
2 * omega1[0] * omega2[0] * L1 * L2 * cos(theta1[0] -
theta2[0]));
PE[0] = - (m1 * g * L1 * cos(theta1[0]) + m2 * g * (L1 * cos(theta1[0]) + L2
* cos(theta2[0])));
TE[0] = KE[0] + PE[0];
// Helper function to compute derivatives (motion equations with damping)
auto calcDerivatives = [&](float theta1, float theta2, float omega1, float
omega2) {
float delta = theta1 - theta2;
// Common denominator
float denom = 2 * m1 + m2 - m2 * cos(2 * delta);
// Angular acceleration of theta1 (with damping)
float num1 = -g * (2 * m1 + m2) * sin(theta1)
- m2 * g * sin(theta1 - 2 * theta2)
- 2 * sin(delta) * m2 * (omega2 * omega2 * L2 + omega1 *
omega1 * L1 * cos(delta))
- b * omega1; // Damping for theta1
float alpha1 = num1 / (L1 * denom);
// Angular acceleration of theta2 (with damping)
float num2 = 2 * sin(delta) * (omega1 * omega1 * L1 * (m1 + m2)
+ g * (m1 + m2) * cos(theta1)
+ omega2 * omega2 * L2 * m2 * cos(delta))
- b * omega2; // Damping for theta2
float alpha2 = num2 / (L2 * denom);
return std::tuple<float, float, float, float>(omega1, alpha1, omega2,
alpha2);
};
// Open the plot window once
cpgopen("/xs");
cpgscr(1, 1.0, 1.0, 1.0); // index 0, background colour
cpgscr(0, .0, .0, .0); // index 1, default fg colour
cpgsci(0);
cpgenv(-L1 - L2 - 1, L1 + L2 + 1,
-L1 - L2 - 1, L1 + L2 + 1, 0, 1);
cpglab("X (meters)", "Y (meters)", "Double Pendulum Trajectory");
// Runge-Kutta integration loop
for (int i = 1; i <= nintervals; i++) {
cpgbbuf();
t[i] = t[i - 1] + dt;
// Undamped
float k1_theta1, k1_omega1, k1_theta2, k1_omega2;
std::tie(k1_theta1, k1_omega1, k1_theta2, k1_omega2) =
calcDerivatives(theta1[i - 1], theta2[i - 1], omega1[i - 1], omega2[i - 1]);
float k2_theta1, k2_omega1, k2_theta2, k2_omega2;
std::tie(k2_theta1, k2_omega1, k2_theta2, k2_omega2) =
calcDerivatives(theta1[i - 1] + 0.5 * k1_theta1 * dt, theta2[i - 1] + 0.5 *
k1_theta2 * dt, omega1[i - 1] + 0.5 * k1_omega1 * dt, omega2[i - 1] + 0.5 *
k1_omega2 * dt);
float k3_theta1, k3_omega1, k3_theta2, k3_omega2;
std::tie(k3_theta1, k3_omega1, k3_theta2, k3_omega2) =
calcDerivatives(theta1[i - 1] + 0.5 * k2_theta1 * dt, theta2[i - 1] + 0.5 *
k2_theta2 * dt, omega1[i - 1] + 0.5 * k2_omega1 * dt, omega2[i - 1] + 0.5 *
k2_omega2 * dt);
float k4_theta1, k4_omega1, k4_theta2, k4_omega2;
std::tie(k4_theta1, k4_omega1, k4_theta2, k4_omega2) =
calcDerivatives(theta1[i - 1] + k3_theta1 * dt, theta2[i - 1] + k3_theta2 * dt,
omega1[i - 1] + k3_omega1 * dt, omega2[i - 1] + k3_omega2 * dt);
theta1[i] = theta1[i - 1] + (dt / 6.0) * (k1_theta1 + 2.0 * k2_theta1 +
2.0 * k3_theta1 + k4_theta1);
omega1[i] = omega1[i - 1] + (dt / 6.0) * (k1_omega1 + 2.0 * k2_omega1 +
2.0 * k3_omega1 + k4_omega1);
theta2[i] = theta2[i - 1] + (dt / 6.0) * (k1_theta2 + 2.0 * k2_theta2 +
2.0 * k3_theta2 + k4_theta2);
omega2[i] = omega2[i - 1] + (dt / 6.0) * (k1_omega2 + 2.0 * k2_omega2 +
2.0 * k3_omega2 + k4_omega2);
// Compute positions of the pendulums
x1[i] = L1 * sin(theta1[i]);
y1[i] = -L1 * cos(theta1[i]);
x2[i] = x1[i] + L2 * sin(theta2[i]);
y2[i] = y1[i] - L2 * cos(theta2[i]);
// Update trajectory
traj_x1[i] = x1[i];
traj_y1[i] = y1[i];
traj_x2[i] = x2[i];
traj_y2[i] = y2[i];
KE[i] = 0.5 * m1 * (omega1[i] * omega1[i] * L1 * L1) +
0.5 * m2 * ((omega1[i] * omega1[i] * L1 * L1) + (omega2[i] *
omega2[i] * L2 * L2) +
2 * omega1[i] * omega2[i] * L1 * L2 * cos(theta1[i] -
theta2[i]));
PE[i] = - (m1 * g * L1 * cos(theta1[i]) + m2 * g * (L1 * cos(theta1[i]) +
L2 * cos(theta2[i])));
TE[i] = KE[i] + PE[i];
if (PE[i] < min_E) {
min_E = PE[i];
}
// Erase the previous frame
cpgeras();
// Reset the axis window and box
cpgswin(-L1 - L2 - 1, L1 + L2 + 1,
-L1 - L2 - 1, L1 + L2 + 1);
cpgbox("ABCNST", 0, 0, "ABCNST", 0, 0);
cpglab("X (meters)", "Y (meters)", "Double Pendulum Trajectory");
// Plot trajectories
cpgsci(2); // Red for trajectory of mass 1
cpgline(i, traj_x1, traj_y1);
cpgsci(3); // Green for trajectory of mass 2
cpgline(i, traj_x2, traj_y2);
// Plot current positions
cpgsci(2); // Red for mass 1
cpgpt(1, &x1[i], &y1[i], 17);
cpgsci(3); // Green for mass 2
cpgpt(1, &x2[i], &y2[i], 17);
// Draw strings
cpgsci(1); // White color for strings
float origin_x = 0.0, origin_y = 0.0;
float line1_x[2] = {origin_x, x1[i]};
float line1_y[2] = {origin_y, y1[i]};
cpgline(2, line1_x, line1_y); // Line from origin to first mass
float line2_x[2] = {x1[i], x2[i]};
float line2_y[2] = {y1[i], y2[i]};
cpgline(2, line2_x, line2_y); // Line from first mass to second mass
cpgebuf();
usleep(20000); // Slow down for visualization
}
float max_E = 0.0;
for (int i = 0; i <= nintervals; i++) {
PE[i] -= min_E;
TE[i] -= min_E;
if(TE[i] > max_E) {
max_E = TE[i];
}
}
// Energy Plot
cpgopen("/xs");
cpgenv(0.0, t[nintervals], -50, max_E + 50, 0, 1); // Adjust axis to include
all energies
cpglab("Time (s)", "Energy (J)", "Energy of Double Pendulum");
// Plot kinetic energy in red
cpgsci(2); // Set color to red
cpgline(nintervals + 1, t, KE);
// Plot potential energy in blue
cpgsci(4); // Set color to blue
cpgline(nintervals + 1, t, PE);
// Plot total energy in green
cpgsci(3); // Set color to green
cpgline(nintervals + 1, t, TE);
cpgclos(); // Close plot window
}
// Function to integrate the double pendulum (damped and undamped) using Runge-
Kutta method
void graphs(float L1, float L2, float m1, float m2, float b, float theta1_0,
float theta2_0, float omega1_0, float omega2_0) {
float t[nintervals + 1];
float theta1[nintervals + 1], theta2[nintervals + 1];
float omega1[nintervals + 1], omega2[nintervals + 1];
// Initial conditions
t[0] = 0.0;
theta1[0] = theta1_0;
omega1[0] = omega1_0;
theta2[0] = theta2_0;
omega2[0] = omega2_0;
// Helper function to compute derivatives (motion equations with damping)
auto calcDerivatives = [&](float theta1, float theta2, float omega1, float
omega2) {
float delta = theta1 - theta2;
// Common denominator
float denom = 2 * m1 + m2 - m2 * cos(2 * delta);
// Angular acceleration of theta1 (with damping)
float num1 = -g * (2 * m1 + m2) * sin(theta1)
- m2 * g * sin(theta1 - 2 * theta2)
- 2 * sin(delta) * m2 * (omega2 * omega2 * L2 + omega1 *
omega1 * L1 * cos(delta))
- b * omega1; // Damping for theta1
float alpha1 = num1 / (L1 * denom);
// Angular acceleration of theta2 (with damping)
float num2 = 2 * sin(delta) * (omega1 * omega1 * L1 * (m1 + m2)
+ g * (m1 + m2) * cos(theta1)
+ omega2 * omega2 * L2 * m2 * cos(delta))
- b * omega2; // Damping for theta2
float alpha2 = num2 / (L2 * denom);
return std::tuple<float, float, float, float>(omega1, alpha1, omega2,
alpha2);
};
// Runge-Kutta integration
for (int i = 1; i <= nintervals; i++) {
t[i] = t[i - 1] + dt;
float k1_theta1, k1_omega1, k1_theta2, k1_omega2;
std::tie(k1_theta1, k1_omega1, k1_theta2, k1_omega2) =
calcDerivatives(theta1[i - 1], theta2[i - 1], omega1[i - 1], omega2[i - 1]);
float k2_theta1, k2_omega1, k2_theta2, k2_omega2;
std::tie(k2_theta1, k2_omega1, k2_theta2, k2_omega2) =
calcDerivatives(theta1[i - 1] + 0.5 * k1_theta1 * dt, theta2[i - 1] + 0.5 *
k1_theta2 * dt, omega1[i - 1] + 0.5 * k1_omega1 * dt, omega2[i - 1] + 0.5 *
k1_omega2 * dt);
float k3_theta1, k3_omega1, k3_theta2, k3_omega2;
std::tie(k3_theta1, k3_omega1, k3_theta2, k3_omega2) =
calcDerivatives(theta1[i - 1] + 0.5 * k2_theta1 * dt, theta2[i - 1] + 0.5 *
k2_theta2 * dt, omega1[i - 1] + 0.5 * k2_omega1 * dt, omega2[i - 1] + 0.5 *
k2_omega2 * dt);
float k4_theta1, k4_omega1, k4_theta2, k4_omega2;
std::tie(k4_theta1, k4_omega1, k4_theta2, k4_omega2) =
calcDerivatives(theta1[i - 1] + k3_theta1 * dt, theta2[i - 1] + k3_theta2 * dt,
omega1[i - 1] + k3_omega1 * dt, omega2[i - 1] + k3_omega2 * dt);
theta1[i] = theta1[i - 1] + (dt / 6.0) * (k1_theta1 + 2.0 * k2_theta1 +
2.0 * k3_theta1 + k4_theta1);
omega1[i] = omega1[i - 1] + (dt / 6.0) * (k1_omega1 + 2.0 * k2_omega1 +
2.0 * k3_omega1 + k4_omega1);
theta2[i] = theta2[i - 1] + (dt / 6.0) * (k1_theta2 + 2.0 * k2_theta2 +
2.0 * k3_theta2 + k4_theta2);
omega2[i] = omega2[i - 1] + (dt / 6.0) * (k1_omega2 + 2.0 * k2_omega2 +
2.0 * k3_omega2 + k4_omega2);
}
// Displacement and Velocity Graphs for Double Pendulum
// Calculate the min and max values for displacement
float min_theta = theta1[0], max_theta = theta1[0];
for (int i = 0; i <= nintervals; i++) {
if (theta1[i] < min_theta) min_theta = theta1[i];
if (theta2[i] < min_theta) min_theta = theta2[i];
if (theta1[i] > max_theta) max_theta = theta1[i];
if (theta2[i] > max_theta) max_theta = theta2[i];
}
// Add some padding for better visualization
float padding_theta = (max_theta - min_theta) * 0.1;
min_theta -= padding_theta;
max_theta += padding_theta;
// -------------------------- Graph 1: Displacement -------------------------
--
cpgopen("/xs");
cpgenv(0.0, t[nintervals], min_theta, max_theta, 0, 1); // Time vs. Angle
cpglab("Time (s)", "Theta (radians)", "Displacement of Pendulums");
// Displacement of Pendulum 1 (Red)
cpgsci(2); // Red color
cpgline(nintervals + 1, t, theta1);
// Displacement of Pendulum 2 (Green)
cpgsci(3); // Green color
cpgline(nintervals + 1, t, theta2);
// Calculate the min and max values for velocity
float min_omega = omega1[0], max_omega = omega1[0];
for (int i = 0; i <= nintervals; i++) {
if (omega1[i] < min_omega) min_omega = omega1[i];
if (omega2[i] < min_omega) min_omega = omega2[i];
if (omega1[i] > max_omega) max_omega = omega1[i];
if (omega2[i] > max_omega) max_omega = omega2[i];
}
// Add some padding for better visualization
float padding_omega = (max_omega - min_omega) * 0.1;
min_omega -= padding_omega;
max_omega += padding_omega;
// -------------------------- Graph 2: Velocity ---------------------------
cpgopen("/xs");
cpgenv(0.0, t[nintervals], min_omega, max_omega, 0, 1); // Time vs. Angular
Velocity
cpglab("Time (s)", "Omega (rad/s)", "Velocity of Pendulums");
// Velocity of Pendulum 1 (Red)
cpgsci(2); // Red color
cpgline(nintervals + 1, t, omega1);
// Velocity of Pendulum 2 (Green)
cpgsci(3); // Green color
cpgline(nintervals + 1, t, omega2);
}
int main() {
// Get double pendulum parameters and integrate
float L1, L2, m1, m2, b, theta1_0, theta2_0, omega1_0, omega2_0;
std::cout << "Enter length of first pendulum (m): "; std::cin >> L1;
std::cout << "Enter length of second pendulum (m): "; std::cin >> L2;
std::cout << "Enter mass of first pendulum (kg): "; std::cin >> m1;
std::cout << "Enter mass of second pendulum (kg): "; std::cin >> m2;
std::cout << "Enter damping coefficient: "; std::cin >> b;
std::cout << "Enter initial angle of first pendulum (rad): "; std::cin >>
theta1_0;
std::cout << "Enter initial angle of second pendulum (rad): "; std::cin >>
theta2_0;
std::cout << "Enter initial angular velocity of first pendulum (rad/s): ";
std::cin >> omega1_0;
std::cout << "Enter initial angular velocity of second pendulum (rad/s): ";
std::cin >> omega2_0;
graphs(L1, L2, m1, m2, b, theta1_0, theta2_0, omega1_0, omega2_0);
simulation(L1, L2, m1, m2, b, theta1_0, theta2_0, omega1_0, omega2_0);
return 0;
}
