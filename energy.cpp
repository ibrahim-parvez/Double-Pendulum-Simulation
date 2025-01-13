#include <cmath>
#include <iostream>
#include "cpgplot.h"
#include <tuple>
#include <unistd.h> // For usleep
// Constants
const int nintervals = 3000; // Number of intervals for integration
const float g = 9.81; // Gravity (m/s^2)
const float dt = 0.05; // Time step for numerical integration (s)
// Function to integrate the double pendulum using the Runge-Kutta method
void integrateDoublePendulum(float L1, float L2, float m1, float m2, float
theta1_0, float theta2_0, float omega1_0, float omega2_0) {
float t[nintervals + 1];
float theta1[nintervals + 1], theta2[nintervals + 1];
float omega1[nintervals + 1], omega2[nintervals + 1];
float x1[nintervals + 1], y1[nintervals + 1], x2[nintervals + 1],
y2[nintervals + 1];
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
// Helper function to compute derivatives (motion equations)
auto calcDerivatives = [&](float theta1, float theta2, float omega1, float
omega2) {
float delta = theta1 - theta2;
// Common denominator
float denom = 2 * m1 + m2 - m2 * cos(2 * delta);
// Angular acceleration of theta1
float num1 = -g * (2 * m1 + m2) * sin(theta1)
- m2 * g * sin(theta1 - 2 * theta2)
- 2 * sin(delta) * m2 * (omega2 * omega2 * L2 + omega1 *
omega1 * L1 * cos(delta));
float alpha1 = num1 / (L1 * denom);
// Angular acceleration of theta2
float num2 = 2 * sin(delta) * (omega1 * omega1 * L1 * (m1 + m2)
+ g * (m1 + m2) * cos(theta1)
+ omega2 * omega2 * L2 * m2 * cos(delta));
float alpha2 = num2 / (L2 * denom);
return std::tuple<float, float, float, float>(omega1, alpha1, omega2,
alpha2);
};
// Open the plot window once
cpgopen("/XWINDOW");
cpgsci(1);
cpgenv(-L1 - L2 - 1, L1 + L2 + 1,
-L1 - L2 - 1, L1 + L2 + 1, 0, 1);
cpglab("X (meters)", "Y (meters)", "Double Pendulum Animation");
// Runge-Kutta integration loop
for (int i = 1; i <= nintervals; i++) {
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
// Clear the previous plot (using background color)
cpgsci(0); // Set color to black (background)
cpgrect(-L1 - L2 - 1, L1 + L2 + 1,
-L1 - L2 - 1, L1 + L2 + 1); // Clear
the screen
// Plot pendulum 1 (mass 1)
cpgsci(2); // Red for pendulum 1
cpgpt(1, &x1[i], &y1[i], 17); // Plot point for mass 1
// Plot pendulum 2 (mass 2)
cpgsci(3); // Green for pendulum 2
cpgpt(1, &x2[i], &y2[i], 17); // Plot point for mass 2
// Draw line from origin to pendulum 1
cpgsci(2); // Red line
float origin_x = 0.0, origin_y = 0.0; // Origin coordinates
float line1_x[2] = {origin_x, x1[i]};
float line1_y[2] = {origin_y, y1[i]};
cpgline(2, line1_x, line1_y);
// Draw line from pendulum 1 to pendulum 2
cpgsci(3); // Green line
float line2_x[2] = {x1[i], x2[i]};
float line2_y[2] = {y1[i], y2[i]};
cpgline(2, line2_x, line2_y);
// Wait for 50 ms before the next update (for smooth animation)
usleep(50000); // 50 ms
}
cpgclos(); // Close the plot after the loop finishes
}
int main() {
float L1, L2, m1, m2, theta1_0, theta2_0, omega1_0, omega2_0;
// User input for initial conditions
std::cout << "Enter length of first pendulum (m): "; std::cin >> L1;
std::cout << "Enter length of second pendulum (m): "; std::cin >> L2;
std::cout << "Enter mass of first pendulum (kg): "; std::cin >> m1;
std::cout << "Enter mass of second pendulum (kg): "; std::cin >> m2;
std::cout << "Enter initial angle of first pendulum (rad): "; std::cin >>
theta1_0;
std::cout << "Enter initial angle of second pendulum (rad): "; std::cin >>
theta2_0;
std::cout << "Enter initial angular velocity of first pendulum (rad/s): ";
std::cin >> omega1_0;
std::cout << "Enter initial angular velocity of second pendulum (rad/s): ";
std::cin >> omega2_0;
integrateDoublePendulum(L1, L2, m1, m2, theta1_0, theta2_0, omega1_0,
omega2_0);
return 0;
}
