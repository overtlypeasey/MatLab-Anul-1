% Fastest Trajectory Using Cycloid (Brachistochrone Problem)
% Parameters:
% - Start Point: (0, 50)
% - End Point: (150, 0)
% - Coefficient of friction (mu) = 0 (Frictionless)
% Goal: Find the cycloid parameters that define the fastest path

% Clear workspace and command window
clear; clc;

% Given parameters
H = 50;          % Vertical height difference in meters (y from 50 to 0)
L = 150;         % Horizontal length in meters
mu = 0;          % Coefficient of friction (set to 0 for cycloid)
g = 9.81;        % Acceleration due to gravity (m/s^2)

% Number of discretization points
N = 100;         % Number of segments

% Define the cycloid parametric equations
% x(theta) = r (theta - sin(theta))
% y(theta) = 50 - r (1 - cos(theta))

% Set up the system of equations to solve for r and theta_max
% We need to satisfy:
% r * (theta_max - sin(theta_max)) = L = 150
% r * (1 - cos(theta_max)) = H = 50

% Initial guesses for r and theta_max
initial_guess = [10, pi];  % [r, theta_max]

% Define the system of nonlinear equations
cycloid_eqns = @(vars) [
    vars(1) * (vars(2) - sin(vars(2))) - L;
    vars(1) * (1 - cos(vars(2))) - H
    ];

% Options for fsolve
options = optimset('Display','iter','TolX',1e-12,'TolFun',1e-12);

% Solve for r and theta_max using fsolve
[solution, fval, exitflag, output] = fsolve(cycloid_eqns, initial_guess, options);

% Extract optimal parameters
r_opt = solution(1);
theta_max = solution(2);

% Check if solution converged
if exitflag <= 0
    error('fsolve did not converge to a solution. Please check initial guesses or equations.');
end

% Display optimal parameters
fprintf('Optimal Cycloid Parameters:\n');
fprintf('Radius (r) = %.6f m\n', r_opt);
fprintf('Theta Max (theta_max) = %.6f radians\n', theta_max);

% Generate the cycloid path
theta = linspace(0, theta_max, N+1);
x = r_opt * (theta - sin(theta));
y = 50 - r_opt * (1 - cos(theta));

% Ensure the final point matches exactly (150, 0)
x(end) = L;
y(end) = 0;

% Plot the cycloid trajectory
figure;
plot(x, y, 'b-', 'LineWidth', 2);
hold on;
plot(0, 50, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');          % Start point (0,50)
plot(L, 0, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');          % End point (150,0)
xlabel('Horizontal Distance (m)');
ylabel('Vertical Distance (m)');
title('Optimal Cycloid Trajectory from (0,50) to (150,0)');
legend('Cycloid Path', 'Start (0,50)', 'End (150,0)', 'Location', 'Best');
grid on;
axis equal;
hold off;

% Function to calculate total descent time for cycloid
% Since the cycloid is the analytical solution, this is for verification
total_time = compute_cycloid_descent_time(r_opt, theta_max, g);

fprintf('Total Descent Time: %.6f seconds\n', total_time);

% Function to compute total descent time for cycloid
function T = compute_cycloid_descent_time(r, theta_max, g)
    % The descent time for a cycloid is known analytically:
    % T = sqrt(r / g) * theta_max
    T = sqrt(r / g) * theta_max;
end
