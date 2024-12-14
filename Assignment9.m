clear
clc
close all

% Define the parameters
h = 50; % height in meters
L = 150; % length in meters
mu = 0.03; % coefficient of friction
g = 9.81; % acceleration due to gravity in m/s^2

% Define the cycloid parameters
R = h / (2 * pi); % radius of the generating circle
theta_max = 2 * pi; % maximum angle for the cycloid

% Define the cycloid equations
x_cycloid = @(theta) R * (theta - sin(theta));
y_cycloid = @(theta) R * (1 - cos(theta));

% Define the range of theta
theta = linspace(0, theta_max, 1000);

% Calculate the cycloid coordinates
x = x_cycloid(theta);
y = y_cycloid(theta);

% Scale the cycloid to fit the given length and height
x = x * (L / max(x));
y = y * (h / max(y));

% Plot the cycloid
figure;
plot(x, y, 'b', 'LineWidth', 2);
hold on;
xlabel('x (m)');
ylabel('y (m)');
title('Brachistochrone Curve with Friction');
grid on;

% Calculate the time taken to travel along the curve with friction
v = sqrt(2 * g * (h - y) ./ (1 + mu * (x ./ sqrt(x.^2 + y.^2))));
t = trapz(x, 1 ./ v);

% Display the total time
disp(['Total time to travel along the curve: ', num2str(t), ' seconds']);

% Mark the start and end points
plot(0, h, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(L, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

hold off;