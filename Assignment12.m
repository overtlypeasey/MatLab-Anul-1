% filepath: /Users/victor/Documents/MATLAB/MatLab Anul 1/Assignment12.m
clc
clear
close all

% Court dimensions (ft)
court_length = 32;
court_width = 21;
court_height = 15;

% Ball properties
ball_diameter = 40/1000; % meters
ball_mass = 0.025; % kg
initial_speed = 100 * 0.3048; % m/s
initial_position = [6, 6, 5] * 0.3048; % meters
target_position = [court_length - 6, court_width/2, court_height/2] * 0.3048;

% Energy return coefficients
coeff_wall = 0.94;
coeff_floor = 0.90;

% Friction coefficients
mu_floor = 0.05;
mu_wall = 0.08;

% Define gravitational acceleration
gravity = [0, 0, -9.81]; % m/s²

% Calculate weight of the ball
weight = ball_mass * gravity; % N

% Simulation parameters
dt = 0.01; % time step (s)
total_time = 10; % total simulation time (s)
num_steps = total_time / dt;

% Initial velocity towards target
direction = (target_position - initial_position);
direction = direction / norm(direction);
velocity = direction * initial_speed;

% Initialize position history
position_history = zeros(num_steps, 3);
position = initial_position;

for i = 1:num_steps
    % Update velocity with gravity
    velocity = velocity + weight * dt;
    
    % Update position
    position = position + velocity * dt;
    
    % Check for collisions with walls
    if position(1) <= 0 || position(1) >= court_length
        velocity(1) = -velocity(1) * coeff_wall;
        velocity(2) = velocity(2) * (1 - mu_wall);
    end
    if position(2) <= 0 || position(2) >= court_width
        velocity(2) = -velocity(2) * coeff_wall;
        velocity(1) = velocity(1) * (1 - mu_wall);
    end
    % Check for collisions with floor
    if position(3) <= 0
        velocity(3) = -velocity(3) * coeff_floor;
        velocity(1:2) = velocity(1:2) * (1 - mu_floor);
    end
    % Check for collisions with ceiling
    if position(3) >= court_height
        velocity(3) = 0;
    end
    
    % Record position
    position_history(i, :) = position;
    
    % Terminate if velocity is negligible
    if norm(velocity) < 0.1
        position_history = position_history(1:i, :);
        break
    end
end

% Plot trajectory
figure;
plot3(position_history(:,1), position_history(:,2), position_history(:,3), 'b');
hold on;

% Plot start point
plot3(initial_position(1), initial_position(2), initial_position(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Plot court boundaries
% Define court corners
corners = [
    0, 0, 0;
    court_length, 0, 0;
    court_length, court_width, 0;
    0, court_width, 0;
    0, 0, 0;
    0, 0, court_height;
    court_length, 0, court_height;
    court_length, court_width, court_height;
    0, court_width, court_height;
    0, 0, court_height
];

% Plot court edges
plot3(corners(:,1), corners(:,2), corners(:,3), 'k--');

% Labels and title
xlabel('Length (m)');
ylabel('Width (m)');
zlabel('Height (m)');
title('Ball Trajectory in Squash Court');
legend('Trajectory', 'Start Point', 'Court Outline');
grid on;
hold off;