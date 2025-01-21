clear
clc
close all

% Rocket parameters
m = 549054; % kg
mb = 22000; % booster mass
thrust = (7600-6000).*1000; % kn.*10000
thrust0 = thrust;
ff = 2200; % fuel flow in kg/s
jettison_point = 162; % in seconds
Area = 10; % cross-section of the rocket

% Earth parameters
g0 = 9.81; % Earth surface gravity
rho = 1.225; % density MSL
G = m.*g0; % Weight formula
Re = 6371000; % m (Earth radius)
CD = 0.85; % coeff of drag
H = 8500; % m (scale height for atmosphere)
G_const = 6.674e-11; % Gravitational constant
M_earth = 5.972e24; % Mass of Earth in kg

% Gravity turn
phi = 0; % Initial angle for gravity turn
orbital_inclination_target = 90;
inclination_rate = 0.3;

% Simulation parameters & initial state
dt = 1; % time increment in s
total_time = 400;
t = 0:dt:total_time;
num_steps = length(t); % Number of simulation steps
v = zeros(1, length(t)); % m/s
vr = zeros(1, length(t)); % radial velocity
vtheta = zeros(1, length(t)); % tangential velocity
r = Re * ones(1, length(t)); % radius from Earth's center
altitude = zeros(1, length(t)); % altitude from Earth's surface
thermal_flux = zeros(1, num_steps); % Thermal flux (W/m²)
orbital_insertion = false;

theta = zeros(1, length(t)); % Initialize theta for storing angular position
drag_force = 1/2.*rho.*v(1)^2.*CD.*Area;

for i = 1:length(t)
    % Thrust spool up
    if thrust < 7600*1000
        thrust = thrust + 1000*1000;
    end

    % Decompose thrust into radial and tangential components
    thrust_r = thrust * cosd(phi);
    thrust_theta = thrust * sind(phi);

    % Radial and tangential accelerations
    ar = (thrust_r - G - drag_force) / m;
    atheta = thrust_theta / m;

    % Use absolute value of radial acceleration if altitude is below a threshold (e.g., 10 meters)
    if altitude(i) < 10
        ar = abs(ar);
    end

    % Update velocities
    if i > 1
        vr(i) = vr(i-1) + ar * dt;
        vtheta(i) = vtheta(i-1) + atheta * dt;
    else
        vr(i) = ar * dt;
        vtheta(i) = atheta * dt;
    end

    % Calculate the resultant velocity
    v(i) = sqrt(vr(i)^2 + vtheta(i)^2);

    % Update positions
    if i > 1
        r(i) = r(i-1) + vr(i) * dt;
        theta(i) = theta(i-1) + (vtheta(i) / r(i)) * dt;
    else
        r(i) = vr(i) * dt;
        theta(i) = (vtheta(i) / r(i)) * dt;
    end
    altitude(i) = r(i) - Re;

    % Opposing forces update
    G = m * g0 * (Re / r(i))^2; % Update gravitational force with altitude
    drag_force = 1/2 * rho * v(i)^2 * CD * Area;

    % Gravity turn evolution
    if phi < orbital_inclination_target && i > 10
        phi = phi + inclination_rate * dt;
    end

    % Jettisoning
    if i == jettison_point
        m = m - mb;
        thrust0 = 981 * 1000;
        thrust = 981 * 1000;
        ff = ff / 9;
    end

    % Thermal Flux Calculation
    % Using the empirical formula: q = 0.5 * rho * v^3 * CD
    thermal_flux(i) = 0.5 * rho * v(i)^3 * CD; % Thermal flux (W/m²)

    % Update atmospheric density
    if altitude(i) < 11000
        rho = 1.225 * (1 - 0.0000226 * altitude(i))^4.255;
    else
        rho = rho * exp(-altitude(i) / H);
    end

    % Orbital velocity calculation
    v_orbital_required = sqrt(G_const * M_earth / r(i)); % Orbital velocity
    
    % Check orbital insertion
    if altitude(i) > 200e3 && v(i) >= v_orbital_required % Minimum 200 km altitude
        orbital_insertion = true;
        fprintf('Orbital insertion achieved at t = %d seconds\n', t(i));
        fprintf('Altitude: %.2f km\n', altitude(i) / 1000);
        fprintf('Velocity: %.2f m/s\n', v(i));
        thrust = 0;
        break; % Stop simulation after orbital insertion
    end
end

disp(m)

figure
plot(t, altitude);
xlabel("Time (seconds)");
ylabel("Altitude (meters)");
hold on 
scatter(jettison_point, altitude(jettison_point), "filled", "d")
hold off

figure
plot(t, v);
xlabel("Time (seconds)");
ylabel("Velocity (m/s)");
hold on 
scatter(jettison_point, v(jettison_point), "filled", "d")
hold off

% Plot Thermal Flux vs. Time
figure
plot(t, thermal_flux, 'k', 'LineWidth', 2);
xlabel('Time (seconds)');
ylabel('Thermal Flux (W/m²)');
title('Rocket Thermal Flux Over Time');
grid on;

% Convert polar coordinates to Cartesian for plotting
x = r .* cos(theta);
y = r .* sin(theta);

% Plot the trajectory of the rocket
figure;
plot(x, y, 'b', 'LineWidth', 2);
hold on;
plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Earth
xlabel('X Position (meters)');
ylabel('Y Position (meters)');
title('Rocket Trajectory');
axis equal;
grid on;

% Plot the curvature of the Earth for reference
theta_earth = linspace(0, 2*pi, 100);
earth_x = Re * cos(theta_earth);
earth_y = Re * sin(theta_earth);
plot(earth_x, earth_y, 'g--');
axis equal;
hold off;