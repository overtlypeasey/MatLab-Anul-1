% Two-Stage Rocket Ascent Simulation with Gravity Turn and Heat Flux
% Simulating Falcon 9's Ascent to Low Earth Orbit (LEO)
% Author: [Your Name]
% Date: [Today's Date]

clear; clc; close all;

%% 1. Constants and Parameters

% Earth Parameters
Re = 6371000; % Earth's radius (m)
mu = 3.986004418e14; % Earth's gravitational parameter (m^3/s^2)
g0 = 9.80665; % Sea-level gravitational acceleration (m/s^2)

% Atmospheric Model Parameters
rho0 = 1.225; % Sea-level atmospheric density (kg/m^3)
H_atm = 8500; % Scale height of atmosphere (m)

% Rocket Parameters
% Stage 1 Parameters
stage1.thrust = 7.607e6; % Thrust (N) [Equivalent to Falcon 9 first stage]
stage1.Isp = 282; % Specific Impulse (s) [Typical Merlin 1D sea level]
stage1.mass_dry = 50000; % Dry mass (kg)
stage1.mass_fuel = 830000; % Fuel mass (kg) [Adjusted for sufficient delta-V]
stage1.burn_time = 162; % Burn time (s)

% Stage 2 Parameters
stage2.thrust = 1.000e6; % Thrust (N) [Approximate Merlin Vacuum]
stage2.Isp = 348; % Specific Impulse (s) [Typical Merlin 1D vacuum]
stage2.mass_dry = 12000; % Dry mass (kg)
stage2.mass_fuel = 219420; % Fuel mass (kg) [Adjusted for sufficient delta-V]
stage2.burn_time = 397; % Burn time (s)

% Combined Parameters
stages = {stage1, stage2}; % Cell array to hold stages

% Rocket Cross-Sectional Area and Drag Coefficient
A = 10; % Cross-sectional area (m^2) [Adjust as per rocket design]
Cd = 0.5; % Drag coefficient (dimensionless)

% Gravity Turn Parameters
% Initial Flight Path Angle (degrees)
initial_theta = 90; % Vertical launch

% Final Flight Path Angle (degrees)
final_theta = 0; % Horizontal launch for orbital insertion

% Time over which to perform the gravity turn (seconds)
gravity_turn_time = 200; % Duration of gravity turn

%% 2. Initial Conditions

% Initial State Vector: [x, y, vx, vy, mass]
% x: Horizontal position (m)
% y: Vertical position (m)
% vx: Horizontal velocity (m/s)
% vy: Vertical velocity (m/s)
% mass: Current mass (kg)

initial_mass = stages{1}.mass_dry + stages{1}.mass_fuel + stages{2}.mass_dry + stages{2}.mass_fuel;
initial_state = [0; 0; 0; 0; initial_mass]; % Starting at rest on the ground

%% 3. Simulation Time

t_start = 0;
t_end = stages{1}.burn_time + stages{2}.burn_time + 600; % Extended extra time after burns for orbital stability
t_span = [t_start t_end];

%% 4. Define the ODE System

% Define the ODE function with stage handling and gravity turn
function dstate_dt = rocket_ascent_ode(t, state, stages, A, Cd, rho0, H_atm, mu, Re, g0, initial_theta, final_theta, gravity_turn_time)
    % Unpack state
    x = state(1); % Horizontal position (m)
    y = state(2); % Vertical position (m)
    vx = state(3); % Horizontal velocity (m/s)
    vy = state(4); % Vertical velocity (m/s)
    mass = state(5); % Current mass (kg)
    
    % Calculate radial distance from Earth's center
    r = Re + y;
    
    % Determine current stage based on time
    if t <= stages{1}.burn_time
        current_stage = stages{1};
        stage_number = 1;
    elseif t <= (stages{1}.burn_time + stages{2}.burn_time)
        current_stage = stages{2};
        stage_number = 2;
    else
        current_stage = [];
        stage_number = 0;
    end
    
    % Calculate flight path angle (theta) for gravity turn
    if t <= gravity_turn_time
        theta = initial_theta - (initial_theta - final_theta) * (t / gravity_turn_time);
    else
        theta = final_theta;
    end
    theta_rad = deg2rad(theta); % Convert to radians
    
    % Calculate Thrust Components
    if ~isempty(current_stage)
        Fthrust = current_stage.thrust;
        Isp = current_stage.Isp;
        mdot = Fthrust / (Isp * g0); % Mass flow rate (kg/s)
        
        % Ensure we do not burn more fuel than available
        if mass <= (stages{1}.mass_dry + stages{2}.mass_dry)
            Fthrust = 0;
            mdot = 0;
        end
    else
        Fthrust = 0;
        mdot = 0;
    end
    
    % Update mass
    dm_dt = -mdot;
    if mass <= (stages{1}.mass_dry + stages{2}.mass_dry)
        dm_dt = 0;
    end
    
    % Atmospheric Density
    if y < 0
        rho = rho0; % Below ground level
    else
        rho = rho0 * exp(-y / H_atm);
    end
    
    % Calculate Velocity Magnitude
    v = sqrt(vx^2 + vy^2);
    
    % Calculate Drag Force
    Fd = 0.5 * rho * v^2 * Cd * A;
    
    % Components of Drag Force
    if v > 0
        Fdx = Fd * (vx / v);
        Fdy = Fd * (vy / v);
    else
        Fdx = 0;
        Fdy = 0;
    end
    
    % Gravitational Acceleration Components
    g = mu / (r)^2; % Gravitational acceleration magnitude (m/s^2)
    g_x = -mu * x / (r^3); % Gravitational acceleration component in x-direction
    g_y = -mu * y / (r^3); % Gravitational acceleration component in y-direction
    
    % Net Forces
    Fnet_x = Fthrust * cos(theta_rad) - Fdx + mass * g_x;
    Fnet_y = Fthrust * sin(theta_rad) - Fdy + mass * g_y;
    
    % Accelerations
    ax = Fnet_x / mass;
    ay = Fnet_y / mass;
    
    % Assign to output
    dx_dt = vx;
    dy_dt = vy;
    dvx_dt = ax;
    dvy_dt = ay;
    dmass_dt = dm_dt;
    
    dstate_dt = [dx_dt; dy_dt; dvx_dt; dvy_dt; dmass_dt];
end

%% 5. Run Simulation Using ode45

% Define the ODE function handle
ode_fun = @(t, y) rocket_ascent_ode(t, y, stages, A, Cd, rho0, H_atm, mu, Re, g0, initial_theta, final_theta, gravity_turn_time);

% Run the ODE solver with tighter tolerances for accuracy
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t_sol, state_sol] = ode45(ode_fun, t_span, initial_state, options);

%% 6. Calculate Heat Flux Post-Simulation

% Extract Trajectory Data
x = state_sol(:,1); % Horizontal position (m)
y = state_sol(:,2); % Vertical position (m)
vx = state_sol(:,3); % Horizontal velocity (m/s)
vy = state_sol(:,4); % Vertical velocity (m/s)
mass = state_sol(:,5); % Mass (kg)

% Calculate Velocity Magnitude
v_mag = sqrt(vx.^2 + vy.^2);

% Calculate Heat Flux
k = 1e-4; % Proportionality constant for heat flux
rho = rho0 * exp(-y / H_atm); % Atmospheric density at each altitude
heat_flux = k .* rho .* v_mag.^3; % Heat flux (W/m^2)

%% 7. Determine Maximum Altitude and Velocity

max_altitude = max(y);
max_velocity = max(v_mag);

% LEO Criteria
LEO_altitude_min = 160000; % 160 km
LEO_velocity_required = 7800; % 7.8 km/s

% Check if LEO Criteria are Met
if (max_altitude >= LEO_altitude_min) && (max_velocity >= LEO_velocity_required)
    LEO_status = 'LEO Achievement Criteria Met!';
else
    LEO_status = 'LEO Achievement Criteria Not Met.';
end

% Display Results
fprintf('Maximum Altitude Achieved: %.2f km\n', max_altitude / 1000);
fprintf('Maximum Velocity Achieved: %.2f km/s\n', max_velocity / 1000);
disp(LEO_status);

%% 8. Visualization

% Ensure that y and heat_flux are non-empty and have valid values
if isempty(y) || isempty(heat_flux) || isempty(v_mag)
    error('Simulation data arrays are empty. Please check the ODE integration.');
end

% Ensure that max(y), max(v_mag), and max(heat_flux) are finite and positive
if ~isfinite(max(y)) || max(y) <= 0
    error('Invalid maximum altitude value. Ensure that altitude data is correct.');
end

if ~isfinite(max(v_mag)) || max(v_mag) <= 0
    error('Invalid maximum velocity value. Ensure that velocity data is correct.');
end

if ~isfinite(max(heat_flux)) || max(heat_flux) < 0
    error('Invalid maximum heat flux value. Ensure that heat flux data is correct.');
end

% Plot Altitude vs. Time
figure;
subplot(3,1,1);
plot(t_sol, y / 1000, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Altitude (km)');
title('Rocket Altitude Over Time');
grid on;
ylim([0, max(y)/1000 * 1.1]);

% Plot Velocity vs. Time
subplot(3,1,2);
plot(t_sol, v_mag / 1000, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Velocity (km/s)');
title('Rocket Velocity Over Time');
grid on;
ylim([0, max(v_mag)/1000 * 1.1]);

% Plot Heat Flux vs. Time
subplot(3,1,3);
plot(t_sol, heat_flux, 'm', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Heat Flux (W/m^2)');
title('Heat Flux Over Time');
grid on;
ylim([0, max(heat_flux)*1.1]);
