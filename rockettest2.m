%--------------------------------------------------------------------------
% Falcon 9 Ascent (1D Vertical) with:
%    - Exponential Atmosphere (for density)
%    - Simplified Speed of Sound vs. Altitude
%    - Mach-Dependent Drag Coefficient
%    - Sutton-Graves Heat Flux (Optional)
%
%  NOTE: This is still a simplified demonstration, not a full flight sim.
%--------------------------------------------------------------------------
clear; clc; close all;

%% -- USER-PARAMS --
thrust   = 7.6e6;      % [N]   Constant thrust (approx for F9 first stage)
mdot     = 2200;       % [kg/s] Propellant mass flow rate
m0       = 550000;     % [kg]  Initial rocket mass
g        = 9.81;       % [m/s^2] Gravitational acceleration
dt       = 1.0;        % [s]   Time step
target_h = 65000;      % [m]   Target altitude to stop simulation (65 km)
useHeatFlux = true;    % Toggle heat flux calculation

%% -- ROCKET GEOMETRY --
rocket_d = 3.7;        % [m]  Approx. rocket body diameter
A        = pi*(rocket_d/2)^2; % Cross-sectional area [m^2]

%% -- ATMOSPHERE PARAMS --
rho0     = 1.225;      % [kg/m^3] Sea-level density
Hscale   = 8500;       % [m] Scale height for exponential model

%% -- HEAT FLUX PARAMS (Sutton-Graves style) --
%   Q_dot = k * sqrt(rho/Rn) * v^3.05  [W/m^2]
k_heat   = 1.83;       % Empirical constant
Rn       = 1.0;        % [m] Nose radius

%% -- INITIAL CONDITIONS --
t = 0.0;               % [s]
m = m0;                % [kg]
v = 0.0;               % [m/s]
h = 0.0;               % [m]

%% -- DATA STORAGE --
% Columns: [time, mass, altitude, velocity, acceleration, Mach, C_d, heatFlux]
simulation_data = [];

%--------------------------------------------------------------------------
% Main loop
%--------------------------------------------------------------------------
while (h < target_h) && (m > 0)
    % 1) Compute local density (exponential approximation)
    rho = rho0 * exp(-h / Hscale);

    % 2) Estimate local temperature (very simplified model)
    %    We'll do a trivial linear drop up to ~11 km, then clamp.
    %    Real standard atmosphere is more nuanced.
    if h < 11000
        T = 288.15 - 0.0065*h;  % Troposphere (ISA)
        if T < 216.65
            T = 216.65;        % Just clamp to ~stratosphere floor
        end
    else
        T = 216.65;            % Stratosphere approximation
    end
    
    % 3) Speed of sound
    gamma = 1.4;       % for air
    R_air = 287;       % [J/(kg*K)]
    aLoc  = sqrt(gamma * R_air * T);   % [m/s]

    % 4) Mach number
    if aLoc < 1e-6
        Mach = 0;  % avoid divide by zero in weird edge cases
    else
        Mach = v / aLoc;
    end
    if Mach < 0
        Mach = abs(Mach);  % If rocket is somehow going downward, still define Mach>0
    end
    
    % 5) Mach-dependent drag coefficient
    %    This is a TOTALLY made-up piecewise function for demonstration.
    %    In reality, you'd have a real C_d(M) curve from wind tunnel or CFD data.
    C_d = getDragCoefficient(Mach);
    
    % 6) Drag force
    %    F_drag = 0.5 * C_d * A * rho * v^2
    %    Direction: opposes motion (assuming upward flight for now).
    F_drag = 0.5 * C_d * A * rho * (v^2);
    
    % 7) Net force
    %    F_net = Thrust - Weight - Drag
    weight = m * g;
    F_net = thrust - weight - F_drag;
    
    % 8) Acceleration
    a = F_net / m;
    
    % 9) Update velocity & altitude
    v = v + a * dt;
    h = h + v * dt;
    
    % 10) Mass depletion
    m = m - mdot * dt;
    
    % 11) Heat flux (if desired)
    heatFlux = 0;
    if useHeatFlux && (v > 0)
        heatFlux = k_heat * sqrt(rho / Rn) * v^3.05; % [W/m^2]
    end
    
    % 12) Increment time
    t = t + dt;
    
    % 13) Store data every 5 seconds
    if abs(mod(t, 5)) < 1e-9
        simulation_data = [simulation_data; ...
                           t, m, h, v, a, Mach, C_d, heatFlux];
    end
end

%--------------------------------------------------------------------------
% Display a portion of the simulation results
%--------------------------------------------------------------------------
disp(' Time (s) |  Altitude (m) | Vel (m/s) |   Mach  |   Cd   | HeatFlux (W/m^2)');
for i = 1:4:size(simulation_data, 1)  % Print every 20s (4 steps if storing every 5s)
    fprintf('%8.1f | %12.1f | %9.1f | %7.2f | %6.3f | %14.2f\n', ...
        simulation_data(i,1), ...
        simulation_data(i,3), ...
        simulation_data(i,4), ...
        simulation_data(i,6), ...
        simulation_data(i,7), ...
        simulation_data(i,8));
end

%--------------------------------------------------------------------------
% Optional: Plot altitude and Mach number vs. time
%--------------------------------------------------------------------------
figure;
subplot(2,1,1);
plot(simulation_data(:,1), simulation_data(:,3)*1e-3, 'b-o','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Altitude (km)');
title('1D Ascent with Mach-Dependent Drag');
grid on;

subplot(2,1,2);
plot(simulation_data(:,1), simulation_data(:,6), 'r-o','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Mach Number');
grid on;

%--------------------------------------------------------------------------
% End of Main Script
%--------------------------------------------------------------------------


%% Local function: Mach-based Drag Coefficient (totally illustrative)
function Cd_out = getDragCoefficient(Mach)
    % A very rough piecewise approximation:
    %   - Subsonic (M < 0.8):  ~0.2
    %   - Transonic (0.8 < M < 1.2): peaks near 0.3
    %   - Supersonic (1.2 < M < 5):  gradually drops
    %   - Hypersonic (M > 5):  not relevant for typical LEO ascent, but just in case
    %
    % In reality, you'd have a detailed table or polynomial fit from data.

    if Mach < 0.8
        Cd_out = 0.20;
    elseif Mach < 1.2
        % Let it peak around Mach=1.0
        Cd_out = 0.30;
    elseif Mach < 5
        % Gradual drop in supersonic range
        % You could do a linear or polynomial from 0.30 down to maybe 0.25
        Cd_out = 0.25 + 0.05*(2.0 - Mach/5);  
        % This is arbitrary, just for shape
        if Cd_out < 0.20, Cd_out = 0.20; end
    else
        % Hypersonic regime (not typical for a LEO rocket's lower atmosphere, but included)
        Cd_out = 0.20;
    end
end
