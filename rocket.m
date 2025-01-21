clear
clc
close all

%rocket parameters
m = 549054; %kg
mb = 22000; %booster mass
thrust = (7600-6000).*1000; %kn.*10000
thrust0 = thrust;
ff = 2200; %fuel flow in kg/s
jettison_point = 162; % in seconds
Area = 10; %cross-section of the rocket

%Earth parameters
g = 9.81; %Earth accel
rho = 1.225; %density MSL
G = m.*g; % Weight formula
Re = 6371000; % m (Earth radius)
CD = 0.85; %coeff of drag
H = 8500; % m (scale height for atmosphere)
G_const = 6.674e-11; % Gravitational constant
M_earth = 5.972e24; % Mass of Earth in kg


%Gravity turn
theta = 0;
orbital_inclination_target = 90;
inclination_rate = 0.3;

%simulation parameters & initial state
dt = 1; %time increment in s
total_time = 400;
t = 0:dt:total_time;
num_steps = length(t); % Number of simulation steps
v = zeros(1, length(t)); %m/s
vy = 0; % vertical speed
drag_force = 1/2.*rho.*v(1)^2.*CD.*Area;
altitude = zeros(1, length(t)); % meters
thermal_flux = zeros(1, num_steps); % Thermal flux (W/m²)


for i = 1:length(t)
    %Thrust spool up
    if thrust < 7600*1000
        thrust = thrust+1000*1000;
    end

    %accel calculation
    a = (thrust - G - drag_force)/m;
    if a<0
        a = 0;
    end
    tana = a .* cosd(theta); % vertical acceleration
    
    %relative velocity increase
    if i > 1
        v(i) = v(i-1)+a;
    else
        v(i) = a;
    end

    %vertical velocity increase
    % vy = vy + tana;
    vy = v(i) .* cosd(theta);

    %mass depletion
    m = m-ff;

    % display(a);
    % display(thrust);
    % display(drag_force)

    %opposing forces update
    G = m.*g;
    drag_force = 1/2.*rho.*v(i)^2.*CD.*Area;

    %altitude increase
    if i > 1
        altitude(i) = altitude(i-1) + vy;
    else
        altitude(i) = vy;
    end
    h = altitude(i);

    %Gravity turn evolution
    if theta < orbital_inclination_target && i > 10
        theta = theta+inclination_rate;
    end

    %jettisoning
    if i == jettison_point
        m = m-mb;
        thrust0 = 981.*1000;
        thrust = 981.*1000;
        ff = ff/9;
    end


    % **Thermal Flux Calculation**
    % Using the empirical formula: q = 0.5 * rho * v^3 * CD
    thermal_flux(i) = 0.5 * rho * v(i)^3 * CD; % Thermal flux (W/m²)

    % Update atmospheric density and gravity
    h = altitude(i);
    if altitude < 11000
        rho = 1.225 * (1 - 0.0000226*h).^4.255;
    else
        rho = rho * exp(-h / H);
    end
    g = g * (Re / (Re + h))^2;

       % Orbital velocity calculation
    r = Re + altitude(i); % Current radius (Earth + altitude)
    v_orbital_required = sqrt(G_const * M_earth / r); % Orbital velocity
    
    %Check orbital insertion
    if altitude(i) > 200e3 && v(i) >= v_orbital_required % Minimum 200 km altitude
        orbital_insertion = true;
        fprintf('Orbital insertion achieved at t = %d seconds\n', t(i));
        fprintf('Altitude: %.2f km\n', altitude(i)/1000);
        fprintf('Velocity: %.2f m/s\n', v(i));
        thrust = 0;
    else
        fprintf('Failure after orbital insertion');
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