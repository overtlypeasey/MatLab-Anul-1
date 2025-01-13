clear
clc
close all

%rocket parameters
m = 549054; %kg
mb = 22000; %booster mass
thrust = (7600-6000).*1000; %kn.*10000
ff = 2200; %fuel flow in kg/s
jettison_point = 162; % in seconds
Area = 10; %cross-section of the rocket

%Earth parameters
g = 9.81; %Earth accel
rho = 1.225; %density MSL
G = m.*g; % Weight formula
CD = 0.85; %coeff of drag

%Gravity turn
theta = 0;
orbital_inclination_target = 90;
inclination_rate = 0.488;

%simulation parameters & initial state
dt = 1; %time increment in s
total_time = 200;
t = 0:dt:total_time;
v = zeros(1, length(t)); %m/s
vy = 0; % vertical speed
drag_force = 1/2.*rho.*v(1)^2.*CD.*Area;
altitude = zeros(1, length(t)); % meters


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
    vy = vy + tana;

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

    %Gravity turn evolution
    if theta < orbital_inclination_target
        theta = theta+inclination_rate;
    end

    %jettisoning
    if i == jettison_point
        m = m-mb;
        thrust = 1000^2;
    end
end
display(v(i).*3.6)
display(altitude(i)./1000)

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

