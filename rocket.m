clear
clc
close all

m = 549054; %kg
mb = 22000;
thrust = (7600-6000).*1000; %kn.*10000
g = 9.655;
rho = 0.3227;
ff = 2200; %kg/s
G = m.*g;
CD = 0.85;
theta = 0;
orbital_inclination_target = 90;
inclination_rate = 0.488;

dt = 1;
total_time = 200;
t = 0:dt:total_time;
v = zeros(1, total_time); %m/s
vy = 0;
drag_force = 1/2.*rho.*v(1)^2.*CD;
altitude = zeros(1, total_time); %km

for i = 1:length(t)
    if thrust < 7600*1000
        thrust = thrust+1000*1000;
    end
    a = (thrust - G - drag_force)/m;
    if a<0
        a = 0;
    end
    tana = a .* cos(theta);
    if i > 1
        v(i) = v(i-1)+a;
    else
        v(i) = a;
    end
    vy = vy + tana;
    m = m-ff;
    % display(a);
    % display(thrust);
    % display(drag_force)
    G = m.*g;
    drag_force = 1/2.*rho.*v(i)^2.*CD;
    if i > 1
        altitude(i) = altitude(i-1) + vy;
    else
        altitude(i) = vy;
    end
    if theta < orbital_inclination_target
        theta = theta+inclination_rate;
    end
    if i == 162
        m = m-mb;
        thrust = 1000^2;
    end
end
display(v(i).*3.6)
display(altitude(i)./1000)

figure
plot(t, altitude);
xlabel("seconds")
ylabel("m");

figure
plot(t, v);
xlabel("seconds")
ylabel("m/s");

