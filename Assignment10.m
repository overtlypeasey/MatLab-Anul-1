N = 1000; 
dt_min = 0; 
dt_max = 0.01; 
p_min = -1e-8; 
p_max = 1e-8; 

radius = 50e-6;
density_particle = 2500; 
Cd = 0.69; 
g = 10; 
density_air = 1.4; 

dt = dt_min + (dt_max - dt_min) * rand(N, 1); 
px = p_min + (p_max - p_min) * rand(N, 1); 
py = p_min + (p_max - p_min) * rand(N, 1); 
pz = p_min + (p_max - p_min) * rand(N, 1); 

t = cumsum(dt); 

x0 = 0;
y0 = 0;
z0 = 10;
vx = 0;
vy = 0;
vz = 0;

volume = (4/3) * pi * radius^3;
mass = density_particle * volume;

x = zeros(N, 1);
y = zeros(N, 1);
z = zeros(N, 1);
x(1) = x0;
y(1) = y0;
z(1) = z0;

for i = 2:N
    vx = vx + px(i) / mass;
    vy = vy + py(i) / mass;
    vz = vz + pz(i) / mass;
    
    v = sqrt(vx^2 + vy^2 + vz^2);
    Fd = 0.5 * Cd * density_air * pi * radius^2 * v^2;
    
    ax = -Fd * vx / (mass * v);
    ay = -Fd * vy / (mass * v);
    az = -Fd * vz / (mass * v) - g;
    
    vx = vx + ax * dt(i);
    vy = vy + ay * dt(i);
    vz = vz + az * dt(i);
    
    x(i) = x(i-1) + vx * dt(i);
    y(i) = y(i-1) + vy * dt(i);
    z(i) = z(i-1) + vz * dt(i);
end

figure;
plot3(x, y, z, 'LineWidth', 1.5);
xlabel('X Position (m)', 'FontSize', 12);
ylabel('Y Position (m)', 'FontSize', 12);
zlabel('Z Position (m)', 'FontSize', 12);
title('Trajectory of the Particle', 'FontSize', 14);
grid on;
axis equal;