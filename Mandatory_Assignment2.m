clear
clc
close all

x = -10:0.1:10;

y = x.^2 + 4*x - 3;

plot(x, y, 'k', 'LineWidth', 2);
hold on;

theta = 33;

theta_rad = deg2rad(theta);

x_rot = x * cos(theta_rad) - y * sin(theta_rad);
y_rot = x * sin(theta_rad) + y * cos(theta_rad);

plot(x_rot, y_rot, 'g', 'LineWidth', 2);

[x_intersect, y_intersect] = polyxpoly(x, y, x_rot, y_rot);

plot(x_intersect, y_intersect, 'r*', 'MarkerSize', 10);

xlabel('x');
ylabel('y');
title('Plot of y = x^2 + 4x - 3 and its 33-degree rotation with intersections');

for i = 1:length(x_intersect)-1
    patch_x = [x_intersect(i), x_intersect(i+1), x_rot(x_rot >= x_intersect(i) & x_rot <= x_intersect(i+1)), x(x >= x_intersect(i) & x <= x_intersect(i+1))];
    patch_y = [y_intersect(i), y_intersect(i+1), y_rot(x_rot >= x_intersect(i) & x_rot <= x_intersect(i+1)), y(x >= x_intersect(i) & x <= x_intersect(i+1))];
    
    fill(patch_x, patch_y, 'y', 'FaceAlpha', 0.5);
end

grid on;

hold off;