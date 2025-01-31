clear
clc
close all

N = 10;
x = 20 * rand(1, N) - 10;
y = 20 * rand(1, N) - 10;

% Compute the centroid of the points
centroid_x = mean(x);
centroid_y = mean(y);

% Calculate angles from the centroid to each point
angles = atan2(y - centroid_y, x - centroid_x);

% Sort the points based on the calculated angles
[~, sorted_idx] = sort(angles);
x = x(sorted_idx);
y = y(sorted_idx);

% Close the path by returning to the first point
x = [x, x(1)];
y = [y, y(1)];

figure;
plot(x, y, '-o');
xlabel('X-axis');
ylabel('Y-axis');
title('Closed Path');
grid on;