clear
clc
close all

N = 50;

x = randi([1, 50], 1, N);
y = randi([1, 50], 1, N);
z = randi([1, 50], 1, N);

center_of_mass_x = mean(x);
center_of_mass_y = mean(y);
center_of_mass_z = mean(z);

x_translated = x - center_of_mass_x;
y_translated = y - center_of_mass_y;
z_translated = z - center_of_mass_z;
center_translated_x = mean(x_translated);
center_translated_y = mean(y_translated);
center_translated_z = mean(z_translated);

r = sqrt(x.^2+y.^2);
theta = atan(y./x);

plot3(x_translated, y_translated, z_translated);
hold on
plot3(center_translated_x, center_translated_y, center_translated_z, 'gd')
text(center_translated_x, center_translated_y, center_translated_z, "CG")
hold off

polar_coords = [r; theta; z]';

sorted_polar_coords = sortrows(polar_coords, 2);

r_sorted = sorted_polar_coords(:, 1);
theta_sorted = sorted_polar_coords(:, 2);
z_sorted = sorted_polar_coords(:, 3);