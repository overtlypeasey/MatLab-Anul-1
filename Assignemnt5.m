clc
clear
close all

x = linspace(-2, 2, 600);
y = linspace(-2, 2, 600);
[X, Y] = meshgrid(x, y);
C = X + i * Y;

% Initialize the iteration count matrix
maxIter = 50;
Z = zeros(size(C));
M = zeros(size(C));

% Mandelbrot iteration
for k = 1:maxIter
    Z = Z.^2 + C;
    M(abs(Z) <= 1000) = k;
end

% Plot the Mandelbrot set
figure
image(x, y, M)
colormap([jet(); flipud(jet()); 0 0 0]);
colorbar
axis off;
title('Mandelbrot Set');

